include("nearbottom.jl")

abstract type AbstractSurveyDomain end

"""
    TransectRibbons([; width=20, buffer=0.1])

Specify how to define a survey domain based on equal-width strips centered on each 
transect. This method corresponds to MACE's usual calculations based on equally-spaced
transects.

# Arguments

- `transect_width` : The nominal spacing between transects, in nautical miles.
- `transect_buffer` : Amount to expand each transect strip side-to-side to ensure they
overlap.  Default is 0.1 (i.e., 10%).
"""
struct TransectRibbons{T} <:AbstractSurveyDomain
    transect_width::T
    buffer::T
end
TransectRibbons(; width=20, buffer=0.1) = TransectRibbons(promote(width, buffer)...)

"""
    SurveyHull([k=10])

Specify how to define a survey domain based on a concave hull that wraps around all the 
acoustic transects. The smoothness of this hull can be adjusted by setting the number of 
neighbors `k`.

"""
struct SurveyHull{T<:Integer} <: AbstractSurveyDomain
    k::T
end
SurveyHull(k=10) = SurveyHull(k)


function transect_ribbon(transect, transect_width, dx, buffer=0.1, ord=:y)
    tr1 = @chain transect begin
        @select(:x, :y, :log)
        stack(Not(ord))
        DataFramesMeta.@transform($ord = round.($ord ./ dx) .* dx)
        @by([ord, :variable], :value = mean(:value))
        unstack()
        @orderby($ord)
        @select(:x, :y)
        unique()
    end

    a = pi/2
    R1 = [cos(a) -sin(a); sin(a) cos(a)]
    R2 = [cos(-a) -sin(-a); sin(-a) cos(-a)]

    w = transect_width / 2 * 1.852 * (1 + buffer)

    X = Array(tr1)'
    v = diff(X, dims=2)
    v = v ./ norm.(eachcol(v))' .* w
    v = [v v[:, end]]
    left = R1 * v .+ X
    right = R2 * v .+ X

    ribbon_bounds = [tuple(x...) for x in [eachcol(left); eachcol(reverse(right, dims=2))]] 
    ribbon_bounds = [ribbon_bounds; (left[:, 1]...,)]
    return PolyArea(ribbon_bounds)
end

function survey_domain(acoustics, method::TransectRibbons, order, dx, dy=dx)
    tr_set = map(unique(acoustics.transect)) do i
        tr = transect_ribbon(@subset(acoustics, :transect .== i),
            method.transect_width, dx, method.buffer, order)
    end
    return GeometrySet(tr_set)
end

function survey_domain(acoustics, method::SurveyHull, order, dx, dy=dx)
    unique_points = @chain acoustics begin
        @select(:x, :y)
        unique()
    end
    v = [[row.x, row.y] for row in eachrow(unique_points)]
    hull = concave_hull(v, method.k)
    return Ngon([Point(x...) for x in hull.vertices]...)
end

"""
    grid_domain(domain, dx[, dy=dx])

Given a `GeometrySet` or `Domain` object, fill it with a dense rectangular grid with 
resolution `dx` and `dy`. The centers of these grid cells are returned in a `DataFrame`.
"""
function grid_domain(domain, dx, dy=dx)
    box = boundingbox(domain)
    xgrid = range(box.min.coords.x.val, box.max.coords.x.val, step=dx)
    ygrid = range(box.min.coords.y.val, box.max.coords.y.val, step=dx)
    grid = DataFrame([(;x, y) for x in xgrid, y in ygrid
        if in(Point(x, y), domain)])
    return grid
end

"""
    get_survey_grid(acoustics[; method=TransectRibbons(), dx=10.0, dy=dx, order=:y]])

Construct a regular grid inside the survey area, defined as the set of ribbon-like regions
with width `transect_width` along each survey transect.

# Arguments
- `acoustics` : `DataFrame` of georeferenced acoustic data. Should have `:x`, `:y`, and 
`:log` columns.
- `method` : How to define the survey domain. Options are `TransectRibbons()` (default) or 
`SurveyHull()`.
- `dx`, `dy` : Grid resolution, in km. Default is 10 km.
- `order` : `Symbol` indicating which column of `acoustics` to use to define the direction
of each transect for the purposes of defining the ribbon's boundaries. Defaults to `:y`,
since most of MACE's surveys have north-south transects. For east-west transects, use `:x`,
and for curving transects use `:log`.

# Returns
A tuple `(grid, domain)`, where `grid` is a `DataFrame` contianing the x and y coordinates
of each grid cell, and `domain` is a `GeometrySet` object that contains the geographic 
boundaries of the survey domain (either as a single polygon hull or a collection of 
transect ribbons).

"""
function get_survey_grid(acoustics; method=TransectRibbons(), dx=10.0, dy=dx, order=:y)
    nrow(acoustics) > 3 || error("Not enough locations to define survey grid.")
    domain = survey_domain(acoustics, method, order, dx, dy)
    grid = grid_domain(domain, dx, dy)
    return grid, domain
end


"""
Merge MACE's "macebase2.scaling_key_source_data" table with GAP's "racebase.specimen"
table to make a combined virtual SKSD table. This assigns all specimens from the GAP 
trawls to a new scaling stratum called "BT", which is applied to the bottom 3 m of the 
water column. GAP haul numbers are multiplied by -1 to make them unique from the MACE 
event_id's.
"""
function merge_scaling(scaling_mace, scaling_gap)
    ts_key = @by(scaling_mace, :species_code, :ts_relationship=first(:ts_relationship))
    
    scaling_mace1 = @select(scaling_mace, :survey, :ship, :event_id, :class, :species_code, 
        :primary_length, :ts_length, :ts_relationship, :catch_sampling_expansion,
        :user_defined_expansion, :sample_correction_scalar, :haul_weight, :w)

    scaling_gap1 = @chain scaling_gap begin
        DataFramesMeta.@transform(
            :survey = :cruise,
            :ship = :vessel,
            :event_id = -:haul,
            :class = "BT",
            :primary_length = :length ./ 10,
            :ts_length = :length ./ 10, # not exactly right
            :catch_sampling_expansion = 1,
            :sample_correction_scalar = 1,
            :haul_weight = 1,
            :w = 1
        )
        leftjoin(nearbottom_coefs, on=:species_code)
        leftjoin(ts_key, on=:species_code)
        DataFramesMeta.@transform(
            :ts_relationship = replace(:ts_relationship, missing => "generic_swimbladder_fish")
        )
        @select(:survey, :ship, :event_id, :class, :species_code,
            :primary_length, :ts_length, :ts_relationship, :catch_sampling_expansion,
            :user_defined_expansion, :sample_correction_scalar, :haul_weight, :w)
        dropmissing()
    end
    return [scaling_mace1; scaling_gap1]
end


function merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
    survey = only(unique(trawl_locations_mace.survey))
    trawl_locations_gap.survey .= survey
    return [trawl_locations_mace; trawl_locations_gap]
end

"""
preprocess_survey_data(surveydir[; ebs=true, log_ranges=nothing,
    grid_method=TransectRibbons(), dx=10.0, dy=dx, transect_order=:y,
    missingstring=[".", "NA"]])

Preprocess the survey data files in directory `surveydir` and save the outputs in 
the same directory as .csv files in standard format. 

# Arguments

## Required arguments

- `surveydir` : Path to the directory where the survey data files are located. This
function expects to find the following files, all of which come from running
"download_survey.R":

- scaling_mace.csv : Specimen data from scaling_key_source_data
- acoustics.csv : Acoustic NASC in 0.5 nmi intervals
- trawl_locations_mace.csv : Lat/Lon locations of all MACE trawl events
- measurements.csv : Length-weight measurements

## Optional arguments

- `ebs` : Boolean flag indicating whether the survey took place in the Eastern Bering Sea.
If `ebs=true`, this function will expect the following two files to be present,
containing data from the Groundfish Assessment Program survey, which are used to scale
acoustic data in the bottom 3 m of the water column:
- trawl_locations_gap.csv : Lat/Lon locations of all bottom trawls
- scaling_gap.csv : Specimen data from racebase.specimen
- `dx` : Set the resolution of the sampling grid; defaults to 10.0 (km).
- `dy` : Set if the N-S resolution is different from the E-W resolution, otherwise will be
equal to `dx`.
- `grid_method` : How to define the survey domain for the purposes of constructing the 
simulation grid. The options are `TransectRibbons()` (the default) or `SurveyHull()`. See
their documentation for details and options.
- `missingstring` : What string(s) in the CSV files should be interpreted as missing data?
Defaults to `[".", "NA"]. (If the files have been downloaded correctly, you should not
need to use this.)
- `transect_order` : `Symbol` indicating which column of `acoustics` to use to define the
direction of each transect for the purposes of defining the domain. Defaults to `:y`,
since most of MACE's surveys have primarily north-south transects. For east-west
transects, use `:x`, and for curving transects use `:log`.

The preprocessing includes the following tasks:

- Excluding NASC from non-scaling classes
- If GAP data are present, merging them with the MACE data into unified `scaling` and
`trawl_locations` data tables.
- Geographic projection of spatial data
- Downscaling acoustic resolution (specified by `dx` and `dy` arguments)
- Calculating survey domain as concave hull of transect ends
- Setting up simulation grid 
- General tidying.

The output from these pre-processing operations is written to new files in the same data
directory:
- scaling.csv : Formatted scaling key data (containing bottom trawls if they were included)
- acoustics_projected.csv : Spatially-projected NASC by interval and scaling class.
"""
function preprocess_survey_data(surveydir; ebs=true, log_ranges=nothing, dx=10.0, dy=dx, 
        missingstring=[".", "NA"],
        grid_method=TransectRibbons(), 
        transect_order=:y)
    scaling_mace = CSV.read(joinpath(surveydir, "scaling_mace.csv"), DataFrame,
        missingstring=missingstring)
    trawl_locations_mace = CSV.read(joinpath(surveydir, "trawl_locations_mace.csv"), DataFrame,
        missingstring=missingstring)

    if ebs
        scaling_gap = CSV.read(joinpath(surveydir, "scaling_gap.csv"), DataFrame,
            missingstring=missingstring)
        scaling = merge_scaling(scaling_mace, scaling_gap)
        trawl_locations_gap = CSV.read(joinpath(surveydir, "trawl_locations_gap.csv"), DataFrame,
            missingstring=missingstring)
        trawl_locations = merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
    else
        scaling = scaling_mace
        trawl_locations = trawl_locations_mace
    end
    
    acoustics = CSV.read(joinpath(surveydir, "acoustics.csv"), DataFrame,
        missingstring=missingstring)

    # Merge "filtered" strata with main ones
    scaling.class .= replace.(scaling.class, "_FILTERED" => "")
    acoustics.class .= replace.(acoustics.class, "_FILTERED" => "")
    scaling_classes = unique(scaling.class)

    if isnothing(log_ranges)
        log_ranges = [tuple(extrema(acoustics.start_vessel_log)...)]
    end

    acoustics = @chain acoustics begin
        @select(:transect,
                :interval, 
                :class, 
                :lat = :start_latitude, 
                :lon = :start_longitude,
                :log = :start_vessel_log,
                :nasc)
        @subset(
            in(scaling_classes).(:class),
            in_intervals(:log, log_ranges),
            abs.(:lon) .< 360,
            abs.(:lat) .< 360
        )
        @by([:transect, :interval, :class, :lat, :lon, :log], :nasc = sum(:nasc))
        unstack([:transect, :interval, :lat, :lon, :log], :class, :nasc, fill=0)
        stack(Not([:transect, :interval, :lat, :lon, :log]), variable_name=:class, value_name=:nasc)
    end
    acoustics.nasc[ismissing.(acoustics.nasc)] .= 0

    utmzone = 3
    lla = LLA.(acoustics.lat, acoustics.lon, 0.0)
    utm = [UTM(x, utmzone, true, wgs84) for x in lla]
    acoustics.x = [u.x / 1e3 for u in utm]
    acoustics.y = [u.y / 1e3 for u in utm]

    surveygrid, surveyhull = get_survey_grid(acoustics, method=grid_method,
        dx=dx, dy=dy, order=transect_order)

    xmin = minimum(acoustics.x)
    ymin = minimum(acoustics.y)

    acoustics = @chain acoustics begin
        # DataFramesMeta.@transform(:x = :x .- xmin, :y = :y .- ymin)
        DataFramesMeta.@transform(
            :x = round.(:x ./ dx) .* dx, 
            :y = round.(:y ./ dy) .* dy
        ) 
        # DataFramesMeta.@transform(:x = :x .+ xmin, :y = :y .+ ymin)
        @by([:transect, :class, :x, :y], 
            :lon = mean(:lon), 
            :lat = mean(:lat), 
            :log = mean(:log),
            :nasc = mean(:nasc)
        )
    end

    trawl_locations = @chain trawl_locations begin
        DataFramesMeta.@transform(:lla = LLA.(:latitude, :longitude, 0.0))
        DataFramesMeta.@transform(:utm = [UTM(x, utmzone, true, wgs84) for x in :lla])
        DataFramesMeta.@transform(
            :x = [u.x / 1e3 for u in :utm], 
            :y = [u.y / 1e3 for u in :utm]
        )
        @subset([in(Point(x, y), surveyhull) for (x, y) in zip(:x, :y)])
    end

    length_weight = CSV.read(joinpath(surveydir, "measurements.csv"), DataFrame,
        missingstring=missingstring)
    rename!(lowercase, length_weight)
    length_weight = @chain length_weight begin
        unstack([:specimen_id, :species_code, :event_id], :measurement_type, :measurement_value)
        dropmissing()
    end


    CSV.write(joinpath(surveydir, "scaling.csv"), scaling)
    CSV.write(joinpath(surveydir, "acoustics_projected.csv"), acoustics)
    CSV.write(joinpath(surveydir, "length_weight.csv"), length_weight)
    CSV.write(joinpath(surveydir, "trawl_locations_projected.csv"), trawl_locations)
    CSV.write(joinpath(surveydir, "surveygrid.csv"), surveygrid)
end
