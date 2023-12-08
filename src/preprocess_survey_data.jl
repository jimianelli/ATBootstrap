include("nearbottom.jl")

function get_survey_grid(acoustics, k=20, ; transect_width=20.0, dx=10.0, dy=dx)
    w = transect_width / 2
    transect_ends = @chain acoustics begin
        @orderby(:y)
        @by(:transect, 
            :x = [first(:x) + w, first(:x) - w, last(:x) + w, last(:x) - w], 
            :y = [first(:y), first(:y), last(:y), last(:y)])
    end
    v = [[row.x, row.y] for row in eachrow(transect_ends)]
    surveyhull = concave_hull(v, k)

    xgrid = range(round.(extrema(transect_ends.x))..., step=dx)
    ygrid = range(round.(extrema(transect_ends.y))..., step=dy)
    surveydomain = DataFrame([(;x, y) for x in xgrid, y in ygrid
        if in_hull([x, y], surveyhull)])
    return surveydomain, surveyhull
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
    end
    return [scaling_mace1; scaling_gap1]
end


function merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
    return [trawl_locations_mace; trawl_locations_gap]
end

function tryparse_missing(type, str)
    try
        return parse(type, str)
    catch
        return missing
    end
end

"""
    preprocess_survey_data(surveydir[, dx=0.0, dy=dx])

Preprocess the survey data files in directory `surveydir` and save the outputs in 
the same directory as .csv files in standard format. The optional arguments `dx` and `dy`
set the resolution of the sampling grid; they default to 10.0 (km).

This function expects to find the following files, all of which come from running
"download_survey.R":

- scaling_mace.csv : Specimen data from scaling_key_source_data
- acoustics.csv : Acoustic NASC in 0.5 nmi intervals
- trawl_locations_mace.csv : Lat/Lon locations of all MACE trawl events
- measurements.csv : Length-weight measurements

Additionally, if the option `ebs=true`, this function will expect the following two files
to be present containing data from the Groundfish Assessment Program survey, which are 
used to scale acoustic data in the bottom 3 m of the water column:

- trawl_locations_gap.csv : Lat/Lon locations of all bottom trawls
- scaling_gap.csv : Specimen data from racebase.specimen

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
function preprocess_survey_data(surveydir; ebs=true, dx=10.0, dy=dx)
    scaling_mace = CSV.read(joinpath(surveydir, "scaling_mace.csv"), DataFrame)
    trawl_locations_mace = CSV.read(joinpath(surveydir, "trawl_locations_mace.csv"), DataFrame)

    if ebs
        scaling_gap = CSV.read(joinpath(surveydir, "scaling_gap.csv"), DataFrame)
        scaling = merge_scaling(scaling_mace, scaling_gap)
        trawl_locations_gap = CSV.read(joinpath(surveydir, "trawl_locations_gap.csv"), DataFrame)
        trawl_locations = merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
    else
        scaling = scaling_mace
        trawl_locations = trawl_locations_mace
    end
    
    scaling_classes = unique(scaling.class)

    acoustics = CSV.read(joinpath(surveydir, "acoustics.csv"), DataFrame)
    acoustics = @chain acoustics begin
        @select(:transect,
                :interval, 
                :class, 
                :lat = :start_latitude, 
                :lon = :start_longitude,
                :nasc)
        @subset(in(scaling_classes).(:class),
                abs.(:lon) .< 360,
                abs.(:lat) .< 360)
        @by([:transect, :interval, :class, :lat, :lon], :nasc = sum(:nasc))
        unstack([:transect, :interval, :lat, :lon], :class, :nasc, fill=0)
        stack(Not([:transect, :interval, :lat, :lon]), variable_name=:class, value_name=:nasc)
    end
    acoustics.nasc[ismissing.(acoustics.nasc)] .= 0

    utmzone = 3
    lla = LLA.(acoustics.lat, acoustics.lon, 0.0)
    utm = [UTM(x, utmzone, true, wgs84) for x in lla]
    acoustics.x = [u.x / 1e3 for u in utm]
    acoustics.y = [u.y / 1e3 for u in utm]

    surveydomain, surveyhull = get_survey_grid(acoustics, 20, transect_width=20.0,
        dx=dx, dy=dy)

    acoustics = @chain acoustics begin
        DataFramesMeta.@transform(
            :x = round.(:x ./ dx) .* dx, 
            :y = round.(:y ./ dy) .* dy
        ) 
        @by([:transect, :class, :x, :y], 
            :lon = mean(:lon), 
            :lat = mean(:lat), 
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
        @subset([in_hull([x, y], surveyhull) for (x, y) in zip(:x, :y)])
    end

    length_weight = CSV.read(joinpath(surveydir, "measurements.csv"), DataFrame)
    rename!(lowercase, length_weight)
    length_weight = @chain length_weight begin
        unstack([:specimen_id, :species_code, :event_id], :measurement_type, :measurement_value)
        DataFramesMeta.@transform(
            :organism_weight = tryparse_missing.(Float64, :organism_weight)
        )
        dropmissing()
    end


    CSV.write(joinpath(surveydir, "scaling.csv"), scaling)
    CSV.write(joinpath(surveydir, "acoustics_projected.csv"), acoustics)
    CSV.write(joinpath(surveydir, "length_weight.csv"), length_weight)
    CSV.write(joinpath(surveydir, "trawl_locations_projected.csv"), trawl_locations)
    CSV.write(joinpath(surveydir, "surveydomain.csv"), surveydomain)
end


# scaling_gap1 = @chain scaling_gap begin
#     DataFramesMeta.@transform(
#         :survey = :cruise,
#         :ship = :vessel,
#         :event_id = -:haul,
#         :class = "BT",
#         :primary_length = :length ./ 10,
#         :catch_sampling_expansion = 1,
#         :user_defined_expansion = 1,
#         :sample_correction_scalar = 1,
#     )
#     @select(:survey, :ship, :event_id, :class, :species_code, :primary_length, 
#         :catch_sampling_expansion, :user_defined_expansion, :sample_correction_scalar)
# end
# scaling1 = @select(scaling, :survey, :ship, :event_id, :class, :species_code, :primary_length, 
#         :catch_sampling_expansion, :user_defined_expansion, :sample_correction_scalar)