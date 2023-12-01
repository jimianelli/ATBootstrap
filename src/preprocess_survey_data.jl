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
    return surveydomain
end

"""
    preprocess_survey_data(surveydir[, dx=0.0, dy=dx])

Preprocess the survey data files in directory `surveydir` and save the outputs in 
the same directory as .csv files in standard format. The optional arguments `dx` and `dy`
set the resolution of the sampling grid; they default to 10.0 (km).

This function expects to find the following files, all of which come from running
"download_survey.R":

- scaling.csv : Specimen data from scaling_key_source_data
- acoustics.csv : Acoustic NASC in 0.5 nmi intervals
- trawl_locations.csv : Lat/Lon locations of all trawl events
- measurements.csv : Length-weight measurements

The preprocessing includes the following tasks:

- Excluding NASC from non-scaling classes
- Geographic projection of spatial data
- Downscaling acoustic resolution (specified by `dx` and `dy` arguments)
- Calculating survey domain as concave hull of transect ends
- Setting up simulation grid 
- General tidying.

"""
function preprocess_survey_data(surveydir, dx=10.0, dy=dx)
    scaling = CSV.read(joinpath(surveydir, "scaling.csv"), DataFrame)
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

    surveydomain = get_survey_grid(acoustics, 20, transect_width=20.0, dx=dx, dy=dy)

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

    trawl_locations = CSV.read(joinpath(surveydir, "trawl_locations.csv"), DataFrame)
    rename!(lowercase, trawl_locations)
    trawl_locations = @chain trawl_locations begin
        DataFramesMeta.@transform(:lla = LLA.(:eqlatitude, :eqlongitude, 0.0))
        DataFramesMeta.@transform(:utm = [UTM(x, utmzone, true, wgs84) for x in :lla])
        DataFramesMeta.@transform(:x = [u.x / 1e3 for u in :utm], :y = [u.y / 1e3 for u in :utm])
    end

    length_weight = CSV.read(joinpath(surveydir, "measurements.csv"), DataFrame)
    rename!(lowercase, length_weight)
    length_weight = @chain length_weight begin
        unstack([:specimen_id, :event_id], :measurement_type, :measurement_value)
        dropmissing()
    end


    CSV.write(joinpath(surveydir, "acoustics_projected.csv"), acoustics)
    CSV.write(joinpath(surveydir, "length_weight.csv"), length_weight)
    CSV.write(joinpath(surveydir, "trawl_locations_projected.csv"), trawl_locations)
    CSV.write(joinpath(surveydir, "surveydomain.csv"), surveydomain)
end