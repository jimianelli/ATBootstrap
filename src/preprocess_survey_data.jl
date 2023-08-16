# using StatsPlots
# using CSV, DataFrames, DataFramesMeta
# using Geodesy
# using ConcaveHull
# using Statistics
# using LinearAlgebra

function preprocess_survey_data(surveydir, dx=10.0, dy=dx)
    # surveydir = joinpath("surveydata", "200707")
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

    # scaling = CSV.read(joinpath(surveydir, "scaling.csv"), DataFrame)
    # event_class = @by(scaling, :event_id, :class = length(unique(class))) # need to figure this out too

    trawl_locations = CSV.read(joinpath(surveydir, "trawl_locations.csv"), DataFrame)
    rename!(lowercase, trawl_locations)
    trawl_locations = @chain trawl_locations begin
        DataFramesMeta.@transform(:lla = LLA.(:eqlatitude, :eqlongitude, 0.0))
        DataFramesMeta.@transform(:utm = [UTM(x, utmzone, true, wgs84) for x in :lla])
        DataFramesMeta.@transform(:x = [u.x / 1e3 for u in :utm], :y = [u.y / 1e3 for u in :utm])
    end
    trawl_locations[!, :edsu_idx] .= 0
    for i in 1:nrow(trawl_locations)
        y = trawl_locations.y[i]
        trawl_locations[i, :edsu_idx] = argmin(norm.(y - s for s in acoustics.y))
    end

    length_weight = CSV.read(joinpath(surveydir, "measurements.csv"), DataFrame)
    rename!(lowercase, length_weight)
    length_weight = @chain length_weight begin
        unstack([:specimen_id, :event_id], :measurement_type, :measurement_value)
        dropmissing()
    end
    # @df length_weight scatter(:fork_length, :organism_weight)

    transect_ends = @chain acoustics begin
        @orderby(:y)
        @by(:transect, 
            :x = [first(:x) + 10, first(:x) - 10, last(:x) + 10, last(:x) - 10], 
            :y = [first(:y), first(:y), last(:y), last(:y)])
    end
    v = [[row.x, row.y] for row in eachrow(transect_ends)]
    surveyhull = concave_hull(v, 20)
    # @df acoustics scatter(:x, :y, aspect_ratio=:equal)
    # plot!(surveyhull)

    dx = 10.0
    dy = 10.0
    xgrid = range(round.(extrema(transect_ends.x))..., step=dx)
    ygrid = range(round.(extrema(transect_ends.y))..., step=dy)
    surveydomain = DataFrame([(;x, y) for x in xgrid, y in ygrid
        if in_hull([x, y], surveyhull)])

    CSV.write(joinpath(surveydir, "acoustics_projected.csv"), acoustics)
    CSV.write(joinpath(surveydir, "length_weight.csv"), length_weight)
    CSV.write(joinpath(surveydir, "trawl_locations_projected.csv"), trawl_locations)
    CSV.write(joinpath(surveydir, "surveydomain.csv"), surveydomain)
end