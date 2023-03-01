using CSV, DataFrames, DataFramesMeta
using GeoStats, GeoStatsPlots
using Statistics, StatsBase
using Distributions
using Random
using ConcaveHull
using StatsPlots, StatsPlots.PlotMeasures

using Revise
includet(joinpath(@__DIR__, "src", "ATBootstrap.jl"))

surveydir = joinpath(@__DIR__, "surveydata", "202207")
resolution = 10.0 # km
const km2nmi = 1 / 1.852

acoustics = CSV.read(joinpath(surveydir, "acoustics_projected.csv"), DataFrame)
trawl_locations = CSV.read(joinpath(surveydir, "trawl_locations_projected.csv"), DataFrame)
scaling = CSV.read(joinpath(surveydir, "scaling.csv"), DataFrame)
surveydomain = CSV.read(joinpath(surveydir, "surveydomain.csv"), DataFrame)
surveydomain = shuffle(surveydomain) # this seems to fix the issue with directional artifacts
surveydomain =  PointSet(Matrix(surveydomain)')
scaling_classes = unique(scaling.class)


acoustics = @chain acoustics begin
    @subset(in(scaling_classes).(:class), :transect .< 200)
    DataFramesMeta.@transform(:x = round.(:x, digits=-1), :y = round.(:y, digits=-1))
    @by([:transect, :class, :x, :y], 
        :lon = mean(:lon), :lat = mean(:lat), :nasc = mean(:nasc))
end
@df acoustics scatter(:x, :y, group=:class, markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

surveydata = ATSurveyData(acoustics, scaling, trawl_locations, surveydomain)

cal_error = 0.1 # dB
dA = (resolution * km2nmi)^2
class_problems = map(scaling_classes) do class
    println(class)
    return ATBootstrapProblem(surveydata, class, cal_error, dA, nreplicates=1000)
end

simdomain = solution_domain(class_problems[1])
sim_fields = [nonneg_lusim(p) for p in class_problems]
sim_plots = map(enumerate(sim_fields)) do (i, x)
    plot(simdomain, zcolor=x, clims=(0, quantile(x, 0.999)), 
        markerstrokewidth=0, markershape=:square, title=string(scaling_classes[i]),
        markersize=1.5, xlabel="Easting (km)", ylabel="Northing (km)")
    df = @subset(acoustics, :class .== scaling_classes[i])
    scatter!(df.x, df.y, color=:white, markersize=df.nasc*1e-3, alpha=0.3,
        markerstrokewidth=0)
end
plot(sim_plots..., size=(1000, 1000))
# i = 1
# x = nonneg_lusim(params[i], zdists[i])
# plot(domain(sol), zcolor=x, markerstrokewidth=0, markershape=:square, clims=(0, quantile(x, 0.999)), 
#     title=string(scaling_classes[i]), markersize=4.2, background_color=:black, 
#     xlabel="Easting (km)", ylabel="Northing (km)", size=(1000, 1000))

simulate(class_problems[1], surveydata)

results = simulate_classes(class_problems, surveydata)


@df results density(:n_age/1e9, group=:age, xlims=(0, 8),
    fill=true, alpha=0.7, ylims=(0, 30), palette=:Paired_10,
    xlabel="Billions of fish", ylabel="Probability density")
@df results density(:n_age/1e9, group=:age, palette=:Paired_10, fill=true, 
    alpha=0.7, ylims=(0, 0.01),
    xlabel="Billions of fish", ylabel="Probability density")

@df results boxplot(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    ylabel="Billions of fish")
@df results boxplot(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    ylims=(0, 10), ylabel="Billions of fish")

@by(results, :age, 
    :n_age = mean(:n_age),
    :std_age = std(:n_age), 
    :cv_age = std(:n_age) / mean(:n_age) * 100)


@df results density(:biomass_age/1e9, group=:age, xlims=(0, 4),
    fill=true, alpha=0.7, ylims=(0, 30), palette=:Paired_10,
    xlabel="Million tons", ylabel="Probability density")
@df results density(:biomass_age/1e9, group=:age, palette=:Paired_10, fill=true, 
    alpha=0.7, ylims=(0, 0.2),
    xlabel="Million tons", ylabel="Probability density")

@df results boxplot(:age, :biomass_age/1e9, group=:age, palette=:Paired_10,
    ylabel="Million tons")
@df results boxplot(:age, :biomass_age/1e9, group=:age, palette=:Paired_10,
    ylims=(0, 10), ylabel="Million tons")

@by(results, :age, 
    :biomass_age = mean(:biomass_age),
    :std_age = std(:biomass_age), 
    :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
