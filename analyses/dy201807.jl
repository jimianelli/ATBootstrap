using CSV, DataFrames, DataFramesMeta
using GeoStats, GeoStatsPlots
using Statistics, StatsBase
using Distributions
using Random
using ConcaveHull
using StatsPlots, StatsPlots.PlotMeasures

using Revise
includet(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))

survey = "201807"
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
resolution = 10.0 # km
const km2nmi = 1 / 1.852

acoustics, scaling, trawl_locations, surveydomain = read_survey_files(surveydir)

scaling_classes = unique(scaling.class)
scaling_classes = ["PK1", "PK1_FILTERED"]

acoustics = @chain acoustics begin
    @subset(in(scaling_classes).(:class), :transect .< 200)
    DataFramesMeta.@transform(:x = round.(:x, digits=-1), :y = round.(:y, digits=-1))
    @by([:transect, :class, :x, :y], 
        :lon = mean(:lon), :lat = mean(:lat), :nasc = mean(:nasc))
end
@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

surveydata = ATSurveyData(acoustics, scaling, trawl_locations, surveydomain)

cal_error = 0.1 # dB
dA = (resolution * km2nmi)^2
class_problems = map(scaling_classes) do class
    println(class)
    return ATBootstrapProblem(surveydata, class, cal_error, dA)
end

simdomain = solution_domain(class_problems[1])
sim_fields = [nonneg_lusim(p) for p in class_problems]
sim_plots = map(enumerate(sim_fields)) do (i, x)
    plot(simdomain, zcolor=x, clims=(0, quantile(x, 0.999)), 
        markerstrokewidth=0, markershape=:square, title=string(scaling_classes[i]),
        markersize=2.4, xlabel="Easting (km)", ylabel="Northing (km)")
    df = @subset(acoustics, :class .== scaling_classes[i])
    scatter!(df.x, df.y, color=:white, markersize=df.nasc*3e-3, alpha=0.3,
        markerstrokewidth=0)
end
plot(sim_plots..., size=(1000, 1000))
# i = 1
# x = nonneg_lusim(params[i], zdists[i])
# plot(domain(sol), zcolor=x, markerstrokewidth=0, markershape=:square, clims=(0, quantile(x, 0.999)), 
#     title=string(scaling_classes[i]), markersize=4.2, background_color=:black, 
#     xlabel="Easting (km)", ylabel="Northing (km)", size=(1000, 1000))
unique(trawl_locations.event_id)
unique(scaling.event_id)

results = simulate_classes(class_problems, surveydata, nreplicates=1000)
# results = simulate(class_problems[1], surveydata)
CSV.write(joinpath(@__DIR__, "results_$(survey).csv"), results)

@df results density(:n_age/1e9, group=:age, #xlims=(0, 8),
    fill=true, alpha=0.7, ylims=(0, 25), palette=:Paired_10,
    xlabel="Billions of fish", ylabel="Probability density",
    title=survey)

@df results boxplot(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    xlabel="Age class", ylabel="Abundance (billions)")

@df results density(:biomass_age/1e9, group=:age, #xlims=(0, 4),
    fill=true, alpha=0.7, ylims=(0, 50), palette=:Paired_10,
    xlabel="Million tons", ylabel="Probability density")

@df results boxplot(:age, :biomass_age/1e9, group=:age, palette=:Paired_10,
    xlabel="Age class", ylabel="Million tons")

results_summary = @chain results begin
    @orderby(:age)
    @by(:age, 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
end



results_step = stepwise_error_removal(class_problems, surveydata; nreplicates = 1000)


stepwise_summary = @chain results_step begin
    @orderby(:age)
    @by([:eliminated_error, :age], 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
    # @select(:age, :eliminated_error, :cv_age)
    # unstack(:age, :eliminated_error, :cv_age)
end

@df stepwise_summary plot(:age, :std_age, group=:eliminated_error, marker=:o, 
    markerstrokewidth=0, size=(800, 600),
    xlabel="Age class", ylabel="C.V. (%)", title=survey)
@df results_summary plot!(:age, :std_age, linewidth=2, marker=:o, label="None", 
    color=:black)


res1 = simulate_classes(class_problems, surveydata, nreplicates=1000,
    bs = BootSpecs(calibration=false))
res1_sumary = @chain res1 begin
    @orderby(:age)
    @by(:age, 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
end

@df res1_sumary plot!(:age, :std_age, linewidth=2, marker=:o, 
    label="trawl_assignments2")