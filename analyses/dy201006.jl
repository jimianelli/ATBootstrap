using CSV, DataFrames, DataFramesMeta
using GeoStats, GeoStatsPlots
using Statistics, StatsBase
using Distributions
using Random
using ConcaveHull
using StatsPlots, StatsPlots.PlotMeasures

using Revise
includet(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
using .ATBootstrap

survey = "201006"
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
resolution = 10.0 # km
const km2nmi = 1 / 1.852

acoustics, scaling, length_weight, trawl_locations, surveydomain = read_survey_files(surveydir)

scaling_classes = unique(scaling.class)

acoustics = @chain acoustics begin
    @subset(in(scaling_classes).(:class))
    DataFramesMeta.@transform(:x = round.(:x, digits=-1), :y = round.(:y, digits=-1))
    @by([:transect, :class, :x, :y], 
        :lon = mean(:lon), :lat = mean(:lat), :nasc = mean(:nasc))
end
@df acoustics scatter(:x, :y, group=:class, markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

@df acoustics plot(:y, :nasc, group=:transect, legend=false)


surveydata = ATSurveyData(acoustics, scaling, length_weight, trawl_locations, surveydomain)

cal_error = 0.1 # dB
dA = (resolution * km2nmi)^2
class_problems = map(scaling_classes) do class
    println(class)
    return ATBootstrapProblem(surveydata, class, dA, nlags=15, weightfunc=h -> 1/h)
end

pp = map(class_problems) do cp
    plot(cp.variogram.empirical, title=cp.class)
    plot!(cp.variogram.model, xlims=(0, 200))
end
plot(pp..., size=(1000, 800))

simdomain = solution_domain(class_problems[1])
sim_fields = [nonneg_lusim(p) for p in class_problems]
sim_plots = map(enumerate(sim_fields)) do (i, x)
    plot(simdomain, zcolor=x, clims=(0, quantile(x, 0.999)), 
        markerstrokewidth=0, markershape=:square, title=string(scaling_classes[i]),
        markersize=2.2, xlabel="Easting (km)", ylabel="Northing (km)")
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

results = simulate_classes(class_problems, surveydata)
results = @subset(results, :age .!= "00")
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

results_step = stepwise_error(class_problems, surveydata; remove=false, nreplicates = 500)
results_step = @subset(results_step, :age .!= "00")

stepwise_summary = @chain results_step begin
    @orderby(:age)
    @by([:added_error, :age], 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
    leftjoin(select(results_summary, [:age, :std_age]), on=:age, makeunique=true)
    DataFramesMeta.@transform(:std_decrease = :std_age ./ :std_age_1)
end

@df stepwise_summary plot(:age, :std_age, group=:added_error, marker=:o, 
    markerstrokewidth=0, size=(800, 600),# yscale=:log10,
    xlabel="Age class", ylabel="C.V. (%)", title=survey)
@df results_summary plot!(:age, :std_age, linewidth=2, marker=:o, label="None", 
    color=:black)

stepwise_totals = @by(results_step, [:added_error, :i], 
    :n = sum(:n_age), 
    :biomass = sum(:biomass_age))
results_totals = @by(results, :i, 
    :n = sum(:n_age), 
    :biomass = sum(:biomass_age),
    :added_error = "All")

results_totals = @chain [results_totals; stepwise_totals] begin
    leftjoin(error_labels, on=:added_error)
end
CSV.write("analyses/stepwise_error_$(survey).csv", results_totals)

@df results_totals boxplot(:error_label, :n)

stds_boot = map(1:1000) do i
    df = resample_df(results_totals)
    @by(df, :error_label, 
        :n_cv = iqr(:n) / mean(:n) ,
        :biomass_cv = iqr(:biomass) / mean(:biomass))
end 
stds_boot = vcat(stds_boot...)

p1 = @df stds_boot boxplot(:error_label, :n_cv, permute=(:x, :y), xflip=true,
    outliers=false, title=survey, ylabel="C.V. (Numbers)");
p2 = @df stds_boot boxplot(:error_label, :biomass_cv, permute=(:x, :y), xflip=true,
    outliers=false, ylabel="C.V. (Biomass)");
plot(p1, p2, layout=(2,1), size=(700, 600), legend=false, xlims=(-0.005, 0.15),
    ylabel="Error source")

