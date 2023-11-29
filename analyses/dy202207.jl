using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using GeoStats, GeoStatsPlots
using Statistics, StatsBase
using Distributions
using Random
using ConcaveHull
using StatsPlots, StatsPlots.PlotMeasures

using Revise
includet(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
using .ATBootstrap

survey = "202207"
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
preprocess_survey_data(surveydir, resolution)
const km2nmi = 1 / 1.852

(; acoustics, scaling, age_length, length_weight, trawl_locations, domain) = read_survey_files(surveydir)

unique(scaling.class)
# Other classes appear to be extra transects...?
scaling_classes = ["SS1", "SS1_FILTERED"]
acoustics = @subset(acoustics, in(scaling_classes).(:class))

acoustics = @chain acoustics begin
    @subset(in(scaling_classes).(:class), :transect .< 200)
    DataFramesMeta.@transform(:x = round.(:x, digits=-1), :y = round.(:y, digits=-1))
    @by([:transect, :class, :x, :y], 
        :lon = mean(:lon), :lat = mean(:lat), :nasc = mean(:nasc))
end
@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")


surveydata = ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    domain)

class_problems = map(scaling_classes) do class
    println(class)
    return ATBootstrapProblem(surveydata, class, dA, nlags=15, weightfunc=h -> 1/h)
end

pp = map(class_problems) do cp
    plot(cp.variogram.empirical, title=cp.class)
    plot!(cp.variogram.model, xlims=(0, 200))
end
plot(pp..., size=(1000, 800))

# Some analysis of variability in NASC
/(quantile(acoustics.nasc[acoustics.nasc .> 0], [0.95, 0.05])...)
x = sort(acoustics.nasc, rev=true)
plot(cumsum(x) / sum(x))
findfirst(cumsum(x) / sum(x) .> 0.5) / length(x)

vgm = class_problems[1].variogram.model
plot(h -> sqrt(vgm(h)) / mean(acoustics.nasc), 0, 100)
vline!([7])
sqrt(vgm(1.852 * 10)) / mean(acoustics.nasc)

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
unique(trawl_locations.event_id)
unique(scaling.event_id)

results = simulate_classes(class_problems, surveydata, nreplicates = 500)
results = @subset(results, :age .!= "00")
CSV.write(joinpath(@__DIR__, "results", "results_$(survey).csv"), results)

p_abundance = @df results violin(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    xlabel="Age class", ylabel="Million tons", legend=false);
p_biomass = @df results violin(:age, :biomass_age/1e9, group=:age, palette=:Paired_10,
    xlabel="Age class", ylabel="Million tons");
plot(p_abundance, p_biomass, size=(900, 400), margin=15px)

results_summary = @chain results begin
    @orderby(:age)
    @by(:age, 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
end

# One-at-a-time error analysis
results_step = stepwise_error(class_problems, surveydata; nreplicates = 500)

stepwise_summary = @chain results_step begin
    @orderby(:age)
    @by([:added_error, :age], 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
    leftjoin(select(results_summary, [:age, :std_age]), on=:age, makeunique=true)
    DataFramesMeta.@transform(:std_decrease = :std_age ./ :std_age_1)
end

@df stepwise_summary plot(:age, :std_age/1e9, group=:added_error, marker=:o, 
    markerstrokewidth=0, size=(1000, 600), margin=20px,# yscale=:log10,
    xlabel="Age class", ylabel="S.D. (Biomass, MT)", title=survey)
@df results_summary plot!(:age, :std_age/1e9, linewidth=2, marker=:o, label="All", 
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
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

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
plot(p1, p2, layout=(2,1), size=(700, 600), legend=false, xlims=(-0.005, 0.20),
    ylabel="Error source")
