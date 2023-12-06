using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using GeoStats
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
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
preprocess_survey_data(surveydir, dx=resolution)

(; acoustics, scaling, age_length, length_weight, trawl_locations, surveydomain) = read_survey_files(surveydir)

unique(scaling.class)
# Other classes appear to be extra transects...?
scaling_classes = ["SS1", "SS1_FILTERED", "BT"]
acoustics = @subset(acoustics,
    in(scaling_classes).(:class),
    in(scaling_classes).(:class), :transect .< 200)

@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")


surveydata = ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    surveydomain, dA)

class_problems = map(scaling_classes) do class
    println(class)
    return ATBootstrapProblem(surveydata, class, nlags=15, weightfunc=h -> 1/h)
end

# Inspect the variograms to make sure they look ok
plot_class_variograms(class_problems, legend=:bottomright)


# Check out a couple of conditional simulations
plot_simulated_nasc(class_problems, surveydata, size=(1000, 600))

# Do the bootstrap uncertainty analysis
results = simulate_classes(class_problems, surveydata, nreplicates = 500)
plot_boot_results(results)
CSV.write(joinpath(@__DIR__, "results", "results_$(survey).csv"), results)

results_summary = @chain results begin
    @orderby(:age)
    @by(:age, 
        :biomass_age = mean(:biomass),
        :std = std(:biomass), 
        :cv = std(:biomass) / mean(:biomass) * 100)
end

# One-at-a-time error analysis
results_step = stepwise_error(class_problems, surveydata; nreplicates = 500)

stepwise_summary = @chain results_step begin
    @orderby(:age)
    @by([:added_error, :age], 
        :biomass = mean(:biomass),
        :std = std(:biomass), 
        :cv = std(:biomass) / mean(:biomass) * 100)
    leftjoin(select(results_summary, [:age, :std]), on=:age, makeunique=true)
    DataFramesMeta.@transform(:std_decrease = :std ./ :std_1)
end

@df stepwise_summary plot(:age, :std/1e9, group=:added_error, marker=:o, 
    markerstrokewidth=0, size=(1000, 600), margin=20px,# yscale=:log10,
    xlabel="Age class", ylabel="S.D. (Biomass, MT)", title=survey)
@df results_summary plot!(:age, :std/1e9, linewidth=2, marker=:o, label="All", 
    color=:black)

    
stepwise_totals = @by(results_step, [:added_error, :i], 
    :n = sum(:n), 
    :biomass = sum(:biomass))
results_totals = @by(results, :i, 
    :n = sum(:n), 
    :biomass = sum(:biomass),
    :added_error = "All")

results_totals = @chain [results_totals; stepwise_totals] begin
    leftjoin(error_labels, on=:added_error)
end
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

stds_boot = map(1:1000) do i
    df = resample_df(results_totals)
    @by(df, :error_label, 
        :n_cv = iqr(:n) / mean(:n) ,
        :biomass_cv = iqr(:biomass) / mean(:biomass))
end 
stds_boot = vcat(stds_boot...)

plot_error_sources(stds_boot, plot_title=survey)
