using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using StatsPlots, StatsPlots.PlotMeasures

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
using .ATBootstrap

survey = "201006"
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
preprocess_survey_data(surveydir, dx=resolution, ebs=true)

(; acoustics, scaling, age_length, length_weight, trawl_locations, surveydomain) = read_survey_files(surveydir)

scaling_classes = unique(scaling.class)
acoustics = @subset(acoustics,
    in(scaling_classes).(:class),
    in(scaling_classes).(:class), :transect .< 200)

@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

surveydata = ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    surveydomain, dA)

atbp = ATBootstrapProblem(surveydata, scaling_classes)

# Inspect the variograms to make sure they look ok
plot_class_variograms(atbp, legend=:bottomright)

# Check out a couple of conditional simulations
plot_simulated_nasc(atbp, surveydata, size=(1000, 600), markersize=2.2)

# Do the bootstrap uncertainty analysis
results = simulate(atbp, surveydata, nreplicates = 500)
plot_boot_results(results)
CSV.write(joinpath(@__DIR__, "results", "results_$(survey).csv"), results)

n_summary = summarize_bootstrap(results, :n)
biomass_summary = summarize_bootstrap(results, :biomass)

# One-at-a-time error analysis
results_step = stepwise_error(atbp, surveydata; nreplicates = 500)

plot_error_source_by_age(results_step, results, :n)
    
results_totals = merge_results(results, results_step)
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

plot_error_sources(results_totals, plot_title=survey, xticks=0:0.01:0.15, size=(800, 500))
