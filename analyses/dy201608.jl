using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using Random
using StatsPlots, StatsPlots.PlotMeasures

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

survey = "201608"
Random.seed!(parse(Int, survey))
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
ATB.preprocess_survey_data(surveydir, dx=resolution, ebs=true)

(; acoustics, scaling, age_length, length_weight, trawl_locations, surveydomain) = ATB.read_survey_files(surveydir)

scaling_classes = unique(scaling.class)
acoustics = @subset(acoustics,
    in(scaling_classes).(:class),
    in(scaling_classes).(:class), :transect .< 200)

@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

surveydata = ATB.ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    surveydomain, dA)

atbp = ATB.ATBootstrapProblem(surveydata, scaling_classes)

# Inspect the variograms to make sure they look ok
ATB.plot_class_variograms(atbp, legend=:bottomright)

# Check out a couple of conditional simulations
ATB.plot_simulated_nasc(atbp, surveydata, size=(1000, 600), markersize=2.2)

cp1 = atbp.class_problems[1]
nasc = ATB.simulate_nasc(cp1)
dom = ATB.solution_domain(cp1)
scatter(dom.x, dom.y, zcolor=nasc, markerstrokewidth=0, legend=false, title="2016")
scatter!(acoustics.x, acoustics.y, color=:white, markerstrokewidth=0, markersize=2)

ATB.plot_geosim_stats(atbp, surveydata)
# Do the bootstrap uncertainty analysis
results = ATB.simulate(atbp, surveydata, nreplicates = 500)
ATB.plot_boot_results(results)
CSV.write(joinpath(@__DIR__, "results", "results_$(survey).csv"), results)

n_summary = ATB.summarize_bootstrap(results, :n)
biomass_summary = ATB.summarize_bootstrap(results, :biomass)

# One-at-a-time error analysis
results_step = ATB.stepwise_error(atbp, surveydata; nreplicates = 500)

ATB.plot_error_source_by_age(results_step, results, :n)
    
results_totals = ATB.merge_results(results, results_step)
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

ATB.plot_error_sources(results_totals, plot_title=survey, xticks=0:0.01:0.15, size=(800, 500))
