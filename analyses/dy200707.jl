using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using Random
using StatsPlots, StatsPlots.PlotMeasures

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

survey = "200707"
Random.seed!(parse(Int, survey))
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
log_ranges = [(0, 8275.49), (8358.5, 9072.5), (9271, 9685)]

ATB.preprocess_survey_data(surveydir, dx=resolution, ebs=true, log_ranges=log_ranges)
(; acoustics, scaling, age_length, length_weight, trawl_locations, surveygrid) = ATB.read_survey_files(surveydir)

@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

surveydata = ATB.ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    surveygrid, dA)

atbp = ATB.ATBootstrapProblem(surveydata)

sim_dists = ATB.zdists(atbp)
sim_dists.survey .= survey
CSV.write(joinpath(@__DIR__, "results", "zdists_$(survey).csv"),
    select(sim_dists, [:survey, :class, :zdist]))
    
# Inspect the variograms to make sure they look ok
ATB.plot_class_variograms(atbp, legend=:bottomright)

# Check out a couple of conditional simulations
ATB.plot_simulated_nasc(atbp, surveydata, size=(1000, 600), markersize=1.3)

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
