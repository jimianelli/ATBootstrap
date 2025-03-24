using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using StableRNGs
using StatsPlots, StatsPlots.PlotMeasures
using GeoStats

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

survey = "202408"
StableRNGs.seed!(parse(Int, survey))
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2

log_ranges = [(1500, 3800), (4000, 4413.25), (4413.25, 4421), (4425, 9999)]
ATB.preprocess_survey_data(surveydir, dx=resolution, ebs=true, log_ranges=log_ranges,
    grid_method=ATB.TransectRibbons(width=40))

(; acoustics, scaling, age_length, length_weight, trawl_locations, surveygrid) = ATB.read_survey_files(surveydir)

_, domain40 = ATB.get_survey_grid(@subset(acoustics, 1500 .<= :log .<= 5290),
    method=ATB.TransectRibbons(width=40, buffer=0.05))
_, domain30 = ATB.get_survey_grid(@subset(acoustics, 5300 .<= :log .<= 5899),
    method=ATB.TransectRibbons(width=30, buffer=0.066))
_, domain20 = ATB.get_survey_grid(@subset(acoustics, 5900 .<= :log .<= 9999),
    method=ATB.TransectRibbons(width=20, buffer=0.1))
surveydomain = GeometrySet(union(domain40, domain30, domain20))
surveygrid = ATB.grid_domain(surveydomain, resolution)
@df surveygrid scatter(:x, :y)
surveygrid = DataFrames.shuffle(surveygrid)
surveygrid = PointSet(Point(x...) for x in eachrow(surveygrid))

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
ATB.plot_geosim_stats(atbp, surveydata, 500)
savefig(joinpath(@__DIR__, "plots", "conditional_nasc_stats_$(survey).png"))
    
# Inspect the variograms to make sure they look ok
ATB.plot_class_variograms(atbp, legend=:bottomright)

# Check out a couple of conditional simulations
ATB.plot_simulated_nasc(atbp, surveydata, size=(1000, 600), markersize=2.5)

# Do the bootstrap uncertainty analysis
results = ATB.simulate(atbp, surveydata, nreplicates = 500)
ATB.plot_boot_results(results)
CSV.write(joinpath(@__DIR__, "results", "results_$(survey).csv"), results)

n_summary = ATB.summarize_bootstrap(results, :n)
biomass_summary = ATB.summarize_bootstrap(results, :biomass)
sum(n_summary.n)
sum(biomass_summary.biomass)

# One-at-a-time error analysis
results_step = ATB.stepwise_error(atbp, surveydata; nreplicates = 500)

ATB.plot_error_source_by_age(results_step, results, :n)
    
results_totals = ATB.merge_results(results, results_step)
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

ATB.plot_error_sources(results_totals, plot_title=survey, xticks=0:0.01:0.15, size=(800, 500))

@chain results_step begin
    @subset(:species_code .== 21740)
    @by([:i, :added_error], :biomass = sum(:biomass))
    @by(:added_error, :biomass = mean(:biomass))
end

@df @subset(results_totals) density(:biomass, group=:added_error)
@df @subset(results_totals) density(:biomass/1e9, group=:added_error, layout=(3, 4), 
    size=(1000, 800))