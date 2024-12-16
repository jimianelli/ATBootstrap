using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using Random
using StatsPlots, StatsPlots.PlotMeasures

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

survey = "201807"
Random.seed!(parse(Int, survey))
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
log_ranges = [(300, 3527.99), (3528, 3864.99), (5145, 7123.99), (7486, 25000)]

ATB.preprocess_survey_data(surveydir, dx=resolution, ebs=true, log_ranges=log_ranges)
(; acoustics, scaling, age_length, length_weight, trawl_locations, surveygrid) = ATB.read_survey_files(surveydir)

scaling_classes = unique(scaling.class)
scaling_classes = ["PK1", "BT"]
acoustics = @subset(acoustics,
    in(scaling_classes).(:class))

@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/200, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

p_xsects = @df acoustics scatter(:x, :y, group=:class, markersize=:nasc/200,
    alpha=0.5, title="(a)", titlealign=:left)
p_trawls = @df @subset(trawl_locations, :event_id .< 0) scatter(:x, :y, label="Bottom",
    markersize=2, title="(b)", titlealign=:left)
@df @subset(trawl_locations, :event_id .> 0) scatter!(p_trawls, :x, :y, label="Midwater",
    markersize=3)
plot(p_xsects, p_trawls, xlabel="Easting (km)", ylabel="Northing (km)", aspect_ratio=:equal,
    markerstrokewidth=0, xlims=(-250, 800), size=(700, 400), dpi=300)
savefig(joinpath(@__DIR__, "plots", "DY201807_maps.png"))

surveydata = ATB.ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    surveygrid, dA)

atbp = ATB.ATBootstrapProblem(surveydata)

sim_dists = ATB.zdists(atbp)
sim_dists.survey .= survey
CSV.write(joinpath(@__DIR__, "results", "zdists_$(survey).csv"),
    select(sim_dists, [:survey, :class, :zdist]))
    
#=
Plotting example conditional simulations
=#
using GeoStats
cp1 = first(atbp.class_problems)
sim_domain = ATB.solution_domain(cp1)

scatter(sim_domain.x, sim_domain.y, markersize=2)
scatter!(acoustics.x, acoustics.y, markersize=2)

nasc_plots = map(1:6) do _ 
    nasc = ATB.simulate_nasc(cp1)
    scatter(sim_domain.x, sim_domain.y, zcolor=nasc, markershape=:square,
        markerstrokewidth=0, markersize=1.7, legend=false, aspect_ratio=:equal,
        clim=(0, 3000), xlabel="Easting (km)", ylabel="Northing (km)")
end
plot(nasc_plots..., layout=(2, 3), size=(1200, 800), margin=15px)
savefig(joinpath(@__DIR__, "plots", "conditional_nasc.png"))

ATB.plot_geosim_stats(atbp, surveydata, 500)
savefig(joinpath(@__DIR__, "plots", "conditional_nasc_stats_$(survey).png"))

#=
Plotting example trawl assignments
=#
tl1 = georef(@subset(trawl_locations, :event_id .> 1), (:x, :y))
sim_coords = ATB.svector_coords.(cp1.geosetup.domain)
trawl_coords = ATB.svector_coords.(domain(tl1))
trawl_assignments_det = ATB.trawl_assignments(sim_coords, trawl_coords, false)
trawl_assignments_rand= ATB.trawl_assignments(sim_coords, trawl_coords, true)

ms = 1.8
pal = :Set1_9
p1 = scatter(first.(sim_coords), last.(sim_coords), color=trawl_assignments_det,
    markershape=:square, markerstrokewidth=0, markersize=ms, legend=false, palette=pal,
    title="(a)", titlealign=:left)
scatter!(p1, first.(trawl_coords), last.(trawl_coords), color=:black, markersize=4)
p2 = scatter(first.(sim_coords), last.(sim_coords), color=trawl_assignments_rand,
    markershape=:square, markerstrokewidth=0, markersize=ms, legend=false, palette=pal,
    title="(b)", titlealign=:left)
scatter!(p2, first.(trawl_coords), last.(trawl_coords), color=:black, markersize=4)
plot(p1, p2, xlabel="Easting (km)", ylabel="Northing (km)", aspect_ratio=:equal,
    markerstrokewidth=0, xlims=(-250, 800), size=(900, 400), dpi=300, 
    bottom_margin=15px)
savefig(joinpath(@__DIR__, "plots", "trawl_assignments.png"))

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

@chain results begin
    stack([:n, :biomass])
    @subset(:species_code .== 21740)
    @by([:i, :variable], :value = sum(:value))
    @by(:variable, :cv = std(:value) / mean(:value))
end

# One-at-a-time error analysis
results_step = ATB.stepwise_error(atbp, surveydata; nreplicates = 500)

ATB.plot_error_source_by_age(results_step, results, :n)
    
results_totals = ATB.merge_results(results, results_step)
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

ATB.plot_error_sources(results_totals, plot_title=survey, xticks=0:0.01:0.15, size=(800, 500))
