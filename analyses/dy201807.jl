using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using StatsPlots, StatsPlots.PlotMeasures

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
using .ATBootstrap

survey = "201807"
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
preprocess_survey_data(surveydir, dx=resolution, ebs=true)

(; acoustics, scaling, age_length, length_weight, trawl_locations, surveydomain) = read_survey_files(surveydir)

scaling_classes = unique(scaling.class)
scaling_classes = ["PK1", "PK1_FILTERED", "BT"]

acoustics = @subset(acoustics, in(scaling_classes).(:class), :transect .< 200)

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

surveydata = ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    surveydomain, dA)

atbp = ATBootstrapProblem(surveydata, scaling_classes)

#=
Plotting example conditional simulations
=#
using GeoStats
cp1 = first(atbp.class_problems)
sim_domain = domain(cp1.geoproblem)
sim_coords = coordinates.(sim_domain)

nasc_plots = map(1:6) do _ 
    nasc = ATBootstrap.simulate_nasc(cp1)
    scatter(first.(sim_coords), last.(sim_coords), zcolor=nasc, markershape=:square,
        markerstrokewidth=0, markersize=1.7, legend=false, aspect_ratio=:equal,
        clim=(0, 3000), xlabel="Easting (km)", ylabel="Northing (km)")
end
plot(nasc_plots..., layout=(2, 3), size=(1200, 800), margin=15px)
savefig(joinpath(@__DIR__, "plots", "conditional_nasc.png"))
#=
Plotting example trawl assignments
=#
tl1 = georef(@subset(trawl_locations, :event_id .> 1), (:x, :y))
trawl_coords = coordinates.(domain(tl1))
trawl_assignments_det = ATBootstrap.trawl_assignments(sim_coords, trawl_coords, false)
trawl_assignments_rand= ATBootstrap.trawl_assignments(sim_coords, trawl_coords, true)

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
plot_class_variograms(atbp, legend=:bottomright)

# Check out a couple of conditional simulations
plot_simulated_nasc(atbp, surveydata, size=(1000, 600), markersize=1.3)

# Do the bootstrap uncertainty analysis
results = simulate(atbp, surveydata, nreplicates = 500)
plot_boot_results(results)
CSV.write(joinpath(@__DIR__, "results", "results_$(survey).csv"), results)

n_summary = summarize_bootstrap(results, :n)
biomass_summary = summarize_bootstrap(results, :biomass)

@chain results begin
    stack([:n, :biomass])
    @subset(:species_code .== 21740)
    @by([:i, :variable], :value = sum(:value))
    @by(:variable, :cv = std(:value) / mean(:value))
end

# One-at-a-time error analysis
results_step = stepwise_error(atbp, surveydata; nreplicates = 500)

plot_error_source_by_age(results_step, results, :n)
    
results_totals = merge_results(results, results_step)
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

plot_error_sources(results_totals, plot_title=survey, xticks=0:0.01:0.15, size=(800, 500))
