using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using StableRNGs
using StatsPlots, StatsPlots.PlotMeasures, ColorSchemes
import GLMakie, GeoMakie
using GADM

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

survey = "201807"
StableRNGs.seed!(parse(Int, survey))
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

i_bt = contains.(trawl_locations.haul_id, "GAP")

geodf = [
    DataFrame(GADM.get("USA", "Alaska"));
    DataFrame(GADM.get("RUS", "Chukot"))
]

fig = GLMakie.Figure()
ga1 = GeoMakie.GeoAxis(fig[1,1], 
    dest="+proj=ortho +lon_0=-170 +lat_0=57",
    # dest = "+proj=utm +zone=3",
    xgridcolor=(:black, 0.25), ygridcolor=(:black, 0.25),
    xlabel = "Longitude", ylabel="Latitude", title="(a)", titlealign=:left,
)
ga2 = GeoMakie.GeoAxis(fig[1,2], 
    dest="+proj=ortho +lon_0=-170 +lat_0=57",
    # dest = "+proj=utm +zone=3",
    xgridcolor=(:black, 0.25), ygridcolor=(:black, 0.25), yaxisposition=:right,
    xlabel = "Longitude", ylabel="Latitude", title="(b)", titlealign=:left,
)
GLMakie.poly!(ga1, geodf.geom, color=:grey)
GLMakie.poly!(ga2, geodf.geom, color=:grey)
colors = [:dodgerblue, :darkorange]
for (i, gdf) in enumerate(groupby(acoustics, :class))
    GLMakie.scatter!(ga1, gdf.lon, gdf.lat, label=first(gdf.class),
        markersize=gdf.nasc/200 .+ 1, alpha=0.5, color=colors[i])
end
GLMakie.scatter!(ga2, trawl_locations[i_bt,:].longitude, trawl_locations[i_bt, :].latitude,
    markersize=3, color=colors[1], label="Bottom")
GLMakie.scatter!(ga2, trawl_locations[.!i_bt,:].longitude, trawl_locations[.!i_bt, :].latitude,
    markersize=5, color=colors[2], label="Midwater")

for ax in [ga1, ga2]
    GLMakie.xlims!(ax, -181, -158)
    GLMakie.ylims!(ax, 53, 63)
    ax.xticks = -185:5:-150
    ax.titlefont = :regular
    # ax.yticks = 50:5:65
end
elem_1 = GLMakie.MarkerElement(color=colors[1], marker=:circle)
elem_2 = GLMakie.MarkerElement(color=colors[2], marker=:circle)
GLMakie.Legend(fig[1,3], 
    [elem_1, elem_2],
    ["Bottom", "Midwater"]
)
GLMakie.save(joinpath(@__DIR__, "plots", "DY201807_maps.png"), fig,
    size=(700, 300), px_per_unit=4)

# p_xsects = @df acoustics scatter(:x, :y, group=:class, markersize=:nasc/200,
#     alpha=0.5, title="(a)", titlealign=:left)
# p_trawls = @df trawl_locations[i_bt, :] scatter(:x, :y, label="Bottom",
#     markersize=2, title="(b)", titlealign=:left)
# @df trawl_locations[.! i_bt, :] scatter!(p_trawls, :x, :y, label="Midwater",
#     markersize=3)
# plot(p_xsects, p_trawls, xlabel="Easting (km)", ylabel="Northing (km)", aspect_ratio=:equal,
#     markerstrokewidth=0, xlims=(-250, 800), size=(700, 400), dpi=300)
# savefig(joinpath(@__DIR__, "plots", "DY201807_maps.png"))

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
cp1 = atbp.class_problems[2]
sim_domain = ATB.solution_domain(cp1)

scatter(sim_domain.x, sim_domain.y, markersize=2)
scatter!(acoustics.x, acoustics.y, markersize=2)

clims=(0, 2000)
cscheme = :matter
p_xsects = @df @subset(acoustics, :class.=="PK1") scatter(:x, :y,
    markersize=:nasc/300 .+ 1, markerstrokewidth=0, label="",
    zcolor=:nasc, clims=clims, legend=false, color=cscheme,
    title="(a)", titlealign=:left, ylabel="Northing (km)")

titles = "(" .* ('b':'f') .* ")"
nasc_plots = map(1:5) do i
    nasc = ATB.simulate_nasc(cp1)
    p = scatter(sim_domain.x, sim_domain.y, zcolor=nasc, markershape=:square,
        markerstrokewidth=0, markersize=1.7, legend=false, aspect_ratio=:equal,
        color=cscheme,
        clim=clims, xlabel="", ylabel="", title=titles[i], titlealign=:left)
    if i == 3
        ylabel!(p, "Northing (km)")
    end
    if 3 .<= i .<= 5
        xlabel!(p, "Easting (km)")
    end
    p
end
h2 = scatter([0,0], [0,0], zcolor=[clims...], clims=clims, color=cscheme,
    xlims=(1,1.1), label="", colorbar_title="NASC (m² nmi⁻²)", framestyle=:none)
l = @layout [grid(2, 3) a{0.01w}]
plot([p_xsects; nasc_plots; h2]..., layout=l, link=:all, size=(1200, 800), margin=15px)
savefig(joinpath(@__DIR__, "plots", "conditional_nasc.png"))

ATB.plot_geosim_stats(atbp, surveydata, 500)
savefig(joinpath(@__DIR__, "plots", "conditional_nasc_stats_$(survey).png"))

#=
Plotting example trawl assignments
=#
tl1 = georef(trawl_locations[.! i_bt, :], (:x, :y))
sim_coords = ATB.svector_coords.(cp1.geosetup.domain)
trawl_coords = ATB.svector_coords.(domain(tl1))
trawl_assignments_det = ATB.trawl_assignments(sim_coords, trawl_coords, false)
trawl_assignments_rand= ATB.trawl_assignments(sim_coords, trawl_coords, true)

ms = 1.63
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

# Inspect the vari3333ograms to make sure they look ok
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
