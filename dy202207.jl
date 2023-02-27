using CSV, DataFrames, DataFramesMeta
using GeoStats, GeoStatsPlots
using NearestNeighbors
using Distances
using Statistics, StatsBase
using KernelDensity
using Distributions
using ConcaveHull
using Random
using LinearAlgebra
using StatsPlots, StatsPlots.PlotMeasures
using ProgressMeter

using Revise
includet(joinpath(@__DIR__, "src", "ATBootstrap.jl"))

surveydir = joinpath(@__DIR__, "surveydata", "202207")
resolution = 10.0 # km
const km2nmi = 1 / 1.852

acoustics = CSV.read(joinpath(surveydir, "acoustics_projected.csv"), DataFrame)
trawl_locations = CSV.read(joinpath(surveydir, "trawl_locations_projected.csv"), DataFrame)
scaling = CSV.read(joinpath(surveydir, "scaling.csv"), DataFrame)
surveydomain = CSV.read(joinpath(surveydir, "surveydomain.csv"), DataFrame)
surveydomain = shuffle(surveydomain) # this seems to fix the issue with directional artifacts
surveydomain =  PointSet(Matrix(surveydomain)')

scaling_classes = unique(scaling.class)

# spatial binning

acoustics = @chain acoustics begin
    @subset(in(scaling_classes).(:class), :transect .< 200)
    DataFramesMeta.@transform(:x = round.(:x, digits=-1), :y = round.(:y, digits=-1))
    @by([:transect, :class, :x, :y], 
        :lon = mean(:lon), :lat = mean(:lat), :nasc = mean(:nasc))
end
@df acoustics scatter(:x, :y, group=:class, markersize=:nasc/500, markerstrokewidth=0, alpha=0.5)
@df trawl_locations scatter!(:x, :y, label="")

surveydata = (;acoustics, scaling, trawl_locations, scaling_classes)

variograms, problems = define_conditional_sims(acoustics, scaling_classes)

p = plot();
for (i, v) in enumerate(variograms)
    plot!(p, v.empirical.abscissa, v.empirical.ordinate, marker=:o, color=i, label=scaling_classes[i])
    plot!(p, v.model, xlims=(0, 200), color=i, label="", )
end
p

params = map(zip(problems, variograms)) do (problem, variogram)
    get_lungs_params(problem, variogram.model)
end

sol = solve(problems[1], LUGS(:nasc => (variogram = variograms[1].model,)))

dists = [Gamma, InverseGamma, InverseGaussian, LogNormal]
dist_fits = compare_distributions(dists, acoustics.nasc, params[4], nreplicates=100, verbose=true)

optimal_dists = [choose_distribution(dists, acoustics.nasc, p) for p in params]
zdists = [get_zdists(d, p) for (d, p) in zip(optimal_dists, params)]

sim_fields = [nonneg_lusim(p, z) for (p, z) in zip(params, zdists)]
sim_plots = map(enumerate(sim_fields)) do (i, x)
    plot(domain(sol), zcolor=x, clims=(0, quantile(x, 0.999)), 
        markerstrokewidth=0, markershape=:square, title=string(scaling_classes[i]),
        markersize=1.5, xlabel="Easting (km)", ylabel="Northing (km)")
    df = @subset(acoustics, :class .== scaling_classes[i])
    scatter!(df.x, df.y, color=:white, markersize=df.nasc*1e-3, alpha=0.3,
        markerstrokewidth=0)
end
plot(sim_plots..., size=(1000, 1000))

i = 1
x = nonneg_lusim(params[i], zdists[i])
plot(domain(sol), zcolor=x, markerstrokewidth=0, markershape=:square, clims=(0, quantile(x, 0.999)), 
    title=string(scaling_classes[i]), markersize=4.2, background_color=:black, 
    xlabel="Easting (km)", ylabel="Northing (km)", size=(1000, 1000))

cal_error = Normal(0, 0.1)
dA = (resolution * km2nmi)^2
nreplicates = 30
class_results = map(scaling_classes) do class
    # class = "SS4"
    println(class)
    println("Setting up LUNGS...")
    acoustics_sub = @subset(acoustics, :class .== class)
    scaling_sub = resample_scaling(@subset(scaling, :class .== class))
    variogram, problem = define_conditional_sim(acoustics_sub)
    params = get_lungs_params(problem, variogram.model)
    optimal_dist = choose_distribution(dists, acoustics_sub.nasc, params)
    zdists = get_zdists(optimal_dist, params)

    println("Bootstrapping...")
    results = @showprogress map(1:nreplicates) do i
        trawl_means = get_trawl_means(scaling_sub, trawl_locations, true)
        popat!(trawl_means, rand(1:nrow(trawl_means)))
        geotrawl_means = @chain trawl_means begin
            @select(:x, :y, :ts, :length, :weight) 
            georef((:x, :y))
        end
        age_comp = proportion_at_age(scaling_sub, stochastic=true)
        ii = trawl_assignments_rand(coordinates.(surveydomain), 
                    coordinates.(domain(geotrawl_means)))

        df = DataFrame(nasc = nonneg_lusim(params, zdists) .* exp10(rand(cal_error) / 10),
                class = class,
                event_id = trawl_means.event_id[ii])
        df = @chain df begin
            leftjoin(@select(trawl_means, :event_id, :sigma_bs, :weight), on=:event_id)
            leftjoin(age_comp, on=:event_id)
            DataFrames.stack(Not([:event_id, :class, :nasc, :sigma_bs, :weight]), variable_name=:age, value_name=:p_age)
            DataFramesMeta.@transform(:n_age = :nasc ./ (4π * :sigma_bs) .* :p_age)
            @by(:age, 
                :n_age = sum(skipmissing(:n_age)) * dA,
                :weight = mean(:weight))
            DataFramesMeta.@transform(:i = i)
        end
        return df
    end
    return vcat(results...)
end

results = @chain vcat(class_results...) begin
    @by([:age, :i],
        :n_age = sum(:n_age), :weight = mean(:weight))
    DataFramesMeta.@transform(:biomass_age = :n_age .* :weight)
end



# nreplicates = 1000
# results = @showprogress map(1:nreplicates) do i
#     trawl_means = get_trawl_means(scaling, trawl_locations, true)
#     popat!(trawl_means, rand(1:nrow(trawl_means)))
#     geotrawl_means = @chain trawl_means begin
#         @select(:x, :y, :ts, :length, :weight) 
#         georef((:x, :y))
#     end
#     age_comp = proportion_at_age(scaling, stochastic=true)
#     dfs = map(enumerate(scaling_classes)) do (i, class)
#         df = DataFrame(nasc = nonneg_lusim(params[i], zdists[i]),
#             class = class,
#             event_id = trawl_assignments_rand(coordinates.(surveydomain), 
#                 coordinates.(domain(geotrawl_means))))
#     end
#     df = vcat(dfs...)
#     df = @chain df begin
#         leftjoin(@select(trawl_means, :event_id, :class, :sigma_bs), on=[:event_id, :class])
#         leftjoin(age_comp, on=[:event_id, :class])
#         DataFrames.stack(Not([:event_id, :class, :nasc, :sigma_bs]), variable_name=:age, value_name=:p_age)
#         DataFramesMeta.@transform(:n_age = :nasc ./ (4π * :sigma_bs) .* :p_age)
#         @by(:age, :n_age = sum(skipmissing(:n_age)) * dA)
#         DataFramesMeta.@transform(:i = i)
#     end
#     return df
# end
# results = vcat(results...)

@df results density(:n_age/1e9, group=:age, xlims=(0, 4),
    fill=true, alpha=0.7, ylims=(0, 30), palette=:Paired_10,
    xlabel="Billions of fish", ylabel="Probability density")
@df results density(:n_age/1e9, group=:age, palette=:Paired_10, fill=true, 
    alpha=0.7, ylims=(0, 0.01),
    xlabel="Billions of fish", ylabel="Probability density")

@df results boxplot(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    ylabel="Billions of fish")
@df results boxplot(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    ylims=(0, 10), ylabel="Billions of fish")

@by(results, :age, 
    :n_age = mean(:n_age),
    :std_age = std(:n_age), 
    :cv_age = std(:n_age) / mean(:n_age) * 100)


@df results density(:biomass_age/1e9, group=:age, xlims=(0, 4),
    fill=true, alpha=0.7, ylims=(0, 30), palette=:Paired_10,
    xlabel="Million tons", ylabel="Probability density")
@df results density(:n_age/1e9, group=:age, palette=:Paired_10, fill=true, 
    alpha=0.7, ylims=(0, 0.01),
    xlabel="Billions of fish", ylabel="Probability density")

@df results boxplot(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    ylabel="Billions of fish")
@df results boxplot(:age, :n_age/1e9, group=:age, palette=:Paired_10,
    ylims=(0, 10), ylabel="Billions of fish")

@by(results, :age, 
    :n_age = mean(:n_age),
    :std_age = std(:n_age), 
    :cv_age = std(:n_age) / mean(:n_age) * 100)
