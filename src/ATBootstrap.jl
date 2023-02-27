# module ATBootstrap
using CSV, DataFrames, DataFramesMeta
using GeoStats, GeoStatsPlots
using Statistics, StatsBase
using Distributions
using LinearAlgebra

# struct SurveyData
#     acoustics
#     scaling
#     trawl_locations
# end



function define_conditional_sim(acoustics, maxlag=200.0)
    geonasc = acoustics[!, [:nasc, :x, :y]]
    geonasc.nasc .+= 1e-3
    geonasc.x .+= 1e-3 .* randn.()
    geonasc = georef(geonasc, (:x, :y))
    evg = EmpiricalVariogram(geonasc, :nasc, maxlag=maxlag)
    tvg = fit(ExponentialVariogram, evg)
    prob = SimulationProblem(geonasc, surveydomain, :nasc => Float64, 1)
    variogram = (empirical = evg, model = tvg)
    return (variogram, prob)
end

# struct LUNGS
#     vparams
# end

# function GeoStats.solve(problem, solver::LUNGS)

# end


"""
Parameters for lower-upper non-gaussian simulation
"""
function get_lungs_params(problem, variogram, variable=:nasc)
    solver = LUGS(variable => (variogram = variogram,))
    # sol = solve(problem, solver)
    preproc = preprocess(problem, solver);
    params = preproc[(variable,)][1]
    return (data=params[1], μx=params[2], L=params[3], μ=params[4], dlocs=params[5], slocs=params[6])
end



dist_params(d::Type{Gamma}, μ, v) = (v / μ, μ^2 / v)
dist_params(d::Type{InverseGaussian}, μ, v) = (μ, μ^3 / v)
dist_params(d::Type{InverseGamma}, μ, v) = (μ^2 / v + 2, μ^3/v + μ)
dist_params(d::Type{LogNormal}, μ, v) = ( log(μ) - log(v/exp(2log(μ)) + 1) / 2, sqrt(log(v/exp(2log(μ)) + 1)) )

function get_zdists(Dist, lungs_params, ϵ=cbrt(eps()))
    data, μx, L, μ, dlocs, slocs = lungs_params
    npts = length(dlocs) + length(slocs)
    μx = copy(μx)
    μx[μx .< ϵ] .= ϵ 
    μx = μx .* mean(data) ./ mean(μx)
    μz = L \ μx
    μz[μz .<= ϵ] .= ϵ
    vz = ones(length(μz))
    return [Dist(p...) for p in dist_params.(Dist, μz, vz)]
end

function nonneg_lusim!(x, params, zdists)
    data, μx, L, μ, dlocs, slocs = params
    x[slocs] = L * rand.(zdists)
    x[dlocs] = data
end

function nonneg_lusim(params, zdists)
    data, μx, L, μ, dlocs, slocs = params
    npts = length(dlocs) + length(slocs)
    x = zeros(npts)
    nonneg_lusim!(x, params, zdists)
    return x
end


function compare_distributions(distributions, nasc, lungs_params; nreplicates=500, verbose=false)
    bin_edges = [0; 2 .^ (0:14)]
    h_nasc = normalize(fit(StatsBase.Histogram, acoustics.nasc, bin_edges), mode=:density)
    fit_list = []
    for Dist in distributions
        if verbose
            println("Comparing with $(Dist)...")
        end
        zdists = get_zdists(Dist, lungs_params) 
        fits = map(1:nreplicates) do i
            x = nonneg_lusim(lungs_params, zdists)
            h_sim = normalize(fit(StatsBase.Histogram, x, bin_edges), mode=:density)
            kld = evaluate(KLDivergence(), h_nasc.weights, h_sim.weights)
            return (distribution = Dist, kld = kld)
        end
        push!(fit_list, DataFrame(fits))
    end
    dist_fits = @chain vcat(fit_list...) begin
        @subset(isfinite.(:kld))
        @by(:distribution, 
            :mean_kld = mean(:kld), 
            :se_kld = std(:kld) / sqrt(length(:kld)))
    end
    return dist_fits
end


function choose_distribution(distributions, nasc, lungs_params; nreplicates=500, verbose=false)
    dist_fits = compare_distributions(distributions, nasc, lungs_params, nreplicates=nreplicates,
        verbose = verbose)
    return dist_fits.distribution[argmin(dist_fits.mean_kld)]    
end

predict_ts(L) = -66 + 20log10(L)
predict_ts_stochastic(L) = -66 + 0.14randn() + 20log10(L) # standard error from Nate's paper (email on 2023-01-31)
predict_weight(L) = 1e-5 * L^2.9
predict_weight_stochastic(L) = predict_weight(L) + 


# Numbers from Matta and Kimura 2012, Age determination manual
const L∞ = 67.33 # mm
const t₀ = -0.205
const K = 0.1937
const AGE_MAX = 10
predict_length(t) = L∞ * (1 - exp(-K * (t - t₀)))
function predict_age(L, age_max=AGE_MAX)
    if L < predict_length(age_max)
        return round(Int, log(1 - L/L∞) / -K + t₀)
    else
       return age_max
   end
end
predict_age_stochastic(L, age_max=AGE_MAX) = max(1, predict_age(L, age_max) + rand([-1, 0, 0, 0, 1]))

const a = 1.9
function trawl_assignments_rand!(assignments, kdtree, pixel_coords, trawl_coords)
    idx, dists = knn(kdtree, pixel_coords, length(trawl_coords))
    idx1 = [sample(i, Weights(1 ./ d.^a)) for (i, d) in zip(idx, dists)]
    assignments .= idx1
end
function trawl_assignments_rand!(assignments, pixel_coords, trawl_coords)
    kdtree = KDTree(trawl_coords)
    trawl_assignments_rand!(assignments, kdtree, pixel_coords, trawl_coords)
end
function trawl_assignments_rand(pixel_coords, trawl_coords)
    assignments = Vector{Int}(undef, length(pixel_coords))
    trawl_assignments_rand!(assignments, pixel_coords, trawl_coords)
    return assignments
end

function resample_df(df)
    n = nrow(df)
    ii = sample(1:n, n)
    return @view df[ii, :]
end

function resample_scaling(df)
    return DataFramesMeta.combine(resample_df, DataFramesMeta.groupby(df, [:event_id, :class]))
end

function get_trawl_means(scaling, trawl_locations, stochastic=false)
    ts = stochastic ? predict_ts_stochastic : predict_ts
    weight = stochastic ? predict_weight_stochastic : predict_weight
    trawl_means = @chain scaling begin
        DataFramesMeta.@transform(:sigma_bs = exp10.(ts.(:ts_length)/10),
                                  :weight = weight.(:primary_length))
        @by(:event_id, 
            :sigma_bs = mean(:sigma_bs, Weights(:w)),
            :sigma_bs_std = std(:sigma_bs, Weights(:w)),
            :length = mean(:primary_length, Weights(:w)),
            :weight = mean(:weight, Weights(:w)))
        DataFramesMeta.@transform(:ts = 10log10.(:sigma_bs))
        innerjoin(trawl_locations, on=:event_id)
    end
    return trawl_means
end

function proportion_at_age(scaling; age_max=AGE_MAX, stochastic=false)
    age = stochastic ? predict_age_stochastic : predict_age
    all_ages = @chain Iterators.product(unique(scaling.event_id), 1:age_max) begin
        DataFrame()
        rename([:event_id, :age])
        sort()
    end
    age_comp = @chain scaling begin
        DataFramesMeta.@transform(:age = age.(:primary_length))
        @by([:event_id, :age], :p_age=sum(:w))
        DataFrames.groupby(:event_id)
        DataFramesMeta.@transform(:p_age = :p_age / sum(:p_age))
        rightjoin(all_ages, on=[:event_id, :age])
        DataFramesMeta.@transform(:p_age = replace(:p_age, missing => 0.0))
        DataFramesMeta.@transform(:age = "age" .* lpad.(string.(:age), 2, "0"))
        @orderby(:event_id, :age)
    end
    # return age_comp
    return DataFrames.unstack(age_comp, :age, :p_age)
end

# end # module