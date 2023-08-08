module ATBootstrap

using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using GeoStats, GeoStatsPlots
using Statistics, StatsBase
using Distributions
using LinearAlgebra
using Distances
using NearestNeighbors
using ProgressMeter

export read_survey_files,
    ATSurveyData,
    ATBootstrapProblem,
    resample_df,
    nonneg_lusim,
    nonneg_lusim!,
    nonneg_lumult,
    nonneg_lumult!,
    solution_domain,
    simulate_classes,
    stepwise_error,
    BootSpecs,
    error_labels

# struct SurveyData
#     acoustics
#     scaling
#     trawl_locations
# end

const zdist_candidates = [Gamma, InverseGamma, InverseGaussian, LogNormal]

function define_conditional_sim(acoustics, surveydomain, maxlag=200.0)
    geonasc = acoustics[!, [:nasc, :x, :y]]
    geonasc.nasc .+= 1e-3
    geonasc.x .+= 1e-3 .* randn.()
    geonasc = georef(geonasc, (:x, :y))
    evg = EmpiricalVariogram(geonasc, :nasc, nlags=20, maxlag=maxlag)
    tvg = fit(ExponentialVariogram, evg, h -> 1/h)
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

function nonneg_lumult!(x, params, z)
    data, μx, L, μ, dlocs, slocs = params
    x[slocs] = L * z
    x[dlocs] = data
end

function nonneg_lumult(params, z)
    data, μx, L, μ, dlocs, slocs = params
    npts = length(dlocs) + length(slocs)
    x = zeros(npts)
    nonneg_lumult!(x, params, z)
    return x
end

function nonneg_lusim(params, zdists)
    z = rand.(zdists)
    return nonneg_lumult(params, z)
end

function nonneg_lusim(atbp)
    z = rand.(atbp.zdists)
    return nonneg_lumult(atbp.params, z)
end


function compare_distributions(distributions, nasc, lungs_params; nreplicates=500, verbose=false)
    bin_edges = [0; 2 .^ (0:14)]
    h_nasc = normalize(fit(StatsBase.Histogram, nasc, bin_edges), mode=:density)
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

include("mace_ts.jl")

function make_weight_function(stochastic=false, sd=0.05)
    err = stochastic ? (1 + sd*randn()) : 1.0
    predict_weight(L) = (1e-5 * L^2.9) * err
    return predict_weight
end

# Numbers from Matta and Kimura 2012, Age determination manual
const L∞ = 67.33 # mm
const t₀ = -0.205
const K = 0.1937
const AGE_MAX = 10
predict_length(t) = L∞ * (1 - exp(-K * (t - t₀)))
function predict_age_deterministic(L, age_max=AGE_MAX)
    if L < predict_length(age_max)
        return round(Int, log(1 - L/L∞) / -K + t₀)
    else
       return age_max
   end
end

function predict_age_stochastic(L, age_max=AGE_MAX)
    return max(0, predict_age_deterministic(L + 2randn(), age_max))# + rand([-1, 0, 0, 0, 1]))
end

function predict_age(L, stochastic=true, age_max=AGE_MAX)
    if stochastic
        return predict_age_stochastic(L, age_max)
    else
        return predict_age_deterministic(L, age_max)
    end
end


const a = 1.9
function trawl_assignments!(assignments, kdtree::KDTree, pixel_coords, trawl_coords, stochastic)
    idx, dists = knn(kdtree, pixel_coords, length(trawl_coords))
    if stochastic
        idx1 = [sample(i, Weights(1 ./ d.^a)) for (i, d) in zip(idx, dists)]
    else
        idx1 = [i[argmin(d)] for (i, d) in zip(idx, dists)]
    end
    assignments .= idx1
end

function trawl_assignments!(assignments, pixel_coords, trawl_coords, stochastic=true)
    kdtree = KDTree(trawl_coords)
    trawl_assignments!(assignments, kdtree, pixel_coords, trawl_coords, stochastic)
end

function trawl_assignments(pixel_coords, trawl_coords, stochastic=true)
    assignments = Vector{Int}(undef, length(pixel_coords))
    trawl_assignments!(assignments, pixel_coords, trawl_coords, stochastic)
    return assignments
end

function resample_df(df, stochastic=true)
    n = nrow(df)
    if stochastic
        ii = sample(1:n, n)
    else
        ii = 1:n
    end
    return @view df[ii, :]
end

function resample_scaling(df, stochastic=true)
    return DataFramesMeta.combine(x -> resample_df(x, stochastic),
        DataFramesMeta.groupby(df, [:event_id, :class]))
end

function get_trawl_means(scaling, trawl_locations)
    trawl_means = @chain scaling begin
        @by(:event_id, 
            :sigma_bs = mean(:sigma_bs, Weights(:w)),
            :sigma_bs_std = std(:sigma_bs, Weights(:w)),
            :length = mean(:primary_length, Weights(:w)))        
            DataFramesMeta.@transform(:ts = 10log10.(:sigma_bs))
        innerjoin(trawl_locations, on=:event_id)
    end
    return trawl_means
end

function make_all_ages(scaling, age_max)
    return allcombinations(DataFrame, 
        event_id = unique(scaling.event_id), 
        age = 0:age_max)
end

function weights_at_age(scaling, all_ages, stochastic=false)
    predict_weight = make_weight_function(stochastic)
    res = @chain scaling begin
        DataFramesMeta.@transform(
            # :age = predict_age.(:primary_length, stochastic), 
            :weight = predict_weight.(:primary_length))
        rightjoin(all_ages, on=[:event_id, :age])
        DataFramesMeta.@transform(:weight = replace(:weight, missing => 0.0))
        DataFramesMeta.@transform(:age = lpad.(string.(:age), 2, "0"))
        @by(:age, :weight = mean(:weight))
    end
    return res
end

function proportion_at_age(scaling, all_ages; age_max=AGE_MAX)#, stochastic=false)
    age_comp = @chain scaling begin
        # DataFramesMeta.@transform(:age = predict_age.(:primary_length, stochastic))
        @by([:event_id, :age], :p_age=sum(:w))
        DataFrames.groupby(:event_id)
        DataFramesMeta.@transform(:p_age = :p_age / sum(:p_age))
        rightjoin(all_ages, on=[:event_id, :age])
        DataFramesMeta.@transform(:p_age = replace(:p_age, missing => 0.0))
        DataFramesMeta.@transform(:age = lpad.(string.(:age), 2, "0"))
        @orderby(:event_id, :age)
    end
    # return age_comp
    return DataFrames.unstack(age_comp, :age, :p_age)
end

function simulate_cal_error(cal_error, stochastic=true) 
    if stochastic
        return exp10(cal_error*randn() / 10)
    else
        return exp10(0.0)
    end
end

function ATSurveyData(acoustics, scaling, trawl_locations, domain)
    return (;acoustics, scaling, trawl_locations, domain)
end

@kwdef struct BootSpecs
    predict_ts::Bool=true
    resample_scaling::Bool=true
    jackknife_trawl::Bool=true
    age_length::Bool=true
    weights_at_age::Bool=true
    trawl_assignments::Bool=true
    nonneg_lusim::Bool=true
    calibration::Bool=true
end
BootSpecs(b::Bool) = BootSpecs(fill(b, length(fieldnames(BootSpecs)))...)

function ATBootstrapProblem(surveydata, class, cal_error, dA;
        zdist_candidates=zdist_candidates)
    acoustics_sub = @subset(surveydata.acoustics, :class .== class)
    variogram, problem = define_conditional_sim(acoustics_sub, surveydata.domain)
    params = get_lungs_params(problem, variogram.model)
    optimal_dist = choose_distribution(zdist_candidates, acoustics_sub.nasc, params)
    zdists = get_zdists(optimal_dist, params)
    return (; class, variogram, problem, params, optimal_dist, zdists,
        cal_error, dA)
end

function solution_domain(atbp, variable=:nasc)
    sol = solve(atbp.problem, LUGS(variable => (variogram = atbp.variogram.model,)))
    return domain(sol)
end

function simulate(atbp, surveydata; nreplicates=500, bs=BootSpecs(), age_max=AGE_MAX)
    acoustics, scaling, trawl_locations, scaling_classes = surveydata
    class, variogram, problem, params, optimal_dist, zdists, cal_error, dA = atbp
    scaling_sub = @subset(scaling, :class .== class)

    z0 = mean.(zdists)

    println("Bootstrapping $(class)...")
    results = @showprogress map(1:nreplicates) do i
        scaling_boot = resample_scaling(scaling_sub, bs.resample_scaling)

        predict_ts = make_ts_function(bs.predict_ts)
        scaling_boot = DataFramesMeta.@transform(scaling_boot,
            :sigma_bs = exp10.(predict_ts.(:ts_relationship, :ts_length)/10),
            :age = predict_age.(:primary_length, bs.age_length))

        trawl_means = get_trawl_means(scaling_boot, trawl_locations)
        if bs.jackknife_trawl
            popat!(trawl_means, rand(1:nrow(trawl_means)))
        end
        geotrawl_means = @chain trawl_means begin
            @select(:x, :y, :ts, :length) 
            georef((:x, :y))
        end
        all_ages = make_all_ages(scaling_boot, age_max)
        age_comp = proportion_at_age(scaling_boot, all_ages)
        age_weights = weights_at_age(scaling_boot, all_ages, bs.weights_at_age)
        @assert ! any(ismissing, age_weights.weight)
        if nrow(age_weights) != AGE_MAX+1 
            println("nrow=$(nrow(age_weights))")
        end
        ii = trawl_assignments(coordinates.(surveydata.domain), 
                    coordinates.(domain(geotrawl_means)), bs.trawl_assignments)

        nasc = bs.nonneg_lusim ? nonneg_lusim(atbp) : nonneg_lumult(params, z0)
        cal_error_sim = simulate_cal_error(cal_error,  bs.calibration)

        df = DataFrame(
            nasc = nasc * cal_error_sim,
            class = class,
            event_id = trawl_means.event_id[ii]
        )
        df = @chain df begin
            leftjoin(@select(trawl_means, :event_id, :sigma_bs), on=:event_id)
            leftjoin(age_comp, on=:event_id)
            DataFrames.stack(Not([:event_id, :class, :nasc, :sigma_bs]), variable_name=:age, value_name=:p_age)
            DataFramesMeta.@transform(:n_age = :nasc ./ (4π * :sigma_bs) .* :p_age)
            @by(:age, 
                :n_age = sum(skipmissing(:n_age)) * dA)
            leftjoin(age_weights, on=:age)
            DataFramesMeta.@transform(:biomass_age = :n_age .* :weight, :i = i)
        end
        return df
    end

    return vcat(results...)
end

function simulate_classes(class_problems, surveydata; nreplicates=500, bs=BootSpecs())
    class_results = map(p -> simulate(p, surveydata; nreplicates, bs), class_problems)
    results = @chain vcat(class_results...) begin
        @by([:age, :i],
            :n_age = sum(:n_age), :biomass_age = sum(:biomass_age))
    end
    return results
end

function stepwise_error(class_problems, surveydata; remove=true, kwargs...)
    error_sources = string.(fieldnames(BootSpecs))
    colname = remove ? :eliminated_error : :added_error
    results = map(eachindex(error_sources)) do i
        prefix = remove ? "Omitting" : "Adding"
        println("\n$(prefix) $(error_sources[i]) ($(i)/$(length(error_sources)))...")
        errs = fill(remove, length(error_sources))
        errs[i] = !remove
        bs = BootSpecs(errs...)
        res = simulate_classes(class_problems, surveydata; bs, kwargs...)
        res[!, colname] .= error_sources[i]
        res
    end
    return vcat(results...)
end

# Ordered, categorical labels for stepwise error results
error_labels = DataFrame(
    added_error = ["calibration", "nonneg_lusim", "jackknife_trawl", "trawl_assignments",
        "resample_scaling", "weights_at_age", "predict_ts", "age_length", "All"],
    error_label = CategoricalArray(
        ["Calibration", "Spatial sampling", "Trawl jackknife", "Trawl assignment", 
        "Resample catches", "Length-weight", "TS models", "Age-length", "All"],
        levels=["Calibration", "Spatial sampling", "Trawl jackknife", "Trawl assignment", 
        "Resample catches", "Length-weight", "TS models", "Age-length", "All"]
    )
)

function read_survey_files(surveydir)
    acoustics = CSV.read(joinpath(surveydir, "acoustics_projected.csv"), DataFrame)
    trawl_locations = CSV.read(joinpath(surveydir, "trawl_locations_projected.csv"), DataFrame)
    scaling = CSV.read(joinpath(surveydir, "scaling.csv"), DataFrame)
    surveydomain = CSV.read(joinpath(surveydir, "surveydomain.csv"), DataFrame)
    surveydomain = DataFrames.shuffle(surveydomain) # this seems to fix the issue with directional artifacts
    surveydomain =  PointSet(Matrix(surveydomain)')
    return (;acoustics, scaling, trawl_locations, surveydomain)
end



end # module