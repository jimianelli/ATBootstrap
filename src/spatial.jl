const zdist_candidates = [Gamma, InverseGamma, InverseGaussian, LogNormal]

function define_conditional_sim(acoustics, surveydomain; maxlag=200.0, 
        nlags=20, weightfunc=h -> 1/h)
    geonasc = acoustics[!, [:nasc, :x, :y]]
    geonasc.nasc .+= 1e-3
    geonasc.x .+= 1e-3 .* randn.()
    geonasc = georef(geonasc, (:x, :y))
    evg = EmpiricalVariogram(geonasc, :nasc, nlags=nlags, maxlag=maxlag)
    tvg = fit(ExponentialVariogram, evg, weightfunc)
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
