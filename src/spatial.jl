# const zdist_candidates = [Gamma, InverseGamma, InverseGaussian, LogNormal]

"""
    define_conditional_sim(acoustics, surveydomain[; maxlag=200.0, nlags=10,
        weigtfunc=h -> 1/h])

Set up a conditional geostatistical simulation
"""
function define_conditional_sim(acoustics, surveydomain; maxlag=200.0, 
        nlags=10, weightfunc=h -> 1/h)
    geonasc = acoustics[!, [:nasc, :x, :y]]
    geonasc.nasc .+= 1e-3 # add a small epsilon so no zeros
    geonasc.x .+= 1e-3 .* randn.()
    geonasc = georef(geonasc, (:x, :y))
    evg = EmpiricalVariogram(geonasc, :nasc, nlags=nlags, maxlag=maxlag)
    tvg = fit(ExponentialVariogram, evg, weightfunc)
    prob = SimulationProblem(geonasc, surveydomain, :nasc => Float64, 1)
    variogram = (empirical = evg, model = tvg)
    return (variogram, prob)
end

"""
Parameters for lower-upper non-gaussian simulation
"""
function get_lungs_params(problem, variogram, variable=:nasc)
    solver = LUGS(variable => (variogram = variogram,))
    preproc = preprocess(problem, solver);
    params = preproc[Set([variable])][Set([variable])]
    return (data=params[1], μx=params[2], L=params[3], μ=params[4], dlocs=params[5], slocs=params[6])
end

dist_params(d::Type{Gamma}, μ, v) = (v / μ, μ^2 / v)
dist_params(d::Type{InverseGaussian}, μ, v) = (μ, μ^3 / v)
dist_params(d::Type{InverseGamma}, μ, v) = (μ^2 / v + 2, μ^3/v + μ)
dist_params(d::Type{LogNormal}, μ, v) = ( log(μ) - log(v/exp(2log(μ)) + 1) / 2, sqrt(log(v/exp(2log(μ)) + 1)) )

function parameterize_zdists(Dist, lungs_params, ϵ=cbrt(eps()))
    data, μx, L, μ, dlocs, slocs = lungs_params
    μx = copy(μx)
    μx[μx .< ϵ] .= ϵ 
    μx = μx .* mean(data) ./ mean(μx)
    μz = L \ μx
    μz[μz .<= ϵ] .= ϵ
    vz = ones(length(μz))
    return [Dist(p...) for p in dist_params.(Dist, μz, vz)]
end


function nonneg_lumult(params, z)
    data, μx, L, μ, data_locs, sim_locs = params
    npts = length(data_locs) + length(sim_locs)
    x = zeros(npts)
    # nonneg_lumult!(x, params, z)
    x[sim_locs] = L * z
    x[data_locs] = data
    return x
end

function nonneg_lusim(params, zdists)
    z = rand.(zdists)
    return nonneg_lumult(params, z)
end

function nonneg_lusim(scp::ScalingClassProblem)
    z = rand.(scp.zdists)
    return nonneg_lumult(scp.params, z)
end

simulate_nasc(scp::ScalingClassProblem) = nonneg_lusim(scp)


function choose_z_distribution(distributions, nasc, lungs_params; nreplicates=500, verbose=false)
    bin_edges = [0; 2 .^ (0:14)]
    h_nasc = normalize(fit(StatsBase.Histogram, nasc, bin_edges), mode=:density)
    fit_list = []
    for Dist in distributions
        if verbose
            println("Comparing with $(Dist)...")
        end
        zdists = parameterize_zdists(Dist, lungs_params) 
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
    return dist_fits.distribution[argmin(dist_fits.mean_kld)]    
end

# function avg_nearest_neighbor_distance(coords)
#     kdtree = KDTree(coords)
#     _, dist = knn(kdtree, coords, 2, true)
#     return mean(last.(dist))
# end

# function calculate_exponent(pixel_coords, trawl_coords, kdtree)
#     center_pixel = mean(pixel_coords)
#     d = avg_nearest_neighbor_distance(trawl_coords) / 2

#     idx, dist = knn(kdtree, center_pixel, length(trawl_coords), true)

# end


"""
    trawl_assignments(pixel_coords, trawl_coords[, stochastic=true, a=1.9])

Assign each acoustic cell to a trawl, either deterministically (i.e. a standard MACE 
nearest-trawl assignment) or probabilistically, weighted proportional to distance^-a.

# Arguments
- `pixel_coords`, `trawl_coords` : Vectors of coordinates for the acoustic cells and 
    trawl locations. Each element of these is a 2D vector with the x, y coordinates of the
    pixel or trawl.
- `stochastic=true` : Whether trawl assignment should be random (the default) or
    deterministic.
- `a=1.9` : Magic number 
"""
function trawl_assignments(pixel_coords, trawl_coords, stochastic=true, a=1.9)
    # Calculate a k-dimensional tree for efficiently finding nearest neighbors
    kdtree = KDTree(trawl_coords)
    #=
    idx and dists are vectors the same length as the number of pixels/acoustic cells.
    Each element of idx and dists is a vector the same length as the number of trawls.
    The ith element of idx is a vector of indices to each of the trawl locations.
    The ith element of dists is a vector of distances to the trawls indexed by idx.
    This means that the jth element of dists[i] is the distance from pixel i to trawl j.
    =#
    idx, dists = knn(kdtree, pixel_coords, length(trawl_coords))
    # allocate an empty vector for the trawl assignments
    assignments = Vector{Int}(undef, length(pixel_coords))
    
    for i in eachindex(assignments)
        if stochastic
            # draw a random trawl index, with probability inversely related to distance
            trawl_idx = sample(idx[i], Weights(dists[i].^-a))
        else
            # assign pixel i to the nearest trawl
            trawl_idx = idx[i][argmin(dists[i])]
        end
        assignments[i] = trawl_idx
    end
    # trawl_assignments!(assignments, pixel_coords, trawl_coords, stochastic)
    return assignments
end
