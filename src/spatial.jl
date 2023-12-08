# const zdist_candidates = [Gamma, InverseGamma, InverseGaussian, LogNormal]

"""
    define_conditional_sim(acoustics, sim_domain[; maxlag=200.0, nlags=10,
        weigtfunc=h -> 1/h])

Set up a conditional geostatistical simulation, to probabilistically interpolate the
observed NASC data in `acoustics` to the simulation points in `sim_domain`.

An empirical variogram with `nlags` bins spaced evenly between 0 and `maxlag` is calculated
for the supplied acoustic data. An exponential variogram model is then fitted to this
empirical variogram via weighted least squares. By default the weighting function is 
inverse to the lag distance.

Returns a tuple with two elements:
- `variogram` : `NamedTuple` containing the empirical and model variograms in fields 
    `empirical` and `model`
- `geoproblem` : `GeoStats.SimulationProblem` containing the specifications for the
    simulation.
"""
function define_conditional_sim(acoustics, sim_domain; maxlag=200.0, 
        nlags=10, weightfunc=h -> 1/h)
    geonasc = acoustics[!, [:nasc, :x, :y]]
    geonasc.nasc .+= 1e-3 # add a small epsilon so no zeros
    geonasc.x .+= 1e-3 .* randn.()
    geonasc = georef(geonasc, (:x, :y))
    evg = EmpiricalVariogram(geonasc, :nasc, nlags=nlags, maxlag=maxlag)
    tvg = fit(ExponentialVariogram, evg, weightfunc)
    geoproblem = SimulationProblem(geonasc, sim_domain, :nasc => Float64, 1)
    variogram = (empirical = evg, model = tvg)
    return (variogram, geoproblem)
end

"""
    get_lungs_params(geoproblem, variogram[, variable=:nasc])

Calculate the parameters for lower-upper non-gaussian simulation from a
`GeoStats.SimulationProblem` and a `GeoStats.Variogram` model. If the variable to be 
simulated is something other than `:nasc`, this can be specified with the optional
`variable` argument. Returns a `NamedTuple` with the following fields:

- `data` : Observed data for the conditional simulation
- `μx` : Mean value at each simulation point
- `L` : Lower-triangular Cholesky factor of the covariance matrix of simulated data points
- `μ` : 
- `dlocs`, `slocs` : Locations of data and simulation points.
"""
function get_lungs_params(geoproblem, variogram, variable=:nasc)
    solver = LUGS(variable => (variogram = variogram,))
    preproc = preprocess(geoproblem, solver);
    pars = preproc[Set([variable])][Set([variable])]
    return (data=pars[1], μx=pars[2], L=pars[3], μ=pars[4], dlocs=pars[5], slocs=pars[6])
end

dist_params(d::Type{Gamma}, μ, v) = (v / μ, μ^2 / v)
dist_params(d::Type{InverseGaussian}, μ, v) = (μ, μ^3 / v)
dist_params(d::Type{InverseGamma}, μ, v) = (μ^2 / v + 2, μ^3/v + μ)
dist_params(d::Type{LogNormal}, μ, v) = ( log(μ) - log(v/exp(2log(μ)) + 1) / 2, sqrt(log(v/exp(2log(μ)) + 1)) )

"""
    parameterize_zdists(Dist, lungs_params[, ϵ=cbrt(eps())])

Calculate the parameters of the white noise distributions driving a lower-upper
nonnegative Gaussian simulation. Returns a vector of parameterized distributions from the 
family specified by `Dist`, each with variance==1 and a mean that satisfies the
requirements of the conditional simulation.

# Arguments
- `Dist` : A non-negative continuous `Distribution` type from Distributions.jl. Currently 
    `Gamma`, `InverseGaussian`, `InverseGamma`, and `LogNormal` are supported.
- `lungs_params` : Tuple of parameters describing the geostatistical conditional problem,
    including the vector of desired mean values at the simulation points `μx` and the 
    Cholesky triangle of their covariance matrix `L`.  These are obtained via
    `get_lungs_params`.
- `ϵ=cbrt(eps())` : Small value to add to zero-valued elements of the mean vectors. 
    Ensures that each z-distribution has a nonzero mean/variance. 
"""
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
    x[sim_locs] = L * z
    x[data_locs] = data
    return x
end

function nonneg_lusim(params, zdists)
    z = rand.(zdists)
    return nonneg_lumult(params, z)
end

"""
    nonneg_lusim(scp::ScalingClassProblem)

Generate a single conditional simulation of NASC based on the parameters defined in `scp`, 
an instance of `ScalingClassProblem`.

This function is an alias of `nonneg_lusim(scp::ScalingClassProblem)`
"""
function nonneg_lusim(scp::ScalingClassProblem)
    z = rand.(scp.zdists)
    return nonneg_lumult(scp.params, z)
end

"""
    simulate_nasc(scp::ScalingClassProblem)

Generate a single conditional simulation of NASC based on the parameters defined in `scp`, 
an instance of `ScalingClassProblem`.

This function is an alias of `nonneg_lusim(scp::ScalingClassProblem)`
"""
simulate_nasc(scp::ScalingClassProblem) = nonneg_lusim(scp)

"""
    choose_z_distribution(candidate_dists, nasc, lungs_params[; nreplicates=500, verbose=false])

Choose the distribution from `candidate_dists` that best approximates the distribution of
the observed data `nasc` when used to drive the conditional spatial simulation defined by
`lungs_params`.

Each of the candidates will be used to generate `nreplicates` conditional simulations 
(default number is 500). A histogram of the simulated values (with logarithmic bins) is 
calculated and compared to the histogram of the observed `nasc` via the Kullback-Liebler
divergence; the candidate distribution family with the lowest average KLD is returned.
"""
function choose_z_distribution(candidate_dists, nasc, lungs_params; nreplicates=500, verbose=false)
    bin_edges = [0; 2 .^ (0:14)]
    h_nasc = normalize(fit(StatsBase.Histogram, nasc, bin_edges), mode=:density)
    fit_list = []
    for Dist in candidate_dists
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
            # draw a random trawl, with probability inversely related to distance
            trawl_idx = sample(idx[i], Weights(dists[i].^-a))
        else
            # assign pixel i to the nearest trawl
            trawl_idx = idx[i][argmin(dists[i])]
        end
        assignments[i] = trawl_idx
    end
    return assignments
end
