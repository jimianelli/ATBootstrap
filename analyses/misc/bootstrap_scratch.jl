using StatsPlots
using CSV, DataFrames, DataFramesMeta
using Geodesy
using GeoStats
using Statistics
using Distributions
using Random
using LinearAlgebra

intervals = @chain CSV.read(joinpath(@__DIR__, "intervals.csv"), DataFrame) begin
    @subset(:survey .==201807)
    @select(:transect, :interval, :lat = :start_latitude, :lon = :start_longitude)
end
scatter(intervals.lon, intervals.lat,
    markerstrokewidth=0, markersize=3)


acoustics = CSV.read(joinpath(@__DIR__, "acoustics_2018.csv"), DataFrame)
acoustics = @chain acoustics begin
    @subset(:class .== "PK1")
    @by(:interval, :nasc = sum(:nasc))
end
acoustics = leftjoin(intervals, acoustics, on=:interval)
acoustics.nasc[ismissing.(acoustics.nasc)] .= 0

scatter(acoustics.lon, acoustics.lat,
    markerstrokewidth=0, markersize=1.5)


utmzone = 3
lla = LLA.(acoustics.lat, acoustics.lon, 0.0)
utm = [UTM(x, utmzone, true, wgs84) for x in lla]
acoustics.x = [u.x / 1e3 for u in utm]
acoustics.y = [u.y / 1e3 for u in utm]



geonasc = acoustics[!, [:nasc, :x, :y]]
geonasc.nasc .+= 1e-3
geonasc = georef(geonasc, (:x, :y))
plot(geonasc)

evg = EmpiricalVariogram(geonasc, :nasc, maxlag=200)
plot(evg)
tvg = fit(ExponentialVariogram, evg)
plot!(tvg)
sill(tvg)
nugget(tvg)
range(tvg)
mean(acoustics.nasc)
quantile(acoustics.nasc, [0.1, 0.5, 0.9])


prob = SimulationProblem(domain(geonasc), :nasc => Float64, 1)
S1 = LUGS(:nasc => (variogram = tvg, mean=mean(acoustics.nasc)))
sol = solve(prob, S1)
plot(sol)
preproc = preprocess(prob, S1)

params = preproc[(:nasc,)]
d = params[1][2]
L = params[1][3]


# w = randn(length(acoustics.nasc))
w = rand(Gamma(0.006, 12.0), length(acoustics.nasc))
y = d .+ L * w
plot(acoustics.nasc)
plot!(y)

mean(y)
mean(acoustics.nasc)

density(y)
density!(acoustics.nasc)

geosim = georef(DataFrame(x = acoustics.x, y=acoustics.y, nasc=y), (:x, :y))
evgsim = EmpiricalVariogram(geosim, :nasc, maxlag=200)
plot(evgsim)
plot!(evg)



###
using ConcaveHull
transect_ends = @chain acoustics begin
    @orderby(:y)
    @by(:transect, :x = [first(:x), last(:x)], :y = [first(:y), last(:y)])
end

plot(geosim, markersize=2, title="", xlabel="Easting (km)", ylabel="Northing (km)")
v = [[row.x, row.y] for row in eachrow(transect_ends)]#coordinates.(domain(pollock0))
surveyhull = concave_hull(v, 6)
plot!(surveyhull)
dx = 5.0
dy = 5.0
xgrid = range(round.(extrema(transect_ends.x))..., step=dx)
ygrid = range(round.(extrema(transect_ends.y))..., step=dy)
surveydomain = PointSet(hcat([[x, y] for x in xgrid, y in ygrid
    if in_hull([x, y], surveyhull)]...))
plot(surveydomain, markersize=1.2, xlabel="Easting (km)", ylabel="Northing (km)")



prob2 = SimulationProblem(geonasc, surveydomain, :nasc => Float64, 1)
sol2 = solve(prob2, S1)
preproc2 = preprocess(prob2, S1)

params2 = preproc2[(:nasc,)]
z2, d2, L2, μ, dlocs2, slocs2 = params2[1]
npts2 = length(dlocs2) + length(slocs2)

# Q₁₁ = sill(tvg) .- pairwise(tvg, domain(geonasc))
# Q₁₂ = sill(tvg) .- pairwise(tvg, domain(geonasc), surveydomain)
# Q₂₂ = sill(tvg) .- pairwise(tvg, surveydomain)

# L₁₁ = cholesky(Symmetric(Q₁₁)).L
# B₁₂ = L₁₁ \ Q₁₂
# A₂₁ = B₁₂'

# d₂ = A₂₁ * (L₁₁ \ acoustics.nasc)
# L₂₂ = cholesky(Symmetric(Q₂₂ - A₂₁*B₁₂)).L

# plot(d2)
# s = diag(L2)
# L = L2 * spdiagm(1 ./ s)
# D = spdiagm(s.^2)


Q = sill(tvg) .- pairwise(tvg, domain(sol2))
# Q11 = C[dlocs2, dlocs2]
# Q12 = C[dlocs2, slocs2]
Q22 = C[slocs2, slocs2]

plot(d2)
mean(d2)
sum(d2 .<= 0)
d21 = copy(d2); d21[d21 .<= 0] .= 1e-2
plot!(d21)
mean(d21)
plot!(sqrt.(diag(Q22)))
d21 .*= mean(acoustics.nasc) / mean(d21)
plot!(d21)
# μw = vec(nonneg_lsq(L21, d21)) .+= 1e-2
μw = L2 \ d21
μw[μw .<= 0] .= 1e-6

# plot(plot(μw), plot(L2 \ d21), layout=(1, 2), size=(1000, 500))


plot(d21, label="d")
plot!(L2 * μw, label="reconstructed")

plot(μw)
# vw = fill(1., length(μw))
# C22 = L2 * L2'
# below calculateds diagonal of C22 = L2 * L2'
dC22 = [sum(row.^2) for row in eachrow(L2)]
plot(dC22)

vw = ones(length(dC22))#(L2.^2 \ dC22) # (mean(d2) / mean(acoustics.nasc))^2
plot(vw)
plot!(L2.^2 \ dC22, ylims=(0, 2))

shape = vw ./ μw
scale = μw.^2 ./ vw

# Could define some convenience functions here to get parameters,
# dispatching on each distribution, given desired mean and variance, e.g.:
# shape, scale = params_from_mean_var(Gamma, mean, var)
# and then do an exhastive search through them

# zdists = Exponential.(μw)
# zdists = Gamma.(shape, scale)
# zdists = Normal.((μw/3).^(1/4))
# zdists = GeneralizedExtremeValue.(μw, sqrt.(vw), .0)
zdists = InverseGaussian.(μw, μw.^3 ./ vw)
w2 = rand.(zdists)
y2 = L2*w2
y = zeros(npts2)
y[dlocs2] = z2
y[slocs2] = y2
plot(domain(sol2), zcolor=y, clims=(0, 2500), 
    markerstrokewidth=0, markershape=:square,
    markersize=2.5, xlabel="Easting (km)", ylabel="Northing (km)", size=(1000, 1000))
scatter!(acoustics.x, acoustics.y, color=:white, markersize=1.2, markerstrokewidth=0)


yy = [L2 * rand.(zdists) for _ in 1:100]
yy = map(1:100) do _
    y = zeros(npts2)
    y[dlocs2] = z2
    y[slocs2] = L2 * rand.(zdists)
    return y
end



pm = histogram(mean.(yy)); vline!(pm, [mean(acoustics.nasc)], title="Mean");
pv = histogram(std.(yy)); vline!(pv, [std(acoustics.nasc)], title="Std. dev.");
p01 = histogram(quantile.(yy, 0.1)); vline!(p01, [quantile(acoustics.nasc, 0.1)]);
p90 = histogram(quantile.(yy, 0.9)); vline!(p90, [quantile(acoustics.nasc, 0.9)]);
plot(pm, pv, p01, p90, legend=false)
plot(pm, pv, legend=false, xlabel="NASC")

pqq = plot()
for y in yy
    qqplot!(pqq, acoustics.nasc, y, color=:black, alpha=0.3)
end
pqq

#############################################################
using StatsBase

predict_ts(L) = -66 * 20log10(L)
predict_weight(L) = 1e-5 * L^2.9

scaling = CSV.read(joinpath(@__DIR__, "scaling_2018.csv"), DataFrame)
trawl_locations = CSV.read(joinpath(@__DIR__, "trawl_locations_2018.csv"), DataFrame)
utmzone = 3
trawl_locations = @chain trawl_locations begin
    @transform(:lla = LLA.(:EQLatitude, :EQLongitude, 0.0))
    @transform(:utm = [UTM(x, utmzone, true, wgs84) for x in :lla])
    @transform(:x = [u.x / 1e3 for u in :utm], :y = [u.y / 1e3 for u in :utm])
end
trawl_locations[!, :edsu_idx] .= 0
for i in 1:nrow(trawl_locations)
    y = trawl_locations.y[i]
    trawl_locations[i, :edsu_idx] = argmin(norm.(y - s for s in acoustics.y))
end
scaling = @chain scaling begin 
    @subset(:class .== "PK1", :analysis_id .== 7, :species_code .== 21740)
    @transform(:weight = predict_weight.(:primary_length))
end

trawl_means = @chain scaling begin
    @by(:event_id, :sigma_bs = mean(:sigma_bs, Weights(:w)),
                   :length = mean(:primary_length, Weights(:w)),
                   :weight = mean(:weight, Weights(:w)))
    @transform(:ts = 10log10.(:sigma_bs))
    leftjoin(trawl_locations, on=:event_id)
end
geotrawl_means = @chain trawl_means begin
    @select(:x, :y, :ts, :length, :weight) 
    @transform(:ts_centered = :ts .- mean(:ts),
               :weight_centered = :weight .- mean(:weight))
    georef((:x, :y))
end
plot(geotrawl_means)

vg_ts = EmpiricalVariogram(geotrawl_means, :ts)
plot(vg_ts)
vg_weight = EmpiricalVariogram(geotrawl_means, :weight)
plot(vg_weight)
vg_length = EmpiricalVariogram(geotrawl_means, :length)
plot(vg_length)


tvg_ts = fit(Variogram, vg_ts)
tvg_ts = SphericalVariogram(nugget=0., sill=0.5, range=50)
prob_ts = SimulationProblem(geotrawl_means, surveydomain, :ts_centered=>Float64, length(yy))
sol_ts = solve(prob_ts, LUGS(:ts_centered => (variogram=tvg_ts,)))
plot(sol_ts[1], markerstrokewidth=0, markershape=:square, markersize=1.5)
ts_mean = mean(trawl_means.ts)

tvg_weight = fit(SphericalVariogram, vg_weight)
tvg_weight = SphericalVariogram(nugget=0., sill=sill(tvg_weight), 
    range=50.)
prob_weight = SimulationProblem(geotrawl_means, surveydomain, :weight_centered=>Float64, length(yy))
sol_weight = solve(prob_weight, LUGS(:weight_centered => (variogram=tvg_weight,)))
plot(sol_weight[1], markerstrokewidth=0, markershape=:square, markersize=1.5)
weight_mean = mean(trawl_means.weight)

sigma_bs_sim = [exp10.((DataFrame(s).ts_centered .+ ts_mean) ./ 10) for s in sol_ts]
weight_sim = [DataFrame(s).weight_centered .+ weight_mean for s in sol_weight]

dA = dx*dy
N_sim = [sum(nasc ./ (4π * σ)) * dA for (nasc, σ) in zip(yy, sigma_bs_sim)]
histogram(N_sim)
std(N_sim) / mean(N_sim)
biomass_sim = [nasc ./ (4π * σ) .* W * dA 
    for (nasc, σ, W) in zip(yy, sigma_bs_sim, weight_sim)]
B_sim = [sum(B) for B in biomass_sim]

histogram(B_sim/1e3, xlabel="Total biomass (T)", legend=false)
std(B_sim) / mean(B_sim)

nasc_std = vec(std(reduce(hcat, yy), dims=2))
weight_std = vec(std(reduce(hcat, weight_sim), dims=2))
sigma_bs_std = vec(std(reduce(hcat, sigma_bs_sim), dims=2))
biomass_std = vec(std(reduce(hcat, biomass_sim), dims=2))

pnasc = scatter(surveydomain, zcolor=nasc_std, title="NASC", clims=(0, 1000),
    markershape=:square, markersize=1.5, markerstrokewidth=0)
psigma = scatter(surveydomain, zcolor=sigma_bs_std, title="σ_bs",
    markershape=:square, markersize=1.5, markerstrokewidth=0)
pweight = scatter(surveydomain, zcolor=weight_std, title="weight",
    markershape=:square, markersize=1.5, markerstrokewidth=0)
pbiomass = scatter(surveydomain, zcolor=biomass_std, clims=(0, 2e6),
    title="Biomass",
    markershape=:square, markersize=1.5, markerstrokewidth=0)

plot(pnasc, psigma, pweight, pbiomass, size=(1000, 1000))

################################
# Age composition
################################
using NearestNeighbors

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

all_ages = @chain Iterators.product(unique(scaling.event_id), 1:AGE_MAX) begin
    DataFrame()
    rename([:event_id, :age])
    sort()
end
age_comp = @chain scaling begin
    @transform(:age = predict_age.(:primary_length))
    @by([:event_id, :age], :p_age=sum(:w))
    DataFrames.groupby(:event_id)
    @transform(:p_age = :p_age / sum(:p_age))
    rightjoin(all_ages, on=[:event_id, :age])
    @transform(:p_age = replace(:p_age, missing => 0.0))
end



tree = KDTree(coordinates.(domain(geotrawl_means)))
idx, dists = knn(tree, coordinates.(surveydomain), 89)


scatter(surveydomain, zcolor=getindex.(idx, argmin.(dists)), markerstrokewidth=0, markersize=2, c=:prism);
scatter!(domain(geotrawl_means))

scatter(surveydomain, zcolor=minimum.(dists), markerstrokewidth=0, markersize=2);
scatter!(domain(geotrawl_means))

a = 5
idx1 = [sample(i, Weights(1 ./ d.^a)) for (i, d) in zip(idx, dists)]
trawl_assignment = trawl_locations.event_id[idx1]
scatter(surveydomain, zcolor=trawl_assignment, markerstrokewidth=0, markersize=2,
    c=:prism);
scatter!(domain(geotrawl_means))

trawl_assignments = map(1:length(yy)) do _ 
    idx1 = [sample(i, Weights(1 ./ d.^a)) for (i, d) in zip(idx, dists)]
    return trawl_locations.event_id[idx1]
end


#=
calculate nasc and mean TS variograms
define simulation grid
define simulation problems for nasc and TS
use nasc covariance matrix to get L, mean, and variance
define white-noise distributions
define KDTree for trawl locations

for i in 1:nreplicates
    nasc = L * rand.(zdists)
    nasc .*= exp10(rand(calibration_error) / 10)
    nasc .*= rand(dead_zone_proportion)

    trawl_ts = georef(bootstrap mean ts from each trawl)
    mean_ts = simulate(trawl_ts, ts_problem)
    Ntotal = nasc ./ (4pi * mean_ts)

    ptrawls = bootstrap one set of trawl proportions (incorporate selectivity bootstrap?)
    trawl_assignments = IDW sample for each sim location
    N_at_age = Ntotal * ptrawls[trawl_assignments]
end
=#