


cals = CSV.read(joinpath(@__DIR__, "../surveydata/calibration_results.csv"), DataFrame)
# using StatsPlots
# @df cals plot(:sa_corr_calc, marker=:o, label="Calculated")
# @df cals plot!(:sa_corr_ek80, marker=:o, label="EK80 wizard")
# @df cals scatter(:sa_corr_calc, :sa_corr_ek80)

# @df @subset(cals, :sphere_id.=="WC38.1_33") histogram(:sa_corr_calc)
# @df @subset(cals, :sphere_id.=="WC38.1_52") vline!(:sa_corr_calc)

# std(cals.sa_corr_calc)
# std(skipmissing(cals.sa_corr_ek80))
# std(cals.sa_corr_calc) / sqrt(nrow(cals))
# std(skipmissing(cals.sa_corr_ek80)) / sqrt(sum(!ismissing, cals.sa_corr_ek80))
# std(cals.sa_corr_calc) / sqrt(3)
# std(skipmissing(cals.sa_corr_ek80)) / sqrt(3)
# 10log10(std(exp10.(cals.sa_corr_calc / 10)) / sqrt(3) + 1)
# 10log10(std(skipmissing(exp10.(cals.sa_corr_ek80 / 10))) / sqrt(3) + 1)

# const CAL_ERROR = 0.1

function simulate_cal_error(cal_error, stochastic=true) 
    if stochastic
        return exp10(cal_error*randn() / 10)
    else
        return exp10(0.0)
    end
end

"""
Plot the empirical and fitted model variogram for each acoustic scaling stratum in an 
`ATBootstrapProblem`. Optional plotting parameters can be passed in as keyword arguments.
"""
function plot_class_variograms(atbp::ATBootstrapProblem; size=(800, 600), kwargs...)
    pp = map(atbp.class_problems) do cp
        vg_emp = cp.variogram.empirical
        vg_mod = cp.variogram.model
        x = [lag.val for lag in vg_emp.abscissas]
        plot(x, vg_emp.ordinates, title=cp.class, marker=:o,
            label="Empirical", xlabel="Lag (km)", ylabel="γ", legend=:bottomright)
        plot!(h -> vg_mod(h), 0, maximum(x), label="Model")
    end
    plot(pp...; size=size, kwargs...)
end

"""
    plot_simulated_nasc(scp, surveydata, simdomain=solution_domain(scp); 
        alpha=0.3, markersize=2.2, max_bubblesize=15, kwargs...)

Plot a NASC field simulated based on the variogram model in `scp` and conditional on the 
observed acoustic data in `surveydata`. Observed NASC values are overplotted as circles,
with radius proprotional to backscatter. The optional arguments `alpha` and `max_bubblesize`
control the transparency and scaling of these bubbles. Optional argument `markersize` 
controls the size of the simulated data points. 
"""
function plot_simulated_nasc(scp::ScalingClassProblem, surveydata::ATSurveyData,
        simdomain=solution_domain(scp); 
        alpha=0.3, markersize=2.2, max_bubblesize=15, kwargs...)
    sim_field = simulate_nasc(scp)
    df = @subset(surveydata.acoustics, :class .== scp.class)
    bubble_factor = max_bubblesize ./ maximum(df.nasc)

    p = scatter(simdomain.x, simdomain.y, zcolor=sim_field, clims=(0, quantile(sim_field, 0.999)), 
        markerstrokewidth=0, markershape=:square, title=string(scp.class),
        aspect_ratio=:equal, markersize=markersize, legend=false,
        xlabel="Easting (km)", ylabel="Northing (km)", kwargs...)
    scatter!(p, df.x, df.y, color=:white, markersize=df.nasc*bubble_factor, alpha=alpha,
        markerstrokewidth=0)
    return p
end

"""
    plot_simulated_nasc(atbp, surveydata[, 
        simdomain=solution_domain(first(atbp.class_problems));
        alpha=0.3, markersize=2.2, max_bubblesize=15, kwargs...])

Plot simulated NASC fields for each scaling stratum specified in `atbp`, conditional on the 
observed acoustic data in `surveydata`. Observed NASC values are overplotted as circles,
with radius proprotional to backscatter. The optional arguments `alpha` and `max_bubblesize`
control the transparency and scaling of these bubbles. Optional argument `markersize` 
controls the size of the simulated data points. 
"""
function plot_simulated_nasc(atbp::ATBootstrapProblem, surveydata::ATSurveyData, 
        simdomain=solution_domain(first(atbp.class_problems));
        alpha=0.3, markersize=2.2, max_bubblesize=15, kwargs...)
    plots = map(atbp.class_problems) do classprob
        plot_simulated_nasc(classprob, surveydata, simdomain; alpha=alpha,
            markersize=markersize, max_bubblesize=max_bubblesize)
    end
    return plot(plots...; kwargs...)
end

"""
    plot_geosim_stats(abtp::ATBootstrapProbelm, surveydata::ATSurveyData[, n=200])

Create diagnostic plots for the conditional geostatistical simulations of each scattering 
class specified by `atbp` and `surveydata`. This will generate `n` realizations of the 
simulated fields and for each class, and will plot:
1. A histogram of the observed NASC overlaid with a spaghetti plot of kernel densities for 
the simulated fields,
2. A histogram of the *means* of the `n` simulated fields along with the observed mean, and
3. A histogram of the standard deviations of the `n` simulated fields along with the
observed standard deviation.
"""
function plot_geosim_stats(atbp, surveydata, n=500)
    hist_plots = []
    mean_plots = []
    sd_plots = []

    for prob in atbp.class_problems
        obs_nasc = @subset(surveydata.acoustics, prob.class .== :class).nasc
        sim_nascs = [simulate_nasc(prob) for _ in 1:n]

        qq = 0.01:0.01:0.99
        q_obs = quantile(obs_nasc, qq)
        q_sims = [quantile(nasc, qq) for nasc in sim_nascs]
        ph = errorline(q_obs, hcat(q_sims...), errortype=:percentile, percentiles=[5, 95],
            marker=:o, markerstrokewidth=0, markersize=2, fillalpha=0.5,
            xlabel="Observed", ylabel="Simulated", label="Average QQ", 
            title=prob.class, foreground_color_legend=nothing,
            background_color_legend=nothing)
        plot!(ph, x -> x, 0, quantile(obs_nasc, 0.99), label="1:1")
        push!(hist_plots, ph)

        hist_means = normalize(fit(Histogram, mean.(sim_nascs), nbins=30))
        dev_means = (mean(mean.(sim_nascs)) .- mean(obs_nasc)) / mean(obs_nasc)
        pm = plot(hist_means, linewidth=0, label="Simulated mean",
            xlabel="Mean NASC (m² nmi²)", ylabel="Probability density")
        vline!(pm, [mean(obs_nasc)], linewidth=3, label="Observed mean")
        plot!(pm, [first(hist_means.edges[1])], [0], alpha=0,
            label="Avg. Δ: $(round(dev_means*100, digits=1))%",
            foreground_color_legend=nothing, background_color_legend=nothing)
        push!(mean_plots, pm)

        hist_stds = normalize(fit(Histogram, std.(sim_nascs), nbins=30))
        dev_stds = (mean(std.(sim_nascs)) .- std(obs_nasc)) / std(obs_nasc)
        ps = plot(hist_stds, linewidth=0, label="Simulated S.D.",
            xlabel="Std. dev. NASC (m² nmi²)", ylabel="Probability density")
        vline!(ps, [std(obs_nasc)], linewidth=3, label="Observed S.D.")
        plot!(ps, [first(hist_means.edges[1])], [0], alpha=0,
            label="Avg. Δ: $(round(dev_stds*100, digits=1))%",
            foreground_color_legend=nothing, background_color_legend=nothing)
        push!(sd_plots, ps)
    end
    nclasses = length(atbp.class_problems)
    return plot(hist_plots..., mean_plots..., sd_plots..., margin=15px,
        layout=(3, length(atbp.class_problems)), size=(400*nclasses, 600))
end

"""
Make violin plots of pollock abundance- and biomass-at-age from the data frame `results`,
the output of running `simulate`. Plotting options can be passed in as keyword arguments.
"""
function plot_boot_results(results; size=(900, 400), margin=15px, palette=:Paired_10,
        kwargs...)
    pk_results = @subset(results, :species_code .== 21740)
    xticks = sort(unique(pk_results.age))
    xticklabels = string.(xticks)
    xticklabels[end] *= "+"
    p_abundance = @df pk_results violin(:age, :n/1e9, group=:age, palette=palette,
        xlabel="Age class", ylabel="Abundance (billions)", legend=false);
    p_biomass = @df pk_results violin(:age, :biomass/1e9, group=:age, palette=palette,
        xlabel="Age class", ylabel="Biomass (Mt)", legend=false);
    plot(p_abundance, p_biomass; xticks=(xticks, xticklabels), size=size, margin=margin, kwargs...)
end

"""

"""
function plot_error_sources(results_totals; xlims=nothing, kwargs...)
    stds_boot = map(1:1000) do i
        df = resample_df(results_totals)
        @by(df, :error_label, 
            :n_cv = std(:n) / mean(:n) ,
            :biomass_cv = std(:biomass) / mean(:biomass))
    end 
    stds_boot = vcat(stds_boot...)
    
    if xlims == nothing
        xmax = max(maximum(stds_boot.n_cv), maximum(stds_boot.biomass_cv)) * 1.05
        xlims = (-0.005, xmax)
    end
    p1 = @df stds_boot boxplot(:error_label, :n_cv, permute=(:x, :y), xflip=true,
        outliers=false, ylabel="C.V. (Numbers)");
    p2 = @df stds_boot boxplot(:error_label, :biomass_cv, permute=(:x, :y), xflip=true,
        outliers=false, ylabel="C.V. (Biomass)");
    plot(p1, p2; layout=(2,1), legend=false, ylabel="Error source",
        xlims=xlims, kwargs...)
end

function summarize_bootstrap(results, variable=:n; species_codes=21740)
    in_spp = in(species_codes)
    @chain results begin
        stack([:n, :biomass])
        @subset(in_spp.(:species_code), :variable .== string(variable))
        @orderby(:age)
        @by(:age, 
            :mean = mean(:value),
            :std = std(:value), 
            :cv = std(:value) / mean(:value) * 100
        )
        rename(:mean => variable)
    end
end

function summarize_stepwise_bootstrap(results, variable=:n; species_codes=21740)
    in_spp = in(species_codes)
    @chain results begin
        stack([:n, :biomass])
        @subset(in_spp.(:species_code), :variable .== string(variable))
        @orderby(:age)
        @by([:added_error, :age], 
            :mean = mean(:value),
            :std = std(:value), 
            :cv = std(:value) / mean(:value) * 100
        )
        rename(:mean => variable)
    end
end

function plot_error_source_by_age(results_step, results, variable=:n; species_codes=21740, kwargs...)
    stepwise_summary = summarize_stepwise_bootstrap(results_step, variable; 
        species_codes=species_codes) 
    results_summary = summarize_bootstrap(results, variable; species_codes=species_codes)
    xticks = sort(unique(stepwise_summary.age))
    xticklabels = string.(xticks)
    xticklabels[end] *= "+"
    @df stepwise_summary plot(:age, :std/1e9, group=:added_error, marker=:o, 
        xticks = (xticks, xticklabels),
        markerstrokewidth=0, xlabel="Age class", ylabel="S.D. (Biomass, MT)", kwargs...)
    @df results_summary plot!(:age, :std/1e9, linewidth=2, marker=:o, label="All", 
        color=:black)    
end

function merge_results(results, results_step, variable=:n)
    stepwise_totals = @by(results_step, [:added_error, :i], 
        :n = sum(:n), 
        :biomass = sum(:biomass))
    results_totals = @by(results, :i, 
        :n = sum(:n), 
        :biomass = sum(:biomass),
        :added_error = "All")

    return @chain [results_totals; stepwise_totals] begin
        leftjoin(error_labels, on=:added_error)
    end
function nearest_length(L, unique_lengths)
    d, i = findmin(x -> abs(L - x), unique_lengths)
    return i
end

function make_age_length_function(age_length, age_max=10, stochastic=true)

    key = @chain age_length begin
        DataFramesMeta.@transform(
            :fork_length = round.(:fork_length),
            :age = ifelse.(:age .> age_max, age_max, :age)
        )
        @by([:fork_length, :age], :n = length(:age))
        unstack(:age, :n, fill=0)
    end

    arr = Array(key[:, 2:end])
    arr = arr ./ sum(arr, dims=2)
    dists = [Distributions.Categorical(arr[i, :]) for i in 1:size(arr, 1)]
    # length_tree = KDTree(key.fork_length')

    func = stochastic ? rand : mode
    predict_age = let fl = key.fork_length
        L -> begin
            # i, _ = nn(length_tree, [L])
            i = nearest_length(L, fl)
            return func(dists[i])
        end
    end
    return predict_age
end
end

function nearest_length(L, unique_lengths)
    d, i = findmin(x -> abs(L - x), unique_lengths)
    return i
end

function make_age_length_function(age_length, age_max=10, stochastic=true)

    key = @chain age_length begin
        DataFramesMeta.@transform(
            :fork_length = round.(:fork_length),
            :age = ifelse.(:age .> age_max, age_max, :age)
        )
        @by([:fork_length, :age], :n = length(:age))
        unstack(:age, :n, fill=0)
    end

    arr = Array(key[:, 2:end])
    arr = arr ./ sum(arr, dims=2)
    dists = [Distributions.Categorical(arr[i, :]) for i in 1:size(arr, 1)]
    # length_tree = KDTree(key.fork_length')

    func = stochastic ? rand : mode
    predict_age = let fl = key.fork_length
        L -> begin
            # i, _ = nn(length_tree, [L])
            i = nearest_length(L, fl)
            return func(dists[i])
        end
    end
    return predict_age
end


# function make_weight_function(stochastic=false, sd=0.05)
#     err = stochastic ? (1 + sd*randn()) : 1.0
#     predict_weight(L) = (1e-5 * L^2.9) * err
#     return predict_weight
# end


using GLM
using StatsBase

function make_weight_function(length_weight, stochastic=true, nmin=5)
    n = nrow(length_weight)
    ii = stochastic ? sample(1:n, n) : 1:n
    length_weight_boot = @view length_weight[ii, :]
    model = lm(@formula(log(organism_weight) ~ log(fork_length)), length_weight_boot)
    c = DataFrame(coeftable(model))
    μloga = c[1,2]
    μb = c[2,2]
    a = exp(μloga)
    b =  μb
    Lmax = 100
    all_lengths = DataFrame(fork_length = 1:Lmax)   
    binned = @chain length_weight_boot begin
        rightjoin(all_lengths, on=:fork_length)
        @by(:fork_length, 
            :mean_weight = mean(:organism_weight),
            :n = length(:organism_weight))
        DataFramesMeta.@transform(
            :weight = ifelse.(:n .>= nmin, :mean_weight, a.*:fork_length.^b)
        )
    end
    weight_dict = Dict(zip(binned.fork_length, binned.weight))

    function weight_function(L)
        L1 = min(max(round(L), 1.0), Lmax)
        return weight_dict[L1]
    end
    
    return weight_function 
end

# wf = make_weight_function(length_weight)
# wf1 = make_weight_function(length_weight, false)
# wf(40)
# wf1(40)

# m = lm(@formula(log(organism_weight) ~ log(fork_length)), length_weight)
# c = DataFrame(coeftable(m))
# a = exp(c[1,2])
# b = c[2,2]

# f(L) = a * L^b

# nmin = 10
# all_lengths = DataFrame(fork_length = 1:maximum(length_weight.fork_length))
# binned = @chain length_weight begin
#     rightjoin(all_lengths, on=:fork_length)
#     @by(:fork_length, 
#         :mean_weight = mean(:organism_weight),
#         :n = length(:organism_weight))
#     DataFramesMeta.@transform(:weight = ifelse.(:n .> nmin, :mean_weight, f.(:fork_length)))
# end

# @df length_weight scatter(:fork_length, :organism_weight, markerstrokewidth=0, alpha=0.1,
#     xlabel="Fork length (cm)", ylabel="Weight (kg)", label="Measurements")
# @df binned scatter!(:fork_length, :mean_weight, label="Binned average")
# @df binned scatter!(:fork_length, :weight, label="Average, filled in")
# plot!(L -> a * L^b, 10, 69, label="Allometric model")
# plot!(wf, 10, 69, label="Combined")

# scatter(length_weight.fork_length, residuals(m), markerstrokewidth=0, alpha=0.1)

# glmm_fit_sel=exp(betas[1]+length_bins*betas[2])/(1+exp(betas[1]+length_bins*betas[2]))

function make_selectivity_function(stochastic=true)
    # LFS curve
    μ_lfs = [-2.4664253, 0.2698528]
    Σ_lfs = [0.41987457 -0.0233237;
            -0.0233237  0.001459591]
    D_lfs = MvNormal(μ_lfs, Σ_lfs)
    β_lfs = stochastic ? rand(D_lfs) : μ_lfs

    # AWT curve
    μ_awt = [-1.0558410, 0.1741619]
    Σ_awt = [0.49771768 -0.01810365
            -0.01810365  0.0008656043]
    D_awt = MvNormal(μ_awt, Σ_awt)
    β_awt = stochastic ? rand(D_awt) : μ_awt

    function selectivity(L, survey)
        if survey < 202001
            return exp(β_awt[1]+L*β_awt[2])/(1+exp(β_awt[1]+L*β_awt[2]))
        else
            return exp(β_lfs[1]+L*β_lfs[2])/(1+exp(β_lfs[1]+L*β_lfs[2]))
        end
    end

    return selectivity
end

function apply_selectivity!(scaling, selectivity_function)
    # w = Array{eltype(scaling.w)}(undef, nrow(scaling))
    # for (i, r) in enumerate(eachrow(scaling))
    #     L = r.primary_length
    #     if r.species_code == 21740 && r.class != "BT"
    #         s = selectivity_function(L, r.survey)
    #         scaling.sample_correction_scalar[i] = 1 / s
    #         w[i] = r.catch_sampling_expansion * r.user_defined_expansion *
    #             r.sample_correction_scalar * r.haul_weight
    #     end
    # end
    # scaling[:, :w] .= w

    @eachrow! scaling begin
        L = :primary_length
        if :species_code == 21740 && :class != "BT"
            s = selectivity_function(L, :survey)
            :sample_correction_scalar = 1 / s
            :w = :catch_sampling_expansion * :user_defined_expansion * :sample_correction_scalar * :haul_weight
        end
    end
end




TSSpec(f, s) = (;f, s)

const TS_SE_DEFAULT = 3.0

function euphausiids_15_65mm_38khz(L)
    A = -9.30429983e2
    B = 3.21027896e0
    C = 1.74003785e0
    D = 1.36133896e-8
    E = -2.26958555e-6
    F = 1.50291244e-4
    G = -4.86306872e-3
    H = 7.38748423e-2
    I = -4.08004891e-1
    J = -7.39078690e1
    Lo = 3.835e-2  
    
    #  soundspeed in m/s
    c = 1470.0
    #  ***frequency in kHz***
    nu = 38.0

    #  wavelength in m
    lam = c / (nu * 10^3)
    #  wavenumber in m^-1
    k = 1 / lam
    #  wavenumber in radians per m
    k = 2 * pi * (nu * 10^3) / c

    if L < 1.5
        TS = -105.0
    elseif L > 6.5
        TS = -73.0
    else
        #  convert input length L in cm to m
        L = L / 100

        #  Demer and Conti, 2006, Vol 63: 928-935, Eq. 10
        TS = (A * (log10(B * k * L) / (B * k * L))^C +
              D * ((k * L)^6) +
              E * ((k * L)^5) +
              F * ((k * L)^4) +
              G * ((k * L)^3) +
              H * ((k * L)^2) +
              I * (k * L) +
              J + 20.0 * log10(L/Lo))
    end
    return TS
end

age0_pollock_ts(L) = 20log10(L) - 64.86
arctic_cod_ts(L) = 8.03log10(L) - 60.78
capelin_ts(L) = 20log10(L) - 70.3
chrysaora_melanaster_ts(L) = 10log10(π * (2L)^2) - 86.8
eulachon_ts(L) = 20log10(L) - 84.5
eulachon_new_ts(L) = 20log10(L) - 84.5
generalized_physoclist_ts(L) = 20log10(L) - 67.4
generic_fish_no_swimbladder_ts(L) = 20log10(L) - 83.2
generic_swimbladder_fish_ts(L) = generalized_physoclist_ts(L)
herring_75m_v2_ts(L) = 20log10(L) - log10(1+75/10) - 65.4
myctophids_sleucopsarus_ts(L) = 32.1log(log10(L)) - 64.1
pacific_hake_ts(L) = 20log10(L) - 68.0
sandlance_ts(L) = 56.5log10(L) - 125.1
squids_ts(L) = 20log10(L) - 75.4
standard_pollock_ts(L) = 20log10(L) - 66

const ts_lookup = Dict(                   # TS function (38 kHz, length in cm)    stdev
    "age0_pollock" =>               TSSpec(age0_pollock_ts,                 TS_SE_DEFAULT),
    "arctic_cod" =>                 TSSpec(arctic_cod_ts,                   TS_SE_DEFAULT),
    "capelin" =>                    TSSpec(capelin_ts,                      TS_SE_DEFAULT),
    "chrysaora_melanaster" =>       TSSpec(chrysaora_melanaster_ts,         TS_SE_DEFAULT),
    "eulachon" =>                   TSSpec(eulachon_ts,                     TS_SE_DEFAULT),
    "eulachon_new" =>               TSSpec(eulachon_new_ts,                 TS_SE_DEFAULT),
    "euphausiids_15_65mm_38khz" =>  TSSpec(euphausiids_15_65mm_38khz,       TS_SE_DEFAULT),
    "generalized_physoclist" =>     TSSpec(generalized_physoclist_ts,       TS_SE_DEFAULT),
    "generic_fish_no_swimbladder" =>TSSpec(generic_fish_no_swimbladder_ts,  TS_SE_DEFAULT),
    "generic_swimbladder_fish" =>   TSSpec(generic_swimbladder_fish_ts,     TS_SE_DEFAULT),
    "herring_75m_v2" =>             TSSpec(herring_75m_v2_ts,               TS_SE_DEFAULT),
    "myctophids_sleucopsarus" =>    TSSpec(myctophids_sleucopsarus_ts,      TS_SE_DEFAULT),
    "pacific_hake" =>               TSSpec(pacific_hake_ts,                 TS_SE_DEFAULT),
    "sandlance" =>                  TSSpec(sandlance_ts,                    TS_SE_DEFAULT),
    "squids" =>                     TSSpec(squids_ts,                       TS_SE_DEFAULT),
    "standard_pollock" =>           TSSpec(standard_pollock_ts,             0.14),
)

function make_ts_function()#stochastic=false)
    error_dict = Dict(relationship => randn() * ts_lookup[relationship].s
       for relationship in keys(ts_lookup))
    
    function predict_ts(relationship, L, stochastic=false)
        f, s = ts_lookup[relationship]
        err = stochastic ? error_dict[relationship] : 0.0
        f(L) + err
    end
    return predict_ts
end

to_linear(x) = exp10(x / 10)

# using DataFrames
# using Statistics
# using StatsPlots
using Distributions

species = CSV.read(joinpath(@__DIR__, "..", "surveydata", "species.csv"), DataFrame)
species.nearbottom_group .= "Misc"
species.nearbottom_group[species.species_code .== 21740] .= "Pollock"
species.nearbottom_group[species.species_code .== 21725] .= "Arctic cod"
species.nearbottom_group[10110 .<= species.species_code .<= 10120] .= "Large flatfish"
species.nearbottom_group[30050 .<= species.species_code .<= 30535] .= "Rockfish"


nearbottom_coefs = DataFrame(
    nearbottom_group = ["Pollock", "Arctic cod", "Large flatfish", "Rockfish", "Misc"],
    A = [2.52, 16.39, 0.85, 93.59, 11.63],
    lower = [2.21, 2.84, 0.16, 8.63, 3.92],
    upper = [2.86, 62.26, 1.67, 343.43, 22.09],
    b = [66, 67.4, 67.4, 67.4, 67.4]
)

nearbottom_coefs = DataFramesMeta.@transform(nearbottom_coefs,
    :sd_A = ((:A .- :lower) .+ (:upper .- :A)) / 2 / 2,
)
const nearbottom_intercept = 3.43
# bar(nearbottom_coefs.nearbottom_group, nearbottom_coefs.a)

nearbottom_coefs = leftjoin(species, nearbottom_coefs, on=:nearbottom_group)

function make_nearbottom_dict(stochastic=true)
    if stochastic
        dists = truncated.(Normal.(nearbottom_coefs.A, nearbottom_coefs.sd_A),
            # 0, Inf)
        nearbottom_coefs.lower, nearbottom_coefs.upper)
        aa = rand.(dists)
    else 
        aa = nearbottom_coefs.A
    end
    aa = aa ./ exp10.(-nearbottom_coefs.b/10) / 1852^2 / 4π
    return Dict(zip(nearbottom_coefs.species_code, aa)) 
end

function apply_nearbottom_coefficient!(scaling, nearbottom_dict)
    @eachrow! scaling begin
        if :species_code in keys(nearbottom_dict)
            :user_defined_expansion = nearbottom_dict[:species_code] 
        else
            # if species code isn't in dict, it's speckled eelpout (i.e.,"Misc")
            :user_defined_expansion = nearbottom_dict[24166]
        end
        :w = :catch_sampling_expansion * :user_defined_expansion * :sample_correction_scalar * :haul_weight
    end
end

include("nearbottom.jl")

abstract type AbstractSurveyDomain end

"""
    TransectRibbons([; width=20, buffer=0.1])

Specify how to define a survey domain based on equal-width strips centered on each 
transect. This method corresponds to MACE's usual calculations based on equally-spaced
transects.

# Arguments

- `transect_width` : The nominal spacing between transects, in nautical miles.
- `transect_buffer` : Amount to expand each transect strip side-to-side to ensure they
overlap.  Default is 0.1 (i.e., 10%).
"""
struct TransectRibbons{T} <:AbstractSurveyDomain
    transect_width::T
    buffer::T
end
TransectRibbons(; width=20, buffer=0.1) = TransectRibbons(promote(width, buffer)...)

"""
    SurveyHull([k=10])

Specify how to define a survey domain based on a concave hull that wraps around all the 
acoustic transects. The smoothness of this hull can be adjusted by setting the number of 
neighbors `k`.

"""
struct SurveyHull{T<:Integer} <: AbstractSurveyDomain
    k::T
end
SurveyHull(k=10) = SurveyHull(k)


function transect_ribbon(transect, transect_width, dx, buffer=0.1, ord=:y)
    tr1 = @chain transect begin
        @select(:x, :y, :log)
        stack(Not(ord))
        DataFramesMeta.@transform($ord = round.($ord ./ dx) .* dx)
        @by([ord, :variable], :value = mean(:value))
        unstack()
        @orderby($ord)
        @select(:x, :y)
        unique()
    end

    a = pi/2
    R1 = [cos(a) -sin(a); sin(a) cos(a)]
    R2 = [cos(-a) -sin(-a); sin(-a) cos(-a)]

    w = transect_width / 2 * 1.852 * (1 + buffer)

    X = Array(tr1)'
    v = diff(X, dims=2)
    v = v ./ norm.(eachcol(v))' .* w
    v = [v v[:, end]]
    left = R1 * v .+ X
    right = R2 * v .+ X

    ribbon_bounds = [tuple(x...) for x in [eachcol(left); eachcol(reverse(right, dims=2))]] 
    ribbon_bounds = [ribbon_bounds; (left[:, 1]...,)]
    return PolyArea(ribbon_bounds)
end

function survey_domain(acoustics, method::TransectRibbons, order, dx, dy=dx)
    tr_set = map(unique(acoustics.transect)) do i
        tr = transect_ribbon(@subset(acoustics, :transect .== i),
            method.transect_width, dx, method.buffer, order)
    end
    return GeometrySet(tr_set)
end

function survey_domain(acoustics, method::SurveyHull, order, dx, dy=dx)
    unique_points = @chain acoustics begin
        @select(:x, :y)
        unique()
    end
    v = [[row.x, row.y] for row in eachrow(unique_points)]
    hull = concave_hull(v, method.k)
    return Ngon([Point(x...) for x in hull.vertices]...)
end

"""
    grid_domain(domain, dx[, dy=dx])

Given a `GeometrySet` or `Domain` object, fill it with a dense rectangular grid with 
resolution `dx` and `dy`. The centers of these grid cells are returned in a `DataFrame`.
"""
function grid_domain(domain, dx, dy=dx)
    box = boundingbox(domain)
    xgrid = range(box.min.coords.x.val, box.max.coords.x.val, step=dx)
    ygrid = range(box.min.coords.y.val, box.max.coords.y.val, step=dx)
    grid = DataFrame([(;x, y) for x in xgrid, y in ygrid
        if in(Point(x, y), domain)])
    return grid
end

"""
    get_survey_grid(acoustics[; method=TransectRibbons(), dx=10.0, dy=dx, order=:y]])

Construct a regular grid inside the survey area, defined as the set of ribbon-like regions
with width `transect_width` along each survey transect.

# Arguments
- `acoustics` : `DataFrame` of georeferenced acoustic data. Should have `:x`, `:y`, and 
`:log` columns.
- `method` : How to define the survey domain. Options are `TransectRibbons()` (default) or 
`SurveyHull()`.
- `dx`, `dy` : Grid resolution, in km. Default is 10 km.
- `order` : `Symbol` indicating which column of `acoustics` to use to define the direction
of each transect for the purposes of defining the ribbon's boundaries. Defaults to `:y`,
since most of MACE's surveys have north-south transects. For east-west transects, use `:x`,
and for curving transects use `:log`.

# Returns
A tuple `(grid, domain)`, where `grid` is a `DataFrame` contianing the x and y coordinates
of each grid cell, and `domain` is a `GeometrySet` object that contains the geographic 
boundaries of the survey domain (either as a single polygon hull or a collection of 
transect ribbons).

"""
function get_survey_grid(acoustics; method=TransectRibbons(), dx=10.0, dy=dx, order=:y)
    nrow(acoustics) > 3 || error("Not enough locations to define survey grid.")
    domain = survey_domain(acoustics, method, order, dx, dy)
    grid = grid_domain(domain, dx, dy)
    return grid, domain
end

function make_haul_id(prefix, ship, event_id)
    return prefix .* "-".* string.(ship) .* "-" .* string.(event_id)
end

"""
Merge MACE's "macebase2.scaling_key_source_data" table with GAP's "racebase.specimen"
table to make a combined virtual SKSD table. This assigns all specimens from the GAP 
trawls to a new scaling stratum called "BT", which is applied to the bottom 3 m of the 
water column.
"""
function merge_scaling(scaling_mace, scaling_gap)
    ts_key = @by(scaling_mace, :species_code, :ts_relationship=first(:ts_relationship))
    
    scaling_mace1 = @chain scaling_mace begin
        DataFramesMeta.@transform(:haul_id = make_haul_id("MACE", 1, :event_id))
        @select(:survey, :ship, :haul_id, :class, :species_code, 
            :primary_length, :ts_length, :ts_relationship, :catch_sampling_expansion,
            :user_defined_expansion, :sample_correction_scalar, :haul_weight, :w)
    end

    scaling_gap1 = @chain scaling_gap begin
        DataFramesMeta.@transform(
            :survey = :cruise,
            :ship = :vessel,
            :haul_id = make_haul_id("GAP", :vessel, :haul),
            :class = "BT",
            :primary_length = :length ./ 10,
            :ts_length = :length ./ 10, # not exactly right
            :catch_sampling_expansion = 1.0,
            :sample_correction_scalar = 1.0,
            :user_defined_expansion = 1.0,
            :haul_weight = 1.0,
            :w = 1.0
        )
        leftjoin(ts_key, on=:species_code)
        DataFramesMeta.@transform(
            :ts_relationship = replace(:ts_relationship, missing => "generic_swimbladder_fish")
        )
        @select(:survey, :ship, :haul_id, :class, :species_code,
            :primary_length, :ts_length, :ts_relationship, :catch_sampling_expansion,
            :user_defined_expansion, :sample_correction_scalar, :haul_weight, :w)
        dropmissing()
    end
    return [scaling_mace1; scaling_gap1]
end


function merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
    survey = only(unique(trawl_locations_mace.survey))
    trawl_locations_gap.survey .= survey
    tl_mace = @chain trawl_locations_mace begin
        DataFramesMeta.@transform(:haul_id = make_haul_id("MACE", 1, :event_id))
        @select(:survey, :haul_id, :latitude, :longitude)
    end
    tl_gap = @chain trawl_locations_gap begin
        DataFramesMeta.@transform(:haul_id = make_haul_id("GAP", :vessel, :event_id))
        @select(:survey, :haul_id, :latitude, :longitude)
    end

    return [tl_mace; tl_gap]
end

"""
preprocess_survey_data(surveydir[; ebs=true, log_ranges=nothing,
    grid_method=TransectRibbons(), dx=10.0, dy=dx, transect_order=:y,
    missingstring=[".", "NA"]])

Preprocess the survey data files in directory `surveydir` and save the outputs in 
the same directory as .csv files in standard format. 

# Arguments

## Required arguments

- `surveydir` : Path to the directory where the survey data files are located. This
function expects to find the following files, all of which come from running
"download_survey.R":

- scaling_mace.csv : Specimen data from scaling_key_source_data
- acoustics.csv : Acoustic NASC in 0.5 nmi intervals
- trawl_locations_mace.csv : Lat/Lon locations of all MACE trawl events
- measurements.csv : Length-weight measurements

## Optional arguments

- `ebs` : Boolean flag indicating whether the survey took place in the Eastern Bering Sea.
If `ebs=true`, this function will expect the following two files to be present,
containing data from the Groundfish Assessment Program survey, which are used to scale
acoustic data in the bottom 3 m of the water column:
- trawl_locations_gap.csv : Lat/Lon locations of all bottom trawls
- scaling_gap.csv : Specimen data from racebase.specimen
- `dx` : Set the resolution of the sampling grid; defaults to 10.0 (km).
- `dy` : Set if the N-S resolution is different from the E-W resolution, otherwise will be
equal to `dx`.
- `grid_method` : How to define the survey domain for the purposes of constructing the 
simulation grid. The options are `TransectRibbons()` (the default) or `SurveyHull()`. See
their documentation for details and options.
- `missingstring` : What string(s) in the CSV files should be interpreted as missing data?
Defaults to `[".", "NA"]. (If the files have been downloaded correctly, you should not
need to use this.)
- `transect_order` : `Symbol` indicating which column of `acoustics` to use to define the
direction of each transect for the purposes of defining the domain. Defaults to `:y`,
since most of MACE's surveys have primarily north-south transects. For east-west
transects, use `:x`, and for curving transects use `:log`.

The preprocessing includes the following tasks:

- Excluding NASC from non-scaling classes
- If GAP data are present, merging them with the MACE data into unified `scaling` and
`trawl_locations` data tables.
- Geographic projection of spatial data
- Downscaling acoustic resolution (specified by `dx` and `dy` arguments)
- Calculating survey domain as concave hull of transect ends
- Setting up simulation grid 
- General tidying.

The output from these pre-processing operations is written to new files in the same data
directory:
- scaling.csv : Formatted scaling key data (containing bottom trawls if they were included)
- acoustics_projected.csv : Spatially-projected NASC by interval and scaling class.
"""
function preprocess_survey_data(surveydir; ebs=true, log_ranges=nothing, dx=10.0, dy=dx, 
        missingstring=[".", "NA"],
        grid_method=TransectRibbons(), 
        transect_order=:y)
    scaling_mace = CSV.read(joinpath(surveydir, "scaling_mace.csv"), DataFrame,
        missingstring=missingstring)
    trawl_locations_mace = CSV.read(joinpath(surveydir, "trawl_locations_mace.csv"), DataFrame,
        missingstring=missingstring)

    if ebs
        scaling_gap = CSV.read(joinpath(surveydir, "scaling_gap.csv"), DataFrame,
            missingstring=missingstring)
        scaling = merge_scaling(scaling_mace, scaling_gap)
        trawl_locations_gap = CSV.read(joinpath(surveydir, "trawl_locations_gap.csv"), DataFrame,
            missingstring=missingstring)
        trawl_locations = merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
    else
        scaling = scaling_mace
        trawl_locations = trawl_locations_mace
    end
    
    acoustics = CSV.read(joinpath(surveydir, "acoustics.csv"), DataFrame,
        missingstring=missingstring)

    # Merge "filtered" strata with main ones
    scaling.class .= replace.(scaling.class, "_FILTERED" => "")
    acoustics.class .= replace.(acoustics.class, "_FILTERED" => "")
    scaling_classes = unique(scaling.class)

    if isnothing(log_ranges)
        log_ranges = [tuple(extrema(acoustics.start_vessel_log)...)]
    end

    acoustics = @chain acoustics begin
        @select(:transect,
                :interval, 
                :class, 
                :lat = :start_latitude, 
                :lon = :start_longitude,
                :log = :start_vessel_log,
                :nasc)
        @subset(
            in(scaling_classes).(:class),
            in_intervals(:log, log_ranges),
            abs.(:lon) .< 360,
            abs.(:lat) .< 360
        )
        @by([:transect, :interval, :class, :lat, :lon, :log], :nasc = sum(:nasc))
        unstack([:transect, :interval, :lat, :lon, :log], :class, :nasc, fill=0)
        stack(Not([:transect, :interval, :lat, :lon, :log]), variable_name=:class, value_name=:nasc)
    end
    acoustics.nasc[ismissing.(acoustics.nasc)] .= 0

    utmzone = 3
    lla = LLA.(acoustics.lat, acoustics.lon, 0.0)
    utm = [UTM(x, utmzone, true, wgs84) for x in lla]
    acoustics.x = [u.x / 1e3 for u in utm]
    acoustics.y = [u.y / 1e3 for u in utm]

    surveygrid, surveyhull = get_survey_grid(acoustics, method=grid_method,
        dx=dx, dy=dy, order=transect_order)

    xmin = minimum(acoustics.x)
    ymin = minimum(acoustics.y)

    acoustics = @chain acoustics begin
        # DataFramesMeta.@transform(:x = :x .- xmin, :y = :y .- ymin)
        DataFramesMeta.@transform(
            :x = round.(:x ./ dx) .* dx, 
            :y = round.(:y ./ dy) .* dy
        ) 
        # DataFramesMeta.@transform(:x = :x .+ xmin, :y = :y .+ ymin)
        @by([:transect, :class, :x, :y], 
            :lon = mean(:lon), 
            :lat = mean(:lat), 
            :log = mean(:log),
            :nasc = mean(:nasc)
        )
    end

    trawl_locations = @chain trawl_locations begin
        DataFramesMeta.@transform(:lla = LLA.(:latitude, :longitude, 0.0))
        DataFramesMeta.@transform(:utm = [UTM(x, utmzone, true, wgs84) for x in :lla])
        DataFramesMeta.@transform(
            :x = [u.x / 1e3 for u in :utm], 
            :y = [u.y / 1e3 for u in :utm]
        )
        @subset([in(Point(x, y), surveyhull) for (x, y) in zip(:x, :y)])
    end

    length_weight = CSV.read(joinpath(surveydir, "measurements.csv"), DataFrame,
        missingstring=missingstring)
    rename!(lowercase, length_weight)
    length_weight = @chain length_weight begin
        unstack([:specimen_id, :species_code, :event_id], :measurement_type, :measurement_value)
        dropmissing()
    end


    CSV.write(joinpath(surveydir, "scaling.csv"), scaling)
    CSV.write(joinpath(surveydir, "acoustics_projected.csv"), acoustics)
    CSV.write(joinpath(surveydir, "length_weight.csv"), length_weight)
    CSV.write(joinpath(surveydir, "trawl_locations_projected.csv"), trawl_locations)
    CSV.write(joinpath(surveydir, "surveygrid.csv"), surveygrid)
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
        DataFramesMeta.groupby(df, [:haul_id, :class]))
end

function get_trawl_category_means(scaling, aged_species, predict_weight)
    use_ages = in(aged_species)
    trawl_means_cat = @chain scaling begin
        DataFramesMeta.@transform(
            :category = _category.(use_ages, :species_code, :age),
            :weight = predict_weight.(:primary_length),
            :p_nasc = :sigma_bs .* :w
        )
        #=
        If a catch filter is applied, all specimens may be listed twice in the same haul, 
        with "user_defined_expansion" set to 0.0 for one listing and 1.0 for the other.
        This is done to split a haul into two parts, each applied to a different scaling
        stratum. Eliminate weighting factors == 0 here, since they don't affect any of the
        numbers and introduce NaNs if all specimens in a haul/category have :w .== 0.
        =# 
        @subset(:w .> 0)
        @by([:haul_id, :category], 
            :sigma_bs = mean(:sigma_bs),#, Weights(:w)),
            :p_nasc = sum(:p_nasc),
            :weight = mean(:weight)#, Weights(:w)),
        )
        # make proportions sum to 1.0 for each haul
        DataFrames.groupby(:haul_id)
        DataFramesMeta.@transform(:p_nasc = :p_nasc ./ sum(:p_nasc))
    end
    return trawl_means_cat
end

function _category(use_ages, species_code, age)
    if use_ages(species_code)
        return string.(species_code) .* "@" .* string.(age)
    else
        return string.(species_code) .* "@-1"
    end
end

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
    # add small epsilon to avoid occasional "duplicate coordinate" warnings
    geonasc.x .+= 1e-3 .* randn.()
    geonasc.y .+= 1e-3 .* randn.()
    geonasc = georef(geonasc, (:x, :y))
    evg = EmpiricalVariogram(geonasc, :nasc, nlags=nlags, maxlag=maxlag)
    tvg = GeoStatsFunctions.fit(ExponentialVariogram, evg, weightfunc)
    setup = GeoStatsProcesses.randsetup(sim_domain, geonasc, 1)
    variogram = (empirical = evg, model = tvg)
    return (variogram, setup)
end

"""
    get_lungs_params(geoproblem, theoretical_variogram[, variable=:nasc])

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
function get_lungs_params(setup, theoretical_variogram, variable=:nasc)
    # solver = LUGS(variable => (variogram = variogram,))
    # preproc = preprocess(geoproblem, solver);
    # pars = preproc[Set([variable])][Set([variable])]
    prep = GeoStatsProcesses.randprep(
        Random.TaskLocalRNG(),
        GeoStatsProcesses.GaussianProcess(theoretical_variogram),
        GeoStatsProcesses.LUMethod(),
        setup
    )
    pars = prep[variable]
    # (z₁, d₂, L₂₂, μ, dlocs, slocs)
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
    (; data, μx, L, μ, dlocs, slocs) = lungs_params
    μx = copy(μx)
    μx[μx .< ϵ] .= ϵ 
    μx = μx .* mean(data) ./ mean(μx)
    μz = L \ μx
    μz[μz .<= ϵ] .= ϵ
    vz = ones(length(μz))
    return [Dist(p...) for p in dist_params.(Dist, μz, vz)]
end


function nonneg_lumult(params, z)
    (; data, μx, L, μ, dlocs, slocs) = params
    npts = length(dlocs) + length(slocs)
    x = zeros(npts)
    x[slocs] = L * z
    x[dlocs] = data
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

function zdists(atbp::ATBootstrapProblem)
    nts = [(class = cp.class, zdist = cp.zfamily) for cp in atbp.class_problems]
    return DataFrame(nts)
end

"""
    trawl_assignments(pixel_coords, trawl_coords[, stochastic=true[; nneighbors=4, a=1.9]])

Assign each acoustic cell to a trawl, either deterministically (i.e. a standard MACE 
nearest-trawl assignment) or randomly. By default, the random trawl is drawn from one of
the 4 nearest trawls to the acoustic cell, with probability proportional to 1/distance^2.
The number of neighboring trawls and the inverse-distance exponent can be set using the
optional `nneighbors` and `a` keyword arguments.

Returns a vector the same length as `pixel_coords`, each element of which is an integer 
index pointing to one of the trawls in `trawl_coords`.

# Arguments
- `pixel_coords`, `trawl_coords` : Vectors of coordinates for the acoustic cells and 
    trawl locations. Each element of these is a 2D vector with the x, y coordinates of the
    pixel or trawl.
- `stochastic=true` : Whether trawl assignment should be random (the default) or
    deterministic.
- `nneighbors = 4` : Number of neighboring trawls to draw from
- `a = 2` : Exponent for inverse distance weights for trawl assignment.
"""
function trawl_assignments(pixel_coords, trawl_coords, stochastic=true; 
        nneighbors=4, a=2)
    nneighbors = min(nneighbors, length(trawl_coords))
    # Calculate a k-dimensional tree for efficiently finding nearest neighbors
    kdtree = KDTree(trawl_coords)
    #=
    idx and dists are vectors the same length as the number of pixels/acoustic cells.
    Each element of idx and dists is a vector with `nneighbors` elements.
    The ith element of idx is a vector of indices to each of the trawl locations.
    The ith element of dists is a vector of distances to the trawls indexed by idx.
    This means that the jth element of dists[i] is the distance from pixel i to trawl j.
    =#
    idx, dists = knn(kdtree, pixel_coords, nneighbors)#length(trawl_coords))
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
#=
Type definitions for data structures
=#

"""
    ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations,
        grid, dA)

Construct an `ATSurveyData` object, encapsulating all the relevant data 
"""
struct ATSurveyData
    acoustics::DataFrame
    scaling::DataFrame
    age_length::DataFrame
    length_weight::DataFrame
    trawl_locations::DataFrame
    grid::Meshes.PointSet
    dA::Number
end

"""
    BootSpecs([; selectivity, predict_ts, resample_scaling, nearbottom_coefs, age_length,
        weights_at_age, trawl_assignments, simulate_nasc, calibration])
    
Construct a `BootSpecs` object, which specifies which error sources to include when running
`simulate` on an `ATBootstrapProblem`. All arguments are Booleans and default to true; they
can be 

# Examples

julia> BootSpecs(calibration=false) # turn off calibration error but include everything else

julia> BootSpecs(false) # turn off all errors--i.e., do a normal deterministic analysis

"""
@kwdef struct BootSpecs
    selectivity::Bool=true
    predict_ts::Bool=true
    resample_scaling::Bool=true
    nearbottom_coefs::Bool=true
    age_length::Bool=true
    weights_at_age::Bool=true
    trawl_assignments::Bool=true
    simulate_nasc::Bool=true
    calibration::Bool=true
end
BootSpecs(b::Bool) = BootSpecs(fill(b, length(fieldnames(BootSpecs)))...)

struct ScalingClassProblem
    class
    variogram
    geosetup
    params
    zfamily
    zdists
    cal_error
    age_max
    aged_species
end

"""
    ScalingClassProblem(surveydata, class[, cal_error=0.1, age_max=10, maxlag=200,
        nlags=20, weightfunc=h -> 1/h,
        zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal],
        aged_species=[21740]])

Set up a `ScalingClassProblem`, specifying how to perform geostatistical simulations of
backscatter in scaling stratum `class`, conditional on the observed NASC values in 
`surveydata`. 
"""
function ScalingClassProblem(surveydata, class; 
        cal_error=0.1, age_max=10, maxlag=200.0, nlags=20, weightfunc=h -> 1/h,
        zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal],
        aged_species=[21740])
    acoustics_sub = @subset(surveydata.acoustics, :class .== class)
    variogram, geosetup = define_conditional_sim(acoustics_sub, surveydata.grid,
        maxlag=maxlag, nlags=nlags, weightfunc=weightfunc)
    params = get_lungs_params(geosetup, variogram.model)
    optimal_dist = choose_z_distribution(zdist_candidates, acoustics_sub.nasc, params)
    zdists = parameterize_zdists(optimal_dist, params)
    return ScalingClassProblem(class, variogram, geosetup, params, optimal_dist, zdists,
        cal_error, age_max, aged_species)
end

"""
Extract the simulation domain from a `ScalingClassProblem`. Returns a `DataFrame` with two
columns containing the `x` and `y` coordinates of each point at which the spatial field 
is to be simulated.
"""
function solution_domain(scp::ScalingClassProblem, variable=:nasc)
    # sol = solve(scp.geosetup, LUGS(variable => (variogram = scp.variogram.model,)))
    # dom = domain(sol)
    dom = scp.geosetup.domain
    x = [p.coords.x for p in dom]
    y = [p.coords.y for p in dom]
    return DataFrame(x=x, y=y)
end

struct ATBootstrapProblem
    class_problems
    scaling_classes
    age_max
    aged_species
end

 

"""
    ATBootstrapProblem(surveydata[; scaling_classes, cal_error=0.1, age_max=10,
        zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal],
        maxlag=200, nlags=10, weightfunc=h -> 1/h])

Set up an `ATBootstrapProblem`, describing how to do bootstrap analyses of the 
acoustic-trawl survey recorded in `surveydata` for the scaling strata specified in 
`scaling_classes`.
"""
function ATBootstrapProblem(surveydata::ATSurveyData; age_max=10,
        aged_species=[21740], scaling_classes=unique(surveydata.acoustics.class),
        cal_error=0.1, zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal],
        maxlag=200.0, nlags=10, weightfunc=h -> 1/h)
    
    class_problems = map(scaling_classes) do class
        println("Preparing $(class)...")
        return ScalingClassProblem(surveydata, class, maxlag=maxlag, nlags=nlags,
            age_max=age_max, cal_error=cal_error, weightfunc=weightfunc, 
            zdist_candidates=zdist_candidates, aged_species=aged_species)
    end
    return ATBootstrapProblem(class_problems, scaling_classes, age_max, aged_species)
end
