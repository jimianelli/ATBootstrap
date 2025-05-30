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
end