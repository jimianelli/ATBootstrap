"""
Plot the empirical and fitted model variogram for each acoustic scaling stratum in an 
`ATBootstrapProblem`. Optional plotting parameters can be passed in as keyword arguments.
"""
function plot_class_variograms(atbp::ATBootstrapProblem; size=(800, 600), kwargs...)
    pp = map(atbp.class_problems) do cp
        vg_emp = cp.variogram.empirical
        vg_mod = cp.variogram.model
        plot(vg_emp.abscissa, vg_emp.ordinate, title=cp.class, marker=:o,
            label="Empirical", xlabel="Lag (km)", ylabel="γ", legend=:bottomright)
        plot!(h -> vg_mod(h), 0, maximum(vg_emp.abscissa), label="Model")
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
        plot_simulated_nasc(classprob, surveydata, simdomain)
    end
    return plot(plots...; kwargs...)
end

"""
Make violin plots of pollock abundance- and biomass-at-age from the data frame `results`,
the output of running `simulate`. Plotting options can be passed in as keyword arguments.
"""
function plot_boot_results(results; size=(900, 400), margin=15px, palette=:Paired_10,
        kwargs...)
    pk_results = @subset(results, :species_code .== 21740)
    xticks = sort(unique(pk_results.age))
    p_abundance = @df pk_results violin(:age, :n/1e9, group=:age, palette=palette,
        xlabel="Age class", ylabel="Abundance (billions)", legend=false);
    p_biomass = @df pk_results violin(:age, :biomass/1e9, group=:age, palette=palette,
        xlabel="Age class", ylabel="Biomass (Mt)");
    plot(p_abundance, p_biomass; xticks=xticks, size=size, margin=margin, kwargs...)
end


function plot_error_sources(stds_boot; xlims=(-0.005, 0.2), size=(700, 600), 
        kwargs...)
    p1 = @df stds_boot boxplot(:error_label, :n_cv, permute=(:x, :y), xflip=true,
        outliers=false, ylabel="C.V. (Numbers)");
    p2 = @df stds_boot boxplot(:error_label, :biomass_cv, permute=(:x, :y), xflip=true,
        outliers=false, ylabel="C.V. (Biomass)");
    plot(p1, p2; layout=(2,1), legend=false, ylabel="Error source",
        size=size, xlims=xlims, kwargs...)
end