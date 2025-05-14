using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using DataFramesMeta: @transform
using Statistics, StatsBase
using GLM
using StatsPlots, StatsPlots.PlotMeasures
using ColorSchemes

ebs_surveys = ["200707", "200809", "200909", "201006", "201207", "201407", "201608",
    "201807", "202207", "202408"]
ebs_result_files = joinpath.(@__DIR__, "results", "results_" .* ebs_surveys .* ".csv")
zdist_files = joinpath.(@__DIR__, "results", "zdists_" .* ebs_surveys .* ".csv")

zdists = map(zdist_files) do f
    CSV.read(f, DataFrame)
end
zdists = vcat(zdists...)
zdists.zdist = replace.(zdists.zdist, "Distributions." => "")
@by(zdists, :zdist, :n = length(:zdist))
@by(zdists, [:class, :zdist], :n = length(:zdist))

# filter(f -> contains(f, "results"), readdir(@__DIR__, join=true))

results = map(ebs_result_files) do f
    survey = basename(f)[9:14]
    println(survey)
    df = @chain CSV.read(f, DataFrame) begin
        DataFramesMeta.@transform(:survey = survey,
                   :year = parse.(Int, first.(survey, 4)))
        @subset(:age .> 0)
        stack([:n, :biomass]) 
    end
    return df
end
results = vcat(results...)

year_age = @chain results begin
    @subset(:variable .== "n")
    @by([:year, :age], :cv = std(:value) / mean(:value))
end
histogram(year_age.cv)
quantile(year_age.cv, [0.1, 0.9])
unstack(year_age, :age, :year, :cv)

# EVA-estimated CVs from cruise reports
eva = DataFrame(
    year = [2007, 2008, 2009, 2010, 2012, 2014, 2016, 2018, 2022, 2024],
    cv_1d =   [3.8, 5.6, 6.9, 5.4, 3.4, 3.4, 1.9, 3.9, 6.8, 5.6],
    n = [10.04, 5.24, 8.63, 12.97, 7.49, 19.51, 12.22, 5.57, 9.67, 11.55],
    biomass = [2.28, 1.404, 1.331, 2.636, 2.279, 4.743, 4.838, 2.497, 3.834, 2.871]
)

totals = @by(results, [:survey, :year, :variable, :i],
    :value = sum(:value) / 1e9)

p1 = @df @subset(totals, :variable.=="n") density(:value, group=:year,
    xlabel="Abundance (billions)", palette=:Paired_9)
p2 = @df @subset(totals, :variable.=="biomass") density(:value, group=:year,
    xlabel="Biomass (MT)", palette=:Paired_9)
plot(p1, p2, fill=true, fillalpha=0.5, size=(1000, 500), margin=20px)
savefig(joinpath(@__DIR__, "plots", "total_cv_distributions.png"))

qqps_abundance = map(unique(totals.year)) do year
    df = @subset(totals, :year .== year)
    @df @subset(df, :variable.=="n") qqnorm(:value, title=year, markerstrokewidth=0,
        markercolor=1, linecolor=:black)
end
plot(qqps_abundance..., size=(800, 800))
savefig(joinpath(@__DIR__, "plots", "qq_normal_plots_abundance.png"))

qqps_biomass = map(unique(totals.year)) do year
    df = @subset(totals, :year .== year)
    @df @subset(df, :variable.=="biomass") qqnorm(:value, title=year, markerstrokewidth=0,
        markercolor=2, linecolor=:black)
end
plot(qqps_biomass..., size=(800, 800))
savefig(joinpath(@__DIR__, "plots", "qq_normal_plots_biomass.png"))

annual = @chain results begin
    @by([:survey, :year, :variable, :i], 
        :value = sum(:value) / 1e9)
    @by([:survey, :year, :variable],
        :upper = quantile(:value, 0.975),
        :lower = quantile(:value, 0.025),
        :std = std(:value),
        :value = mean(:value))
    @transform(:cv = round.(:std ./ :value * 100, digits=1))
    leftjoin(eva, on=:year)
    @transform(:cvstring = " " .* string.(:cv) .* " (" .* string.(:cv_1d) .* ")")
    @orderby(:variable, :year)
end
CSV.write(
    joinpath(@__DIR__, "results", "annual_uncertainty.csv"),
    @select(annual, :year, :variable, :mean=:value, :std, :cv, :upper, :lower)
)

@by(annual, :variable, 
    :mean = mean(:cv),
    :min = minimum(:cv),
    :max = maximum(:cv))
mean(eva.cv_1d)

p_n = @df @subset(annual, :variable .== "n") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), 
    series_annotation=text.(first.(split.(:cvstring, "(")), :left, :bottom, 9),
    marker=:o, color=1, label="", ylims=(0, 30),
    xlabel="", ylabel="Abundance (billions)")
p_b = @df @subset(annual, :variable .== "biomass") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), marker=:o, color=2, label="",
    series_annotation=text.(:cvstring, :left, :bottom, 9), ylims=(0, 6),
    xlabel="Year", ylabel="Biomass (MT)")
# plot!(p_n, eva.year, eva.n, linestyle=:dash, color=1, label="Survey report")
# plot!(p_b, eva.year, eva.biomass, linestyle=:dash, color=2, label="Survey report")
plot(p_n, p_b, layout=(2, 1), size=(800, 600), margin=20px, dpi=300,
    xticks=2008:2:2024, xlims=(2006.5, 2026))
savefig(joinpath(@__DIR__, "plots", "timeseries.png"))


# 1D vs Bootstrap
# Plot time series of CV
p_cv = @df annual plot(:year, :cv, group=:variable, ylims=(0, 35),
    label=["Bootstrap (biomass)" "Bootstrap (abundance)"], color=[2 1],
    marker=:o, xticks=2008:2:2024, xlabel="Year", ylabel="CV (%)",
    title="(a)", titlealign=:left)
@df eva plot!(p_cv, :year, :cv_1d, label="1D geostatistical", marker=:o)

# Fit bootstrap vs. EVA linear models
m0 = lm(@formula(cv ~ 1), @subset(annual, :variable .== "biomass"))
m_biomass = lm(@formula(cv ~ cv_1d), @subset(annual, :variable .== "biomass"))
r2(m_biomass)
ftest(m0.model, m_biomass.model)
b0, b1 = coef(m_biomass)

df_pred = DataFrame(
    cv_1d = 0:8.0,
)
df_pred = [df_pred predict(m_biomass, df_pred, interval=:confidence)]

# plot regressions
p_reg = @df df_pred plot(:cv_1d, :prediction, 
    ribbon=(:prediction .- :lower, :upper .- :prediction), 
    color=[2], fillalpha=0.2, label="");
@df @subset(annual, :variable .== "biomass") scatter!(p_reg, :cv_1d, :cv, label="",
    color=[2], xlabel="1D CV (%)", ylabel="Bootstrap CV (%)",
    title="(b)", titlealign=:left, legend=:bottomright);
annotate!(p_reg, [5], [0],
    text("y = $(round(b0, sigdigits=2)) + $(round(b1, sigdigits=2))x", pointsize=12,
    color=palette(:default)[2]))
plot(p_cv, p_reg, size=(900, 350), margin=15px)
savefig(joinpath(@__DIR__, "plots", "bootstrap_vs_1d.png"))

@by(annual, :variable, :ratio = mean(:cv ./ :cv_1d))

annual_age = @chain results begin
    @subset(:species_code .== 21740)
    @by([:year, :variable, :age], :mean = mean(:value), :std=std(:value))
    @transform(:cv = :std ./ :mean)
    @transform(:upper = :mean .+ 2 * :std, :lower = :mean .- 2 * :std)
end
CSV.write(joinpath(@__DIR__, "results", "age_classes_uncertainty.csv"), annual_age)

@by(annual_age, :variable, 
    :mean = mean(:cv),
    :lower = quantile(:cv, 0.25),
    :upper = quantile(:cv, 0.75))

plots_n = map(unique(results.year)) do year
    df = @subset(results, :variable.=="n", :year .== year)
    df1 = @subset(annual_age, :variable .== "n", :year .== year)
    cv_labels = text.(string.(round.(Int, df1.cv * 100)), :grey, :bottom, 8)
    p = @df df violin(:age, :value / 1e9, group=:age, color=1, linewidth=0, label="",
        xticks=false, foreground_color_legend = nothing)
    @df df1 annotate!(p, :age .+ 0.5, :mean / 1e9, cv_labels, label="")
    scatter!(p, [1], [0], alpha=0, label=year)
    ylims!(p, (0, quantile(df.value, 0.999)/1e9))
    if year == 2012
        ylabel!(p, "Abundance (billions)")
    end
    if year == 2024
        xlabel!(p, "Age class")
        xticks!(p, 1:10, [string.(1:9); "10+"])
    end
    p
end
pn = plot(plots_n..., layout=(10, 1), size=(300, 900), left_margin=20px)

plots_b = map(unique(results.year)) do year
    df = @subset(results, :variable.=="biomass", :year .== year)
    df1 = @subset(annual_age, :variable .== "biomass", :year .== year)
    cv_labels = text.(string.(round.(Int, df1.cv * 100)), :grey, :bottom, 8)
    p = @df df violin(:age, :value / 1e9, group=:age, color=2, linewidth=0, label="",
        xticks=false, foreground_color_legend = nothing)
    @df df1 annotate!(p, :age .+ 0.5, :mean / 1e9, cv_labels, label="")
    scatter!(p, [1], [0], alpha=0, label=year)
    ylims!(p, (0, quantile(df.value, 0.999)/1e9))
    if year == 2012
        ylabel!(p, "Biomas (Mt)")
    end
    if year == 2024
        xlabel!(p, "Age class")
        xticks!(p, 1:10, [string.(1:9); "10+"])
    end
    p
end
pb = plot(plots_b..., layout=(10, 1), size=(300, 900), left_margin=20px)

plot(pn, pb, size=(1000, 1200))
savefig(joinpath(@__DIR__, "plots", "age_classes.png"))


pn = @df @subset(annual_age, :variable .== "n") scatter(:age, :year, 
    markersize=:upper/5e8, legend=false, title="(a)", titlealign=:left)
@df @subset(annual_age, :variable .== "n") scatter!(pn, :age, :year, 
    markersize=:lower/5e8, color=:white)
pb = @df @subset(annual_age, :variable .== "biomass") scatter(:age, 
    :year, markersize=:upper/1e8, color=2, legend=false, title="(b)", titlealign=:left)
@df @subset(annual_age, :variable .== "biomass") scatter!(pb, :age, :year, markersize=:lower/1e8,
    color=:white)
plot(pn, pb, size=(1000, 500), xticks=(1:10, [string.(1:9); "10+"] ),
    yticks=2008:2:2024, xlims=(0, 11), ylims=(2006, 2025), margin=15px,
    xlabel="Age class", ylabel="Year")
savefig(joinpath(@__DIR__, "plots", "age_classes_bubbles.png"))


ebs_error_files = joinpath.(@__DIR__, "results", "stepwise_error_" .* ebs_surveys .* ".csv")
errors = map(ebs_error_files) do f
    survey = basename(f)[16:21]
    println(survey)
    df = @chain CSV.read(f, DataFrame) begin
        DataFramesMeta.@transform(:survey = survey,
                    :year = parse.(Int, first.(survey, 4)))
        # @subset(:age .> 0)
        stack([:n, :biomass]) 
    end
    return df
end
errors = vcat(errors...)

errors = @transform(errors, 
    :variable=replace.(:variable, "n" => "Abundance", "biomass" => "Biomass"))


function resample_df(df, stochastic=true)
    n = nrow(df)
    if stochastic
        ii = sample(1:n, n)
    else
        ii = 1:n
    end
    return @view df[ii, :]
end

stds_boot = map(1:1000) do i
    df = resample_df(errors)
    @by(df, [:year, :error_label, :variable],
        :cv = std(:value) / mean(:value)
    )
end 
stds_boot = vcat(stds_boot...)
summary_boot = @chain stds_boot begin
    @by([:year, :error_label, :variable],
        :cv = mean(:cv),
        :cv_se = std(:cv)
    )
    @transform(:error_label = CategoricalArray(:error_label, 
        # levels=["Calibration", "Spatial sampling", "Trawl jackknife", "Trawl assignment", 
        # "Selectivity", "Resample catches", "Length-weight", "TS models", "Age-length", "All"])
        levels=["Calibration", "Spatial sampling", "Selectivity", "Resample catches", 
        "Nearbottom coefs", "Trawl assignment", "TS models", "Age-length", "Length-weight", "All"])
    )
    @orderby(:year)
end


error_series = @chain errors begin
    @by([:year, :added_error, :error_label, :variable],
        :std = std(:value),
        :mean = mean(:value))
    DataFramesMeta.@transform(:cv = :std ./ :mean)
    # @subset(:error_label .!= "All")
end

error_series = DataFramesMeta.@transform(error_series, 
    :error_label = CategoricalArray(:error_label, 
        levels=["Calibration", "Spatial sampling", "Resample catches", "Selectivity",
        "Nearbottom coefs", "TS models",  "Age-length", "Trawl assignment",
        "Length-weight", "All"])
)

pe1 = @df @subset(error_series, :variable.=="Abundance") plot(:year, :cv*100, group=:error_label,
    ylabel="Abundance CV (%)", palette=:tableau_10)
pe2 = @df @subset(error_series, :variable.=="Biomass") plot(:year, :cv*100, group=:error_label,
    ylabel="Biomass CV (%)", palette=:tableau_10)
plot(pe1, pe2, marker=:o, layout=(2,1), legend=:outerright,
    linewidth=3, xticks=2008:2:2024, size=(800, 500), margin=20px)
savefig(joinpath(@__DIR__, "plots", "error_timeseries.png"))

@chain error_series begin
    @subset(:added_error .== "simulate_nasc")
    @orderby(:variable, :cv)
end

@chain error_series begin
    # @subset(:added_error .== "simulate_nasc")
    @by([:variable, :error_label], 
        :mean = mean(:cv),
        :min = minimum(:cv),
        :max = maximum(:cv)
    )
    @orderby(:variable, :mean)
end

error_summary = @chain error_series begin
    @by([:error_label, :variable], 
        :cv_std = std(:cv),
        :cv_max = quantile(:cv, 0.75),
        :cv_min = quantile(:cv, 0.25),
        :cv_se = std(:cv) / sqrt(length(:cv)),
        :cv_mean = mean(:cv),
        :cv_med = median(:cv))
end 

@df @orderby(error_series, :error_label) groupedboxplot(:error_label, :cv * 100, group=:variable, 
    outliers=false, permute=(:x, :y), xflip=true, legend=:right, yminorgrid=true,
    ylabel="CV (%)", size=(500, 400), dpi=300)
savefig(joinpath(@__DIR__, "plots", "error_sources.png"))
