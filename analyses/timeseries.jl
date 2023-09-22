using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using Statistics, StatsBase
using StatsPlots, StatsPlots.PlotMeasures
using ColorSchemes

ebs_surveys = ["200707", "200809", "200909", "201006", "201207", "201407", "201608",
    "201807", "202207"]
ebs_result_files = joinpath.(@__DIR__, "results_" .* ebs_surveys .* ".csv")
# filter(f -> contains(f, "results"), readdir(@__DIR__, join=true))

results = map(ebs_result_files) do f
    survey = basename(f)[9:14]
    println(survey)
    df = @chain CSV.read(f, DataFrame) begin
        @transform(:survey = survey,
                   :year = parse.(Int, first.(survey, 4)))
        @subset(:age .> 0)
        stack([:n_age, :biomass_age]) 
    end
    return df
end
results = vcat(results...)

year_age = @chain results begin
    @subset(:variable .== "n_age")
    @by([:year, :age], :cv = std(:value) / mean(:value))
end
histogram(year_age.cv)
quantile(year_age.cv, [0.1, 0.9])
unstack(year_age, :age, :year, :cv)

eva = DataFrame(
    year = [2007, 2008, 2009, 2010, 2012, 2014, 2016, 2018, 2022],
    cv_1d =   round.([.045, .076, .088, .060, .042, .046, .021, .044, .068] * 100, digits=1)
)

totals = @by(results, [:survey, :year, :variable, :i],
    :value = sum(:value) / 1e9)

p1 = @df @subset(totals, :variable.=="n_age") density(:value, group=:year,
    xlabel="Abundance (billions)", palette=:Paired_9)
p2 = @df @subset(totals, :variable.=="biomass_age") density(:value, group=:year,
    xlabel="Biomass (MT)", palette=:Paired_9)
plot(p1, p2, fill=true, fillalpha=0.5, size=(1000, 500), margin=20px)
savefig(joinpath(@__DIR__, "total_cv_distributions.png"))

qqps = map(unique(totals.year)) do year
    df = @subset(totals, :year .== year)
    # p1 = @df @subset(df, :variable.=="n_age") qqnorm(:value)
    @df @subset(df, :variable.=="biomass_age") qqnorm(:value, title=year)
end
plot(qqps..., size=(800, 800))

annual = @chain results begin
    @by([:survey, :year, :variable, :i], 
        :value = sum(:value) / 1e9)
    @by([:survey, :year, :variable],
        :upper = quantile(:value, 0.95),
        :lower = quantile(:value, 0.05),
        :std = std(:value),
        :value = mean(:value))
    @transform(:cv = round.(:std ./ :value * 100, digits=1))
    leftjoin(eva, on=:year)
    @transform(:cvstring = " " .* string.(:cv) .* " (" .* string.(:cv_1d) .* ")")
end

p_n = @df @subset(annual, :variable .== "n_age") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), 
    series_annotation=text.(:cvstring, :left, :bottom, 10),
    # series_annotation=text.(" " .* string.(:cv), :left, :bottom, 10),
    marker=:o, color=1, label="",
    xticks=2007:2022, xlims=(2006.5, 2024), ylims=(0, 32),
    xlabel="Year", ylabel="Abundance (billions)")
p_b = @df @subset(annual, :variable .== "biomass_age") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), marker=:o, color=2, label="",
    series_annotation=text.(:cvstring, :left, :bottom, 10),
    # series_annotation=text.(" " .* string.(:cv), :left, :bottom, 10),
    xticks=2007:2022, xlims=(2006.5, 2024), ylims=(0, 7.5),
    xlabel="Year", ylabel="Biomass (MT)")
plot(p_n, p_b, layout=(2, 1), size=(800, 600), margin=20px, dpi=300)
savefig(joinpath(@__DIR__, "timeseries.png"))


@df annual scatter(:cv, :cv_1d, group=:variable, label=["Biomass" "Abundance"],
    color=[2 1], xlims=(0, 30), ylims=(0, 30))

@df annual plot(:year, :cv, group=:variable, label=["Biomass" "Abundance"], color=[2 1],
    marker=:o, xticks=2007:2022, ylabel="CV (%)", size=(800, 600), margin=20px)
@df eva plot!(:year, :cv_1d, label="1D geostats", marker=:o)

@by(annual, :variable, :ratio = median(:cv ./ :cv_1d))



plots_n = map(unique(results.year)) do year
    xt = year == 2022 ? collect(1:10) : false
    df = @subset(results, :variable.=="n_age", :year .== year)
    p = @df df violin(:age, :value / 1e9, group=:age, color=1, linewidth=0, label="",
        xticks=xt, xlabel="", ylabel="")
end
plot(plots_n..., layout=(10, 1), size=(600, 900), left_margin=20px)
savefig(joinpath(@__DIR__, "abundance.png"))

plots_b = map(unique(results.year)) do year
    xt = year == 2022 ? collect(1:10) : false
    df = @subset(results, :variable.=="biomass_age", :year .== year)
    p = @df df violin(:age, :value / 1e9, group=:age, color=2, linewidth=0, legend=false,
        xticks=xt)#xlabel="Age class", ylabel="Million tons")
    # scatter!(p, [1], [0], alpha=0, label=year)
end
plot(plots_b..., layout=(10, 1), size=(600, 900), left_margin=20px)
savefig(joinpath(@__DIR__, "biomass.png"))


ebs_error_files = joinpath.(@__DIR__, "stepwise_error_" .* ebs_surveys .* ".csv")
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
    :variable=replace.(:variable, "n" => "Abundance", "biomass" => "Biomass"),
    :error_label = replace.(:error_label, "Trawl jackknife" => "Trawl dropping"))


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
        "Trawl dropping", "Trawl assignment", "TS models", "Age-length", "Length-weight", "All"])
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
        # levels=["Calibration", "Spatial sampling", "Trawl jackknife", "Trawl assignment", 
        # "Selectivity", "Resample catches", "Length-weight", "TS models", "Age-length", "All"])
        levels=["Calibration", "Spatial sampling", "Selectivity", "Resample catches", 
        "Trawl dropping", "Trawl assignment", "TS models", "Age-length", "Length-weight", "All"])
)


pe1 = @df @subset(error_series, :variable.=="Abundance") plot(:year, :cv, group=:error_label,
    ylabel="C.V. (Numbers)", palette=:tableau_10)
pe2 = @df @subset(error_series, :variable.=="Biomass") plot(:year, :cv, group=:error_label,
    ylabel="C.V. (Biomass)", palette=:tableau_10)
plot(pe1, pe2, marker=:o, layout=(2,1), legend=:outerright, ylims=(-0.01, 0.35),
    linewidth=3, xticks=2007:2022, size=(1000, 600), margin=20px)
savefig(joinpath(@__DIR__, "error_timeseries.png"))


error_summary = @chain error_series begin
    @by([:error_label, :variable], 
        :cv_std = std(:cv),
        :cv_max = quantile(:cv, 0.75),
        :cv_min = quantile(:cv, 0.25),
        :cv_se = std(:cv) / sqrt(length(:cv)),
        :cv = mean(:cv),
        :cv_med = median(:cv))
end 
@df error_summary groupedbar(:error_label, :cv, yerror=(:cv_min, :cv_max), group=:variable,
    xrotation=45, xlabel="Individual uncertainty source", ylabel="C.V.",
    size=(600, 500))


@df @orderby(error_series, :error_label) groupedboxplot(:error_label, :cv * 100, group=:variable, 
    outliers=false, permute=(:x, :y), xflip=true, legend=:right,
    ylabel="C.V. (%)", size=(600, 500), dpi=300)
savefig(joinpath(@__DIR__, "error_sources.png"))
