using CSV, DataFrames, DataFramesMeta
using Statistics
using StatsPlots, StatsPlots.PlotMeasures
using ColorSchemes

result_files = filter(f -> endswith(f, ".csv"), readdir(@__DIR__, join=true))

results = map(result_files) do f
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

annual = @chain results begin
    @by([:survey, :year, :variable, :i], 
        :value = sum(:value) / 1e9)
    @by([:survey, :year, :variable],
        :upper = quantile(:value, 0.95),
        :lower = quantile(:value, 0.05),
        :std = std(:value),
        :value = mean(:value))
    @transform(:cv = round.(:std ./ :value * 100, digits=1))
end

p_n = @df @subset(annual, :variable .== "n_age") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), 
    series_annotation=text.(" " .* string.(:cv), :left, :bottom, 10),
    marker=:o, color=1, label="",
    xticks=2012:2022, ylims=(0, 32),
    xlabel="Year", ylabel="Abundance (billions)")
p_b = @df @subset(annual, :variable .== "biomass_age") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), marker=:o, color=2, label="",
    series_annotation=text.(" " .* string.(:cv), :left, :bottom, 10),
    xticks=2012:2022, ylims=(0, 7.5),
    xlabel="Year", ylabel="Biomass (MT)")
plot(p_n, p_b, layout=(2, 1), size=(800, 600), margin=20px)
savefig(joinpath(@__DIR__, "timeseries.png"))


plots_n = map(unique(results.year)) do year
    df = @subset(results, :variable.=="n_age", :year .== year)
    p = @df df violin(:age, :value / 1e9, group=:age, color=1, linewidth=0, label="",
        xticks=1:10, xlabel="", ylabel="")
    # scatter!(p, [1], [0], alpha=0, label=year)
end
plot(plots_n..., layout=(5, 1), size=(600, 600), left_margin=20px)
savefig(joinpath(@__DIR__, "abundance.png"))

plots_b = map(unique(results.year)) do year
    df = @subset(results, :variable.=="biomass_age", :year .== year)
    p = @df df violin(:age, :value / 1e9, group=:age, color=2, linewidth=0, legend=false,
        xticks=1:10)#xlabel="Age class", ylabel="Million tons")
    # scatter!(p, [1], [0], alpha=0, label=year)
end
plot(plots_b..., layout=(5, 1), size=(600, 600), left_margin=20px)
savefig(joinpath(@__DIR__, "biomass.png"))

