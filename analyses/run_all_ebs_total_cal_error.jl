using CSV
using DataFrames, DataFramesMeta
using StatsPlots, StatsPlots.PlotMeasures
using Statistics
using Random
using StableRNGs
include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
cal_error = 0.647

ebs_surveys = [
    "200707",
    "200809",
    "200909",
    "201006",
    "201207",
    "201407",
    "201608",
    "201807", 
    "202207",
    "202408" 
]

for survey in ebs_surveys
    println("Running analysis for survey $(survey)")
    surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
    StableRNGs.seed!(parse(Int, survey))
    (; acoustics, scaling, age_length, length_weight, trawl_locations, surveygrid) = ATB.read_survey_files(surveydir)

    if survey == "201807"
        scaling_classes = ["PK1", "BT"] 
    elseif survey == "202207"
        scaling_classes = ["SS1", "SS2", "BT"] 
    else 
        scaling_classes = unique(scaling.class)
    end

    acoustics = @subset(acoustics,
        in(scaling_classes).(:class))
    if survey == "202207"
        acoustics = @subset(acoustics, iseven.(:transect))
    end

    surveydata = ATB.ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
        surveygrid, dA)

    atbp = ATB.ATBootstrapProblem(surveydata, cal_error=cal_error)
    results = ATB.simulate(atbp, surveydata, nreplicates = 500)
    CSV.write(joinpath(@__DIR__, "results", "results_total_cal_uncertianty_$(survey).csv"), results)
end

#=
Plotting results
=#

ebs_result_files = joinpath.(@__DIR__, "results", "results_total_cal_uncertianty_" .* ebs_surveys .* ".csv")
results = map(ebs_result_files) do f
    survey = last(split(f, "_"))[1:6]
    println(survey)
    df = @chain CSV.read(f, DataFrame) begin
        @transform(:survey = survey,
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
    biomass = [2.28, 1.404, 1.331, 2.636, 2.279, 4.743, 4.838, 2.497, 3.834, 2.871]
)

totals = @by(results, [:survey, :year, :variable, :i],
    :value = sum(:value) / 1e9)

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

@by(annual, :variable, :cv = mean(:cv))

official_biomass = CSV.read(joinpath(@__DIR__, "../surveydata/official_biomass.csv"), DataFrame)
official_biomass = @subset(official_biomass, :year .!= 2020)

p_n = @df @subset(annual, :variable .== "n") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), 
    series_annotation=text.(first.(split.(:cvstring, "(")), :left, :bottom, 9),
    marker=:o, color=1, label="",
    xticks=2007:2024, xlims=(2006.5, 2026), ylims=(0, 40),
    xlabel="Year", ylabel="Abundance (billions)")
p_b = @df @subset(annual, :variable .== "biomass") plot(:year, :value, 
    ribbon = (:value .- :lower, :upper .- :value), marker=:o, color=2, label="",
    series_annotation=text.(:cvstring, :left, :bottom, 9),
    xticks=2007:2024, xlims=(2006.5, 2026), ylims=(0, 7.5),
    xlabel="Year", ylabel="Biomass (MT)")
# plot!(p_b, official_biomass.year, official_biomass.biomass,
    # linestyle=:dash, label="Survey report")
plot(p_n, p_b, layout=(2, 1), size=(800, 600), margin=20px, dpi=300)
savefig(joinpath(@__DIR__, "plots", "timeseries_total_cal_error.png"))


annual_age = @chain results begin
    @subset(:species_code .== 21740)
    @by([:year, :variable, :age], :mean = mean(:value), :std=std(:value))
    @transform(:cv = :std ./ :mean)
    @transform(:upper = :mean .+ 2 * :std, :lower = :mean .- 2 * :std)
end

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
savefig(joinpath(@__DIR__, "plots", "age_classes_total_cal_error.png"))