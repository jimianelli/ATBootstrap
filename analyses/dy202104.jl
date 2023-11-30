using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using GeoStats
using Statistics, StatsBase
using Distributions
using Random
using ConcaveHull
using StatsPlots, StatsPlots.PlotMeasures

using Revise
includet(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
using .ATBootstrap

survey = "202104"
surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0 # km
dA = (resolution * km2nmi)^2
preprocess_survey_data(surveydir, resolution)

(; acoustics, scaling, age_length, length_weight, trawl_locations, domain) = read_survey_files(surveydir)


# Only do the outer shelf transects for now
acoustics = @chain acoustics begin
    @subset(:transect .<= 30)
    @by([:transect, :class, :x, :y], 
        :lon = mean(:lon), :lat = mean(:lat), :nasc = mean(:nasc))
end
@df acoustics scatter(:x, :y, group=:class, aspect_ratio=:equal,
    markersize=:nasc/500, markerstrokewidth=0, alpha=0.5, size=(1000, 1000))
@df trawl_locations scatter!(:x, :y, label="", markerstrokewidth=0,
    text=Plots.text.(string.(:event_id), 10))

# Eliminate trawls that were in Shumagins or Shelikof
outer_shelf_trawls = [1:8; 15:37; 47:68]
trawl_locations = @subset(trawl_locations, in(outer_shelf_trawls).(:event_id))
@df trawl_locations scatter!(:x, :y, label="Outer shelf trawls")

scaling = @subset(scaling, in(outer_shelf_trawls).(:event_id))
scaling_classes = ["SS1"]#sort(unique(scaling.class))
acoustics = @subset(acoustics, in(scaling_classes).(:class))


# re-calculate simulation grid after excluding inside transects
domain = get_survey_grid(acoustics, 10, dx=resolution)
@df domain scatter!(:x, :y, markerstrokewidth=0, markersize=1, aspect_ratio=:equal,
    label="Simulation grid", size=(800, 800))
domain = DataFrames.shuffle(domain) # this seems to fix the issue with directional artifacts
domain =  PointSet(Matrix(domain)')

# Set up problem
surveydata = ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations, 
    domain)

class_problems = map(scaling_classes) do class
    println(class)
    return ATBootstrapProblem(surveydata, class, dA, nlags=15, weightfunc=h -> 1/h)
end

# Inspect the variograms to make sure they look ok
plot_class_variograms(class_problems, legend=:bottomright)

# Check out conditional simulations
plot_simulated_nasc(class_problems, surveydata, size=(1000, 600))

# Do the bootstrap uncertainty analysis
results = simulate_classes(class_problems, surveydata, nreplicates = 500)
plot_boot_results(results)
CSV.write(joinpath(@__DIR__, "results", "results_$(survey).csv"), results)

results_summary = @chain results begin
    @orderby(:age)
    @by(:age, 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
end

# One-at-a-time error analysis
results_step = stepwise_error(class_problems, surveydata; nreplicates = 500)

stepwise_summary = @chain results_step begin
    @orderby(:age)
    @by([:added_error, :age], 
        :biomass_age = mean(:biomass_age),
        :std_age = std(:biomass_age), 
        :cv_age = std(:biomass_age) / mean(:biomass_age) * 100)
    leftjoin(select(results_summary, [:age, :std_age]), on=:age, makeunique=true)
    DataFramesMeta.@transform(:std_decrease = :std_age ./ :std_age_1)
end

@df stepwise_summary plot(:age, :std_age/1e9, group=:added_error, marker=:o, 
    markerstrokewidth=0, size=(1000, 600), margin=20px,# yscale=:log10,
    xlabel="Age class", ylabel="S.D. (Biomass, MT)", title=survey)
@df results_summary plot!(:age, :std_age/1e9, linewidth=2, marker=:o, label="All", 
    color=:black)

    
stepwise_totals = @by(results_step, [:added_error, :i], 
    :n = sum(:n_age), 
    :biomass = sum(:biomass_age))
results_totals = @by(results, :i, 
    :n = sum(:n_age), 
    :biomass = sum(:biomass_age),
    :added_error = "All")

results_totals = @chain [results_totals; stepwise_totals] begin
    leftjoin(error_labels, on=:added_error)
end
CSV.write(joinpath(@__DIR__, "results", "stepwise_error_$(survey).csv"), results_totals)

stds_boot = map(1:1000) do i
    df = resample_df(results_totals)
    @by(df, :error_label, 
        :n_cv = iqr(:n) / mean(:n) ,
        :biomass_cv = iqr(:biomass) / mean(:biomass))
end 
stds_boot = vcat(stds_boot...)

plot_error_sources(stds_boot, plot_title=survey, xlims=(-0.005, 0.4))
