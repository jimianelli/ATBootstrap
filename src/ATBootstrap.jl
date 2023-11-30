module ATBootstrap

using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using GeoStats, Geodesy, ConcaveHull
using Statistics, StatsBase
using Distributions
using LinearAlgebra
using Distances
using NearestNeighbors
using ProgressMeter
using RCall
using StatsPlots, StatsPlots.PlotMeasures
using ColorSchemes

export preprocess_survey_data,
    get_survey_grid,
    read_survey_files,
    ATSurveyData,
    ATBootstrapProblem,
    resample_df,
    nonneg_lusim,
    nonneg_lusim!,
    nonneg_lumult,
    nonneg_lumult!,
    solution_domain,
    simulate_classes,
    stepwise_error,
    BootSpecs,
    error_labels,
    plot_class_variograms,
    plot_simulated_nasc,
    plot_boot_results,
    plot_error_sources

# struct SurveyData
#     acoustics
#     scaling 
#     trawl_locations
# end

include("preprocess_survey_data.jl")
include("spatial.jl")
include("mace_ts.jl")
include("mace_age_length.jl")
include("mace_selectivity.jl")
include("mace_length_weight.jl")
include("calibration.jl")
include("scaling.jl")

struct ATSurveyData
    acoustics
    scaling
    age_length
    length_weight
    trawl_locations
    domain
end

@kwdef struct BootSpecs
    selectivity::Bool=true
    predict_ts::Bool=true
    resample_scaling::Bool=true
    drop_trawl::Bool=true
    age_length::Bool=true
    weights_at_age::Bool=true
    trawl_assignments::Bool=true
    nonneg_lusim::Bool=true
    calibration::Bool=true
end
BootSpecs(b::Bool) = BootSpecs(fill(b, length(fieldnames(BootSpecs)))...)

struct ATBootstrapProblem
    class
    variogram
    problem
    params
    optimal_dist
    zdists
    cal_error
    age_max
    dA
end

function ATBootstrapProblem(surveydata, class, dA; cal_error=CAL_ERROR, age_max=10,
        zdist_candidates=zdist_candidates, maxlag=200.0, nlags=20, weightfunc=h -> 1/h)
    acoustics_sub = @subset(surveydata.acoustics, :class .== class)
    variogram, problem = define_conditional_sim(acoustics_sub, surveydata.domain,
        maxlag=maxlag, nlags=nlags, weightfunc=weightfunc)
    params = get_lungs_params(problem, variogram.model)
    optimal_dist = choose_distribution(zdist_candidates, acoustics_sub.nasc, params)
    zdists = get_zdists(optimal_dist, params)
    return ATBootstrapProblem(class, variogram, problem, params, optimal_dist, zdists,
        cal_error, age_max, dA)
end

function solution_domain(atbp::ATBootstrapProblem, variable=:nasc)
    sol = solve(atbp.problem, LUGS(variable => (variogram = atbp.variogram.model,)))
    dom = domain(sol)
    x = [p.coords[1] for p in dom]
    y = [p.coords[2] for p in dom]
    return DataFrame(x=x, y=y)
end

"""
    simulate(atbp, surveydata; nreplicates=500, bs=BootSpecs())

Do a bootstrap analysis for one acoustic class/scaling stratum using the problem defined
in `atbp` and the data in `surveydata`. The number of bootstrap replicates can optionally
be set in `nreplicates`, and the `BootSpecs` object `bs` can be used to specify if any of
the error sources should be omitted (by default all are included).
"""
function simulate(atbp::ATBootstrapProblem, surveydata; nreplicates=500, bs=BootSpecs())

    scaling_sub = @subset(surveydata.scaling, :class .== atbp.class)

    z0 = mean.(atbp.zdists)

    println("Bootstrapping $(atbp.class)...")
    results = @showprogress map(1:nreplicates) do i
        selectivity_function = make_selectivity_function(bs.selectivity)
        scaling_boot = resample_scaling(scaling_sub, bs.resample_scaling)
        apply_selectivity!(scaling_boot, selectivity_function)

        predict_ts = make_ts_function(bs.predict_ts)
        predict_age = make_age_length_function(surveydata.age_length, atbp.age_max,
            bs.age_length)
        scaling_boot = DataFramesMeta.@transform(scaling_boot,
            :sigma_bs = exp10.(predict_ts.(:ts_relationship, :ts_length)/10),
            :age = predict_age.(:primary_length))
        
        trawl_means = get_trawl_means(scaling_boot, surveydata.trawl_locations)
        if bs.drop_trawl
            popat!(trawl_means, rand(1:nrow(trawl_means)))
            # trawl_means = resample_df(trawl_means)
        end
        geotrawl_means = @chain trawl_means begin
            @select(:x, :y, :ts, :length) 
            georef((:x, :y))
        end
        all_ages = make_all_ages(scaling_boot, atbp.age_max)
        age_comp = proportion_at_age(scaling_boot, all_ages)
        age_weights = weights_at_age(scaling_boot, surveydata.length_weight, all_ages, 
            bs.weights_at_age)
        @assert ! any(ismissing, age_weights.weight)
        if nrow(age_weights) != atbp.age_max + 1 
            println("nrow=$(nrow(age_weights))")
        end
        ii = trawl_assignments(coordinates.(surveydata.domain), 
                    coordinates.(domain(geotrawl_means)), bs.trawl_assignments)

        nasc = bs.nonneg_lusim ? nonneg_lusim(atbp) : nonneg_lumult(atbp.params, z0)
        cal_error_sim = simulate_cal_error(atbp.cal_error,  bs.calibration)

        df = DataFrame(
            nasc = nasc * cal_error_sim,
            class = atbp.class,
            event_id = trawl_means.event_id[ii]
        )
        df = @chain df begin
            leftjoin(@select(trawl_means, :event_id, :sigma_bs), on=:event_id)
            leftjoin(age_comp, on=:event_id)
            DataFrames.stack(Not([:event_id, :class, :nasc, :sigma_bs]), 
                variable_name=:age, value_name=:p_age)
            DataFramesMeta.@transform(:n_age = :nasc ./ (4π * :sigma_bs) .* :p_age)
            @by(:age, 
                :n_age = sum(skipmissing(:n_age)) * atbp.dA)
            leftjoin(age_weights, on=:age)
            DataFramesMeta.@transform(:biomass_age = :n_age .* :weight, :i = i)
        end
        return df
    end

    return vcat(results...)
end

"""
    simulate_classes(class_problems, surveydata[; nreplicates=500, bs=BootSpecs(),
        exclude_ages=[0]])

Run a bootstrap uncertainty estimation for each scaling class, based on the
`ATBootstrapProblem` definitions in `class_problems` and the data in `surveydata`. The
number of bootstrap replicates can optionally be set in `nreplicates`, and the `BootSpecs`
object `bs` can be used to specify if any of the error sources should be omitted (by
default all are included).
"""
function simulate_classes(class_problems, surveydata; nreplicates=500, bs=BootSpecs(),
        exclude_ages=[0])
    age_strings = lpad.(string.(exclude_ages), 2)
    class_results = map(p -> simulate(p, surveydata; nreplicates, bs), class_problems)
    results = @chain vcat(class_results...) begin
        @subset(.! in(age_strings).(:age))
        @by([:age, :i],
            :n_age = sum(:n_age), :biomass_age = sum(:biomass_age))
    end
    return results
end

"""
    stepwise_errors(class_problems, surveydata[; nreplicates=500, remove=false])

Run a stepwise quantification of each error source for the survey analysis defined by 
`class_problems` and `surveydata`. If `remove=false` (the default), this will add one 
error source at a time, while fixing all others at zero, i.e. treating them as error-free
or deterministic. Use `nreplicates` to set the number of bootstrap replicates.
"""
function stepwise_error(class_problems, surveydata; nreplicates=500, remove=false)
    error_sources = string.(fieldnames(BootSpecs))
    colname = remove ? :eliminated_error : :added_error
    results = map(eachindex(error_sources)) do i
        prefix = remove ? "Omitting" : "Adding"
        println("\n$(prefix) $(error_sources[i]) ($(i)/$(length(error_sources)))...")
        errs = fill(remove, length(error_sources))
        errs[i] = !remove
        bs = BootSpecs(errs...)
        res = simulate_classes(class_problems, surveydata; bs=bs, nreplicates=nreplicates)
        res[!, colname] .= error_sources[i]
        res
    end
    return vcat(results...)
end

# Ordered, categorical labels for stepwise error results
error_labels = DataFrame(
    added_error = ["calibration", "nonneg_lusim", "selectivity", "resample_scaling",
        "drop_trawl", "trawl_assignments", "predict_ts", "age_length", "weights_at_age", "All"],
    error_label = CategoricalArray(
        # ["Calibration", "Spatial sampling", "Trawl jackknife", "Trawl assignment", 
        # "Selectivity", "Resample catches", "Length-weight", "TS models", "Age-length",
        # "All"],
        # levels=["Calibration", "Spatial sampling", "Trawl jackknife", "Trawl assignment", 
        # "Selectivity", "Resample catches", "Length-weight", "TS models", "Age-length",
        # "All"]
        ["Calibration", "Spatial sampling", "Selectivity", "Resample catches", 
        "Trawl dropping", "Trawl assignment", "TS models", "Age-length", "Length-weight", "All"],
        levels=["Calibration", "Spatial sampling", "Selectivity", "Resample catches", 
        "Trawl dropping", "Trawl assignment", "TS models", "Age-length", "Length-weight", "All"]
    )
)


"""
    `read_survey_files(surveydir)`

Convenience function to read all files in a survey data directory `surveydir`. Assumes
that the directory contains files with the following names:
- acoustics_projected.csv
- trawl_locations_projected.csv
- scaling.csv
- age_length.csv
- length_weight.csv
- surveydomain.csv
These files are produced by running `download_survey.R`, which gets the files from
Macebase and saves them to the survey data directory, followed by `preprocess_survey_data`,
which projects everything geographically and averages it into spatial bins.
"""
function read_survey_files(surveydir)
    acoustics = CSV.read(joinpath(surveydir, "acoustics_projected.csv"), DataFrame)
    trawl_locations = CSV.read(joinpath(surveydir, "trawl_locations_projected.csv"), DataFrame)
    scaling = CSV.read(joinpath(surveydir, "scaling.csv"), DataFrame)
    scaling = DataFramesMeta.@transform(scaling, :sample_correction_scalar = float(:sample_correction_scalar))
    age_length = CSV.read(joinpath(surveydir, "age_length.csv"), DataFrame)
    length_weight = CSV.read(joinpath(surveydir, "length_weight.csv"), DataFrame)
    surveydomain = CSV.read(joinpath(surveydir, "surveydomain.csv"), DataFrame)
    surveydomain = DataFrames.shuffle(surveydomain) # this seems to fix the issue with directional artifacts
    surveydomain =  PointSet(Matrix(surveydomain)')
    # return (;acoustics, scaling, age_length, length_weight, trawl_locations, surveydomain)
    return ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations,
        surveydomain)
end

function plot_class_variograms(class_problems; size=(800, 600), kwargs...)
    pp = map(class_problems) do cp
        vg_emp = cp.variogram.empirical
        vg_mod = cp.variogram.model
        plot(vg_emp.abscissa, vg_emp.ordinate, title=cp.class, marker=:o,
            label="Empirical", xlabel="Lag (km)", ylabel="γ", legend=:bottomright)
        plot!(h -> vg_mod(h), 0, maximum(vg_emp.abscissa), label="Model")
    end
    plot(pp...; size=size, kwargs...)
end

function plot_simulated_nasc(atbp::ATBootstrapProblem, surveydata::ATSurveyData,
        simdomain=solution_domain(atbp))
    sim_field = nonneg_lusim(atbp)
    p = scatter(simdomain.x, simdomain.y, zcolor=sim_field, clims=(0, quantile(sim_field, 0.999)), 
        markerstrokewidth=0, markershape=:square, title=string(atbp.class),
        aspect_ratio=:equal, markersize=2.2, legend=false,
        xlabel="Easting (km)", ylabel="Northing (km)")
    df = @subset(surveydata.acoustics, :class .== atbp.class)
    scatter!(p, df.x, df.y, color=:white, markersize=df.nasc*3e-3, alpha=0.3,
        markerstrokewidth=0)
    return p
end

function plot_simulated_nasc(class_problems::Vector{<:ATBootstrapProblem},
        surveydata::ATSurveyData, simdomain = solution_domain(first(class_problems));
        kwargs...)
    plots = [plot_simulated_nasc(p, surveydata, simdomain) for p in class_problems]
    return plot(plots...; kwargs...)
end

function plot_boot_results(results; size=(900, 400), margin=15px, palette=:Paired_10,
        kwargs...)
    p_abundance = @df results violin(:age, :n_age/1e9, group=:age, palette=palette,
    xlabel="Age class", ylabel="Abundance (billions)", legend=false);
    p_biomass = @df results violin(:age, :biomass_age/1e9, group=:age, palette=palette,
        xlabel="Age class", ylabel="Biomass (Mt)");
    plot(p_abundance, p_biomass; size=size, margin=margin, kwargs...)
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

end # module