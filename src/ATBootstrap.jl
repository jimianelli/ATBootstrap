module ATBootstrap

using CSV, DataFrames, DataFramesMeta, CategoricalArrays
using Random
using GeoStats, Geodesy, ConcaveHull
using Statistics, StatsBase
using Distributions
using LinearAlgebra
using Distances
using NearestNeighbors
using StaticArrays
using ProgressMeter
using RCall
using StatsPlots, StatsPlots.PlotMeasures
using ColorSchemes

export preprocess_survey_data,
    AbstractSurveyDomain,
    TransectRibbons,
    SurveyHull,
    get_survey_grid,
    read_survey_files,
    in_intervals,
    ATSurveyData,
    ScalingClassProblem,
    ATBootstrapProblem,
    resample_df,
    nonneg_lusim,
    nonneg_lusim!,
    nonneg_lumult,
    nonneg_lumult!,
    zdists,
    simulate_nasc,
    solution_domain,
    simulate_class,
    simulate,
    stepwise_error,
    BootSpecs,
    error_labels,
    plot_class_variograms,
    plot_simulated_nasc,
    plot_geosim_stats,
    plot_boot_results,
    plot_error_sources,
    plot_error_source_by_age,
    summarize_bootstrap,
    summarize_stepwise_bootstrap,
    merge_results

include("types.jl")
include("preprocess_survey_data.jl")
include("spatial.jl")
include("mace_ts.jl")
include("mace_age_length.jl")
include("mace_selectivity.jl")
include("mace_length_weight.jl")
include("calibration.jl")
include("scaling.jl")
include("display.jl")

"""
    Extract the bare x and y values from a `Point` as an `SVector`.
"""
svector_coords(pt::Point) = SVector(pt.coords.x.val, pt.coords.y.val)

function simulate_class_iteration(scp::ScalingClassProblem, surveydata::ATSurveyData,
        bs=BootSpecs(), i=1, scaling_sub = @subset(surveydata.scaling, :class .== scp.class),
        surveygrid_coords = svector_coords.(surveydata.grid),
        z0 = mean.(scp.zdists))

    selectivity_function = make_selectivity_function(bs.selectivity)
    scaling_boot = resample_scaling(scaling_sub, bs.resample_scaling)
    apply_selectivity!(scaling_boot, selectivity_function)

    predict_ts = make_ts_function()
    predict_age = make_age_length_function(surveydata.age_length, scp.age_max,
        bs.age_length)
    
    scaling_boot = DataFramesMeta.@transform!(scaling_boot,
        :sigma_bs = to_linear.(predict_ts.(:ts_relationship, :ts_length, bs.predict_ts)),
        :age = predict_age.(:primary_length))
    
    geotrawls = @chain scaling_boot begin
        @by(:event_id, :class = first(:class))
        innerjoin(surveydata.trawl_locations, on=:event_id)
        @select(:x, :y)
        georef((:x, :y))
    end

    all_ages = make_all_ages(scaling_boot, scp.age_max)
    age_weights = pollock_weights_at_age(scaling_boot, surveydata.length_weight,
        all_ages, bs.weights_at_age)
    
    ii = trawl_assignments(surveygrid_coords, 
                svector_coords.(domain(geotrawls)), bs.trawl_assignments)

    nasc = bs.simulate_nasc ? simulate_nasc(scp) : nonneg_lumult(scp.params, z0)
    cal_error_sim = simulate_cal_error(scp.cal_error,  bs.calibration)
    nasc_df = DataFrame(
        nasc = nasc * cal_error_sim,
        event_id = surveydata.trawl_locations.event_id[ii]
    )

    # Special-case processing for bottom-trawl stratum
    if scp.class == "BT"
        # weight scaling data with Nate's nearbottom coefficients
        nearbottom_dict = make_nearbottom_dict(bs.nearbottom_coefs)
        apply_nearbottom_coefficient!(scaling_boot, nearbottom_dict)
        # If species is not pollock, make TS deterministic (variability in backscatter
        # allocation has already been accounted for with nearbottom coefficients)
        scaling_boot = DataFramesMeta.@transform(scaling_boot,
            @byrow :sigma_bs = :species_code == 21740 ? 
                :sigma_bs : to_linear(predict_ts(:ts_relationship, :ts_length, false))
        )
        # remove nearbottom intercept from nasc (from Nate's paper)
        # nasc_df.nasc .-= nearbottom_intercept
        nasc_df.nasc .= max.(nasc_df.nasc, 0) # make sure we don't end up with negative backscatter
    end

    trawl_means_cat = get_trawl_category_means(scaling_boot, scp.aged_species)
    category_nasc = unstack(trawl_means_cat, :event_id, :category, :p_nasc, fill=0)
    category_sigma = @select(trawl_means_cat, :event_id, :category, :sigma_bs)

    df = @chain nasc_df begin
        leftjoin(category_nasc, on=:event_id)
        DataFrames.stack(Not([:event_id, :nasc]),
            variable_name=:category, value_name=:p_nasc)
        leftjoin(category_sigma, on=[:event_id, :category])
        dropmissing()
        DataFramesMeta.@transform(:n = :nasc .* :p_nasc ./ (4π .* :sigma_bs))
        @by(:category,
            :n = sum(skipmissing(:n)) * surveydata.dA)
        DataFramesMeta.@transform(
            :species_code = parse.(Int, first.(split.(:category, "@"))),
            :age = parse.(Int, last.(split.(:category, "@")))
        )
        leftjoin(age_weights, on=[:species_code, :age])
        DataFramesMeta.@transform(
            :biomass = :n .* :weight, 
            :i = i
        )
        @select(:i, :species_code, :age, :category, :n, :biomass)
        @orderby(:i, :species_code, :age, :n, :biomass)
    end
    return df
end


"""
    simulate_class(scp, surveydata; nreplicates=500, bs=BootSpecs())

Do a bootstrap analysis for one acoustic class/scaling stratum using the problem defined
in `scp` and the data in `surveydata`. The number of bootstrap replicates can optionally
be set in `nreplicates`, and the `BootSpecs` object `bs` can be used to specify if any of
the error sources should be omitted (by default all are included).
"""
function simulate_class(scp::ScalingClassProblem, surveydata::ATSurveyData; nreplicates=500, 
        bs=BootSpecs())

    scaling_sub = @subset(surveydata.scaling, :class .== scp.class)
    surveygrid_coords = svector_coords.(surveydata.grid)
    z0 = mean.(scp.zdists)

    println("Bootstrapping $(scp.class)...")
    results = @showprogress map(1:nreplicates) do i
        simulate_class_iteration(scp, surveydata, bs, i, scaling_sub, surveygrid_coords, z0)
    end

    return vcat(results...)
end

"""
    simulate(atbp, surveydata[; nreplicates=500, bs=BootSpecs(),
        report_species=[21740], report_ages=1:first(atbp.class_problems).age_max)

Run a bootstrap uncertainty estimation for each scaling class, based on the
`ATBootstrapProblem` definition in `atbp` and the data in `surveydata`. The
number of bootstrap replicates can optionally be set in `nreplicates`, and the `BootSpecs`
object `bs` can be used to specify if any of the error sources should be omitted (by
default all are included).

By default, abundance and biomass are reported for pollock by age class. Abundances are
calculated for other species (but not biomasses, and not by age, since this code 
doesn't put together the length-age and length-weight tables for them). If you would like
to return the abundance results for other species, you can list their species codes in the
optional `report_species` argument. For *all* species and ages, pass in empty arrays (`[]`)
for `report_species` and `report_ages`.

For pollock, you can limit the ages reported by passing them to the `report_ages` argument.
Note that this just subsets the results post-simulation; to set the maximum age see the 
documentation for `ScalingClassProblem.`
"""
function simulate(atbp::ATBootstrapProblem, surveydata::ATSurveyData; nreplicates=500,
        bs=BootSpecs(), report_species=[21740], report_ages=1:first(atbp.class_problems).age_max)
    
    class_results = map(p -> simulate_class(p, surveydata; nreplicates, bs), atbp.class_problems)
    if length(report_species) == 0
        report_species = unique(surveydata.scaling.species_code)
    end
    if length(report_ages) == 0
        report_ages = [-1; 1:first(atbp.class_problems).age_max]
    end
    results = @chain vcat(class_results...) begin
        @subset(
            in(report_ages).(:age),
            in(report_species).(:species_code)
        )
        @by([:i, :species_code, :age], :n = sum(:n), :biomass = sum(:biomass))
        DataFramesMeta.@transform(:age = replace(:age, -1 => missing))
    end
    return results
end

"""
    stepwise_errors(atbp, surveydata[; nreplicates=500, remove=false])

Run a stepwise quantification of each error source for the survey analysis defined by 
`atbp` and `surveydata`. If `remove=false` (the default), this will add one 
error source at a time, while fixing all others at zero, i.e. treating them as error-free
or deterministic. Use `nreplicates` to set the number of bootstrap replicates.
"""
function stepwise_error(atbp, surveydata; nreplicates=500, remove=false)
    error_sources = string.(fieldnames(BootSpecs))
    colname = remove ? :eliminated_error : :added_error
    results = map(eachindex(error_sources)) do i
        prefix = remove ? "Omitting" : "Adding"
        println("\n$(prefix) $(error_sources[i]) ($(i)/$(length(error_sources)))...")
        errs = fill(remove, length(error_sources))
        errs[i] = !remove
        bs = BootSpecs(errs...)
        res = simulate(atbp, surveydata; bs=bs, nreplicates=nreplicates)
        res[!, colname] .= error_sources[i]
        res
    end
    return vcat(results...)
end

# Ordered, categorical labels for stepwise error results
error_labels = DataFrame(
    added_error = ["calibration", "simulate_nasc", "selectivity", "resample_scaling",
        "drop_trawl", "trawl_assignments", "predict_ts", "age_length", "weights_at_age", "All"],
    error_label = CategoricalArray(
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
- surveygrid.csv
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
    surveygrid = CSV.read(joinpath(surveydir, "surveygrid.csv"), DataFrame)
    surveygrid = DataFrames.shuffle(surveygrid) # this seems to fix the issue with directional artifacts
    surveygrid =  PointSet(Point(x...) for x in eachrow(surveygrid))
    return (;acoustics, scaling, age_length, length_weight, trawl_locations, surveygrid)
end

function in_intervals(x::Real, intervals)
    for int in intervals
        if int[1] <= x <= int[2]
            return true
        end
    end
    return false
end

function in_intervals(x::AbstractVector, intervals)
    return [in_intervals(xi, intervals) for xi in x]
end


end # module