using CSV, DataFrames, DataFramesMeta
using Statistics, StatsBase
using StableRNGs

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

survey = get(ENV, "ATB_SURVEY", "202207")
nsim = parse(Int, get(ENV, "ATB_SPATIAL_NSIM", "40"))
out_path = get(
    ENV,
    "ATB_SPATIAL_PROBE_OUT",
    joinpath(@__DIR__, "..", "docs", "artifacts", "spatial_backend_julia.csv"),
)

StableRNGs.seed!(parse(Int, survey))

surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0
dA = (resolution * km2nmi)^2
scaling_classes = ["SS1", "SS2", "BT"]
lag_points = collect(0.0:20.0:200.0)

(; acoustics, scaling, age_length, length_weight, trawl_locations, surveygrid) = ATB.read_survey_files(surveydir)
acoustics = @subset(acoustics, in(scaling_classes).(:class))

surveydata = ATB.ATSurveyData(
    acoustics,
    scaling,
    age_length,
    length_weight,
    trawl_locations,
    surveygrid,
    dA,
)

atbp = ATB.ATBootstrapProblem(surveydata; scaling_classes=scaling_classes)

det_bs = ATB.BootSpecs(
    selectivity=false,
    predict_ts=false,
    resample_scaling=false,
    nearbottom_coefs=false,
    age_length=false,
    weights_at_age=false,
    trawl_assignments=false,
    simulate_nasc=false,
    calibration=false,
)

spatial_only_bs = ATB.BootSpecs(
    selectivity=false,
    predict_ts=false,
    resample_scaling=false,
    nearbottom_coefs=false,
    age_length=false,
    weights_at_age=false,
    trawl_assignments=false,
    simulate_nasc=true,
    calibration=false,
)

metrics = map(atbp.class_problems) do cp
    obs_nasc = @subset(surveydata.acoustics, :class .== cp.class).nasc
    det_field = ATB.nonneg_lumult(cp.params, mean.(cp.zdists))
    sim_fields = [ATB.simulate_nasc(cp) for _ in 1:nsim]
    sim_means = mean.(sim_fields)
    sim_sds = std.(sim_fields)

    det_res = ATB.simulate_class_iteration(cp, surveydata, det_bs, 1)
    sim_res = [ATB.simulate_class_iteration(cp, surveydata, spatial_only_bs, i) for i in 1:nsim]
    sim_total_n = [sum(df.n) for df in sim_res]
    sim_total_biomass = [sum(df.biomass) for df in sim_res]

    row = DataFrame(
        backend = "julia",
        survey = survey,
        class = cp.class,
        zdist = string(cp.zfamily),
        obs_mean_nasc = mean(obs_nasc),
        obs_sd_nasc = std(obs_nasc),
        field_det_mean = mean(det_field),
        field_det_sd = std(det_field),
        field_sim_mean_mean = mean(sim_means),
        field_sim_mean_sd = std(sim_means),
        field_sim_sd_mean = mean(sim_sds),
        field_sim_sd_sd = std(sim_sds),
        det_total_n = sum(det_res.n),
        det_total_biomass = sum(det_res.biomass),
        sim_total_n_mean = mean(sim_total_n),
        sim_total_n_sd = std(sim_total_n),
        sim_total_biomass_mean = mean(sim_total_biomass),
        sim_total_biomass_sd = std(sim_total_biomass),
    )

    for lag in lag_points
        row[!, Symbol("variogram_gamma_", Int(round(lag)))] = [cp.variogram.model(lag)]
    end

    row
end

mkpath(dirname(out_path))
CSV.write(out_path, vcat(metrics...))
println("Wrote Julia spatial backend metrics to ", out_path)
