using CSV, DataFrames, DataFramesMeta
using LinearAlgebra
using Statistics, StatsBase
using StableRNGs

include(joinpath(@__DIR__, "..", "src", "ATBootstrap.jl"))
import .ATBootstrap as ATB

survey = get(ENV, "ATB_SURVEY", "202207")
out_path = get(
    ENV,
    "ATB_SPATIAL_CONDITIONAL_OUT",
    joinpath(@__DIR__, "..", "docs", "artifacts", "spatial_conditional_julia.csv"),
)

StableRNGs.seed!(parse(Int, survey))

surveydir = joinpath(@__DIR__, "..", "surveydata", survey)
const km2nmi = 1 / 1.852
resolution = 10.0
dA = (resolution * km2nmi)^2
scaling_classes = ["SS1", "SS2", "BT"]

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

head_string(x, n=8) = isempty(x) ? "" : join(first(collect(x), min(length(x), n)), " ")

function summary_stats(x, prefix)
    x = collect(Float64, x)
    DataFrame(
        Symbol(prefix * "_n") => [length(x)],
        Symbol(prefix * "_mean") => [mean(x)],
        Symbol(prefix * "_sd") => [std(x)],
        Symbol(prefix * "_min") => [minimum(x)],
        Symbol(prefix * "_max") => [maximum(x)],
        Symbol(prefix * "_nonpos_n") => [count(<=(0), x)],
    )
end

function derive_mu_z(params; ϵ=cbrt(eps()))
    shift = copy(params.μx)
    shift[shift .< ϵ] .= ϵ
    target_mean = mean(params.data)
    if !isfinite(target_mean) || target_mean <= 0
        target_mean = mean(shift)
    end
    scaled_shift = shift .* target_mean ./ mean(shift)
    μz = params.L \ scaled_shift
    μz[μz .<= ϵ] .= ϵ
    (; scaled_shift, μz)
end

metrics = map(atbp.class_problems) do cp
    params = cp.params
    obs_nasc = @subset(surveydata.acoustics, :class .== cp.class)
    derived = derive_mu_z(params)
    L = Matrix(params.L)
    ldiag = diag(L)
    implied_var = vec(sum(abs2, L; dims=2))

    hcat(
        DataFrame(
            backend = "julia",
            survey = survey,
            class = cp.class,
            domain_n = length(params.dlocs) + length(params.slocs),
            obs_n = nrow(obs_nasc),
            mapped_dloc_n = length(params.dlocs),
            mapped_collision_n = nrow(obs_nasc) - length(params.dlocs),
            params_data_n = length(params.data),
            params_dloc_n = length(params.dlocs),
            params_sloc_n = length(params.slocs),
            dloc_head = head_string(params.dlocs),
            sloc_head = head_string(params.slocs),
            data_mean = mean(params.data),
            data_sd = std(params.data),
            mu = params.μ,
        ),
        summary_stats(params.μx, "prep_shift"),
        summary_stats(derived.scaled_shift, "scaled_shift"),
        summary_stats(derived.μz, "mu_z"),
        summary_stats(ldiag, "L_diag"),
        summary_stats(implied_var, "implied_var"),
    )
end

mkpath(dirname(out_path))
CSV.write(out_path, vcat(metrics...))
println("Wrote Julia spatial conditional metrics to ", out_path)
