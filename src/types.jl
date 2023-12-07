#=
Type definitions for data structures
=#

struct ATSurveyData
    acoustics
    scaling
    age_length
    length_weight
    trawl_locations
    domain
    dA
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

struct ScalingClassProblem
    class
    variogram
    problem
    params
    zfamily
    zdists
    cal_error
    age_max
end

"""
    ScalingClassProblem(surveydata, class[, cal_error=0.1, age_max=10,
        zdist_candidates=zdist_candidates, maxlag=200.0, nlags=20, weightfunc=h -> 1/h])

Set up a `ScalingClassProblem`, specifying how to simulate bootstrap analyses of the 
scaling stratum `class` from the data in `surveydata`.
"""
function ScalingClassProblem(surveydata, class; 
        cal_error=0.1, age_max=10, maxlag=200.0, nlags=20, weightfunc=h -> 1/h,
        zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal])
    acoustics_sub = @subset(surveydata.acoustics, :class .== class)
    variogram, problem = define_conditional_sim(acoustics_sub, surveydata.domain,
        maxlag=maxlag, nlags=nlags, weightfunc=weightfunc)
    params = get_lungs_params(problem, variogram.model)
    optimal_dist = choose_z_distribution(zdist_candidates, acoustics_sub.nasc, params)
    zdists = parameterize_zdists(optimal_dist, params)
    return ScalingClassProblem(class, variogram, problem, params, optimal_dist, zdists,
        cal_error, age_max)
end

struct ATBootstrapProblem{TP<:ScalingClassProblem, TS<:AbstractString}
    class_problems::Vector{TP}
    scaling_classes::Vector{TS}
end

"""
    ATBootstrapProblem(surveydata, scaling_classes[; cal_error=0.1, age_max=10,
        zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal],
        maxlag=200, nlags=10, weightfunc=h -> 1/h])

Set up an `ATBootstrapProblem`, describing how to do bootstrap analyses of the 
acoustic-trawl survey recorded in `surveydata` for the scaling strata specified in 
`scaling_classes`.
"""
function ATBootstrapProblem(surveydata::ATSurveyData, scaling_classes::Vector{<:AbstractString};
        cal_error=0.1, age_max=10,  
        zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal],
        maxlag=200.0, nlags=10, weightfunc=h -> 1/h)
    
    class_problems = map(scaling_classes) do class
        println(class)
        return ScalingClassProblem(surveydata, class, maxlag=maxlag, nlags=nlags,
            age_max=age_max, cal_error=cal_error, weightfunc=weightfunc, 
            zdist_candidates=zdist_candidates)
    end
    return ATBootstrapProblem(class_problems, scaling_classes)
end