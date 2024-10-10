#=
Type definitions for data structures
=#

"""
    ATSurveyData(acoustics, scaling, age_length, length_weight, trawl_locations,
        domain, dA)

Construct an `ATSurveyData` object, encapsulating all the relevant data 
"""
struct ATSurveyData
    acoustics::DataFrame
    scaling::DataFrame
    age_length::DataFrame
    length_weight::DataFrame
    trawl_locations::DataFrame
    domain::Meshes.Domain
    dA::Number
end

"""
    BootSpecs([; selectivity, predict_ts, resample_scaling, drop_trawl, age_length,
        weights_at_age, trawl_assignments, simulate_nasc, calibration])
    
Construct a `BootSpecs` object, which specifies which error sources to include when running
`simulate` on an `ATBootstrapProblem`. All arguments are Booleans and default to true; they
can be 

# Examples

julia> BootSpecs(calibration=false) # turn off calibration error but include everything else

julia> BootSpecs(false) # turn off all errors--i.e., do a normal deterministic analysis

"""
@kwdef struct BootSpecs
    selectivity::Bool=true
    predict_ts::Bool=true
    resample_scaling::Bool=true
    drop_trawl::Bool=true
    age_length::Bool=true
    weights_at_age::Bool=true
    trawl_assignments::Bool=true
    simulate_nasc::Bool=true
    calibration::Bool=true
end
BootSpecs(b::Bool) = BootSpecs(fill(b, length(fieldnames(BootSpecs)))...)

struct ScalingClassProblem
    class
    variogram
    geosetup
    params
    zfamily
    zdists
    cal_error
    age_max
end

"""
    ScalingClassProblem(surveydata, class[, cal_error=0.1, age_max=10,
        zdist_candidates=zdist_candidates, maxlag=200.0, nlags=20, weightfunc=h -> 1/h])

Set up a `ScalingClassProblem`, specifying how to perform geostatistical simulations of
backscatter in scaling stratum `class`, conditional on the observed NASC values in 
`surveydata`. 
"""
function ScalingClassProblem(surveydata, class; 
        cal_error=0.1, age_max=10, maxlag=200.0, nlags=20, weightfunc=h -> 1/h,
        zdist_candidates=[Gamma, InverseGamma, InverseGaussian, LogNormal])
    acoustics_sub = @subset(surveydata.acoustics, :class .== class)
    variogram, geosetup = define_conditional_sim(acoustics_sub, surveydata.domain,
        maxlag=maxlag, nlags=nlags, weightfunc=weightfunc)
    params = get_lungs_params(geosetup, variogram.model)
    optimal_dist = choose_z_distribution(zdist_candidates, acoustics_sub.nasc, params)
    zdists = parameterize_zdists(optimal_dist, params)
    return ScalingClassProblem(class, variogram, geosetup, params, optimal_dist, zdists,
        cal_error, age_max)
end

"""
Extract the simulation domain from a `ScalingClassProblem`. Returns a `DataFrame` with two
columns containing the `x` and `y` coordinates of each point at which the spatial field 
is to be simulated.
"""
function solution_domain(scp::ScalingClassProblem, variable=:nasc)
    # sol = solve(scp.geosetup, LUGS(variable => (variogram = scp.variogram.model,)))
    # dom = domain(sol)
    dom = scp.geosetup.domain
    x = [p.coords.x for p in dom]
    y = [p.coords.y for p in dom]
    return DataFrame(x=x, y=y)
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
        println("Preparing $(class)...")
        return ScalingClassProblem(surveydata, class, maxlag=maxlag, nlags=nlags,
            age_max=age_max, cal_error=cal_error, weightfunc=weightfunc, 
            zdist_candidates=zdist_candidates)
    end
    return ATBootstrapProblem(class_problems, scaling_classes)
end