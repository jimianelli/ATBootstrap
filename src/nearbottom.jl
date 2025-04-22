# using DataFrames
# using Statistics
# using StatsPlots
using Distributions

species = CSV.read(joinpath(@__DIR__, "..", "surveydata", "species.csv"), DataFrame)
species.nearbottom_group .= "Misc"
species.nearbottom_group[species.species_code .== 21740] .= "Pollock"
species.nearbottom_group[species.species_code .== 21725] .= "Arctic cod"
species.nearbottom_group[10110 .<= species.species_code .<= 10120] .= "Large flatfish"
species.nearbottom_group[30050 .<= species.species_code .<= 30535] .= "Rockfish"


nearbottom_coefs = DataFrame(
    nearbottom_group = ["Pollock", "Arctic cod", "Large flatfish", "Rockfish", "Misc"],
    A = [2.52, 16.39, 0.85, 93.59, 11.63],
    lower = [2.21, 2.84, 0.16, 8.63, 3.92],
    upper = [2.86, 62.26, 1.67, 343.43, 22.09],
    b = [66, 67.4, 67.4, 67.4, 67.4]
)

nearbottom_coefs = DataFramesMeta.@transform(nearbottom_coefs,
    :sd_A = ((:A .- :lower) .+ (:upper .- :A)) / 2 / 2,
)
const nearbottom_intercept = 3.43
# bar(nearbottom_coefs.nearbottom_group, nearbottom_coefs.a)

nearbottom_coefs = leftjoin(species, nearbottom_coefs, on=:nearbottom_group)

function make_nearbottom_dict(stochastic=true)
    if stochastic
        dists = truncated.(Normal.(nearbottom_coefs.A, nearbottom_coefs.sd_A),
            # 0, Inf)
        nearbottom_coefs.lower, nearbottom_coefs.upper)
        aa = rand.(dists)
    else 
        aa = nearbottom_coefs.A
    end
    aa = aa ./ exp10.(-nearbottom_coefs.b/10) / 1852^2 / 4π
    return Dict(zip(nearbottom_coefs.species_code, aa)) 
end

function apply_nearbottom_coefficient!(scaling, nearbottom_dict)
    @eachrow! scaling begin
        if :species_code in keys(nearbottom_dict)
            :user_defined_expansion = nearbottom_dict[:species_code] 
        else
            # if species code isn't in dict, it's speckled eelpout (i.e.,"Misc")
            :user_defined_expansion = nearbottom_dict[24166]
        end
        :w = :catch_sampling_expansion * :user_defined_expansion * :sample_correction_scalar * :haul_weight
    end
end
