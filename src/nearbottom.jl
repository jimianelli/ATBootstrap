# using DataFrames
# using Statistics
# using StatsPlots

species = CSV.read(joinpath(@__DIR__, "..", "surveydata", "species.csv"), DataFrame)
species.nearbottom_group .= "Misc"
species.nearbottom_group[species.species_code .== 21740] .= "Pollock"
species.nearbottom_group[species.species_code .== 21725] .= "Arctic cod"
species.nearbottom_group[10110 .<= species.species_code .<= 10120] .= "Large flatfish"
species.nearbottom_group[30050 .<= species.species_code .<= 30535] .= "Rockfish"


nearbottom_coefs = DataFrame(
    nearbottom_group = ["Pollock", "Arctic cod", "Large flatfish", "Rockfish", "Misc"],
    A = [2.52, 16.39, 0.85, 93.59, 11.63],
    b = [66, 67.4, 67.4, 67.4, 67.4]
)

nearbottom_coefs.a = nearbottom_coefs.A ./ exp10.(-nearbottom_coefs.b/10) ./ 1852^2
const nearbottom_intercept = 3.43
# bar(nearbottom_coefs.nearbottom_group, nearbottom_coefs.a)

nearbottom_coefs = @chain leftjoin(species, nearbottom_coefs, on=:nearbottom_group) begin
    @select(:species_code, :user_defined_expansion = :a)
end

