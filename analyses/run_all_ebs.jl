
ebs_surveys = [
    "200707",
    "200809",
    "200909",
    "201006",
    "201207",
    "201407",
    "201608",
    "201807",
    "202207"
]

for survey in ebs_surveys
    println("Running analysis for survey $(survey)")
    include(joinpath(@__DIR__, "dy$(survey).jl"))
end
