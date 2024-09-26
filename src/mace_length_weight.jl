
# function make_weight_function(stochastic=false, sd=0.05)
#     err = stochastic ? (1 + sd*randn()) : 1.0
#     predict_weight(L) = (1e-5 * L^2.9) * err
#     return predict_weight
# end


using GLM
using StatsBase

function make_weight_function(length_weight, stochastic=true, nmin=5)
    n = nrow(length_weight)
    ii = stochastic ? sample(1:n, n) : 1:n
    length_weight_boot = @view length_weight[ii, :]
    model = lm(@formula(log(organism_weight) ~ log(fork_length)), length_weight_boot)
    c = DataFrame(coeftable(model))
    μloga = c[1,2]
    μb = c[2,2]
    a = exp(μloga)
    b =  μb
    Lmax = round(Int, maximum(length_weight.fork_length))
    all_lengths = DataFrame(fork_length = 1:Lmax)   
    binned = @chain length_weight_boot begin
        rightjoin(all_lengths, on=:fork_length)
        @by(:fork_length, 
            :mean_weight = mean(:organism_weight),
            :n = length(:organism_weight))
        DataFramesMeta.@transform(
            :weight = ifelse.(:n .> nmin, :mean_weight, a.*:fork_length.^b)
        )
    end
    weight_dict = Dict(zip(binned.fork_length, binned.weight))

    function weight_function(L)
        L1 = min(max(round(L), 1.0), Lmax)
        return weight_dict[L1]
    end
    
    return weight_function 
end

# wf = make_weight_function(length_weight)
# wf1 = make_weight_function(length_weight, false)
# wf(40)
# wf1(40)

# m = lm(@formula(log(organism_weight) ~ log(fork_length)), length_weight)
# c = DataFrame(coeftable(m))
# a = exp(c[1,2])
# b = c[2,2]

# f(L) = a * L^b

# nmin = 10
# all_lengths = DataFrame(fork_length = 1:maximum(length_weight.fork_length))
# binned = @chain length_weight begin
#     rightjoin(all_lengths, on=:fork_length)
#     @by(:fork_length, 
#         :mean_weight = mean(:organism_weight),
#         :n = length(:organism_weight))
#     DataFramesMeta.@transform(:weight = ifelse.(:n .> nmin, :mean_weight, f.(:fork_length)))
# end

# @df length_weight scatter(:fork_length, :organism_weight, markerstrokewidth=0, alpha=0.1,
#     xlabel="Fork length (cm)", ylabel="Weight (kg)", label="Measurements")
# @df binned scatter!(:fork_length, :mean_weight, label="Binned average")
# @df binned scatter!(:fork_length, :weight, label="Average, filled in")
# plot!(L -> a * L^b, 10, 69, label="Allometric model")
# plot!(wf, 10, 69, label="Combined")

# scatter(length_weight.fork_length, residuals(m), markerstrokewidth=0, alpha=0.1)

