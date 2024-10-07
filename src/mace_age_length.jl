function nearest_length(L, unique_lengths)
    d, i = findmin(x -> abs(L - x), unique_lengths)
    return i
end

function make_age_length_function(age_length, age_max=10, stochastic=true)

    key = @chain age_length begin
        DataFramesMeta.@transform(
            :fork_length = round.(:fork_length),
            :age = ifelse.(:age .> age_max, age_max, :age)
        )
        @by([:fork_length, :age], :n = length(:age))
        unstack(:age, :n, fill=0)
    end

    arr = Array(key[:, 2:end])
    arr = arr ./ sum(arr, dims=2)
    dists = [Distributions.Categorical(arr[i, :]) for i in 1:size(arr, 1)]
    # length_tree = KDTree(key.fork_length')

    func = stochastic ? rand : mode
    predict_age = let fl = key.fork_length
        L -> begin
            # i, _ = nn(length_tree, [L])
            i = nearest_length(L, fl)
            return func(dists[i])
        end
    end
    return predict_age
end
