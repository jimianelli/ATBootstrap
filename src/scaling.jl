# Numbers from Matta and Kimura 2012, Age determination manual
const L∞ = 67.33 # mm
const t₀ = -0.205
const K = 0.1937
const AGE_MAX = 10
predict_length(t) = L∞ * (1 - exp(-K * (t - t₀)))
function predict_age_deterministic(L, age_max=AGE_MAX)
    if L < predict_length(age_max)
        return round(Int, log(1 - L/L∞) / -K + t₀)
    else
       return age_max
   end
end

function predict_age_stochastic(L, age_max=AGE_MAX)
    return max(0, predict_age_deterministic(L + 2randn(), age_max))# + rand([-1, 0, 0, 0, 1]))
end

function predict_age(L, stochastic=true, age_max=AGE_MAX)
    if stochastic
        return predict_age_stochastic(L, age_max)
    else
        return predict_age_deterministic(L, age_max)
    end
end



function resample_df(df, stochastic=true)
    n = nrow(df)
    if stochastic
        ii = sample(1:n, n)
    else
        ii = 1:n
    end
    return @view df[ii, :]
end

function resample_scaling(df, stochastic=true)
    return DataFramesMeta.combine(x -> resample_df(x, stochastic),
        DataFramesMeta.groupby(df, [:event_id, :class]))
end

function get_trawl_means(scaling, trawl_locations)
    trawl_means = @chain scaling begin
        @by(:event_id, 
            :sigma_bs = mean(:sigma_bs, Weights(:w)),
            :sigma_bs_std = std(:sigma_bs, Weights(:w)),
            :length = mean(:primary_length, Weights(:w)))        
            DataFramesMeta.@transform(:ts = 10log10.(:sigma_bs))
        innerjoin(trawl_locations, on=:event_id)
    end
    return trawl_means
end

function make_all_ages(scaling, age_max)
    return allcombinations(DataFrame, 
        event_id = unique(scaling.event_id), 
        age = 0:age_max)
end

function weights_at_age(scaling, length_weight, all_ages, stochastic=false)
    predict_weight = make_weight_function(length_weight, stochastic)
    res = @chain scaling begin
        DataFramesMeta.@transform(
            :weight = predict_weight.(:primary_length))
        rightjoin(all_ages, on=[:event_id, :age])
        DataFramesMeta.@transform(:weight = replace(:weight, missing => 0.0))
        DataFramesMeta.@transform(:age = lpad.(string.(:age), 2, "0"))
        @by(:age, :weight = mean(:weight))
    end
    return res
end

function proportion_at_age(scaling, all_ages; age_max=AGE_MAX)#, stochastic=false)
    age_comp = @chain scaling begin
        # DataFramesMeta.@transform(:age = predict_age.(:primary_length, stochastic))
        @by([:event_id, :age], :p_age=sum(:w))
        DataFrames.groupby(:event_id)
        DataFramesMeta.@transform(:p_age = :p_age / sum(:p_age))
        rightjoin(all_ages, on=[:event_id, :age])
        DataFramesMeta.@transform(:p_age = replace(:p_age, missing => 0.0))
        DataFramesMeta.@transform(:age = lpad.(string.(:age), 2, "0"))
        @orderby(:event_id, :age)
    end
    # return age_comp
    return DataFrames.unstack(age_comp, :age, :p_age)
end
