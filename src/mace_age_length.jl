# const AGE_MAX = 10

function make_age_length_function(age_length, age_max=10, stochastic=true)
    #=
    Means, standard deviations, and extrema for age classes 1 through 25. Copied from 
    MACE Applications\MacebaseAnalysis\AL_relationships\data_matrix_gaussian_summer_ebs_2018.py
    =#
    μa=[15.185,24.67,32.525,39.263,43.459,46.961,48.92,51.052,53.433,55.359,56.938,57.871,
        59.854,61.209,60.487,62.529,62.667,70.778,59.5,68.5,73.5,67.333,67.333,67.333,67.333]
    σa=[1.8679,2.8813,3.3516,3.1506,3.1858,3.3012,3.8953,4.4303,4.9431,4.9798,4.987,5.3796,
        5.9018,5.52,5.6787,7.4417,3.7238,7.4963,4.2032,14.849,7.7782,11.59,11.59,11.59,11.59]
    min_a=[8,15,19,22,31,34,39,40,41,43,45,46,46,51,53,49,57,60,54,58,68,54,54,54,54]
    max_a=[21,37,49,54,56,60,66,66,70,69,73,76,76,74,75,79,67,82,64,79,79,75,75,75,75]
    dists = truncated.(Normal.(μa, σa), min_a, max_a)
    all_ages = DataFrame(age = 1:25)

    age_length = stochastic ? resample_df(age_length) : age_length
    p_age = @chain age_length begin
        @by(:age, :p = length(:length))
        DataFramesMeta.@transform(:p = :p / sum(:p))
        rightjoin(all_ages, on=:age)
        DataFramesMeta.@transform(:p = replace(:p, missing => 0))
        select(:p)
        Array()
        vec()
    end

    if stochastic 
        predict_age = L -> begin
            A = sample(1:25, Weights(pdf.(dists, L) .* p_age))
            return min(A, age_max)
        end
        return predict_age
    else
        predict_age = L -> begin
            A = argmax(pdf.(dists, L) .* p_age)
            return min(A, age_max)
        end
        return predict_age
    end
end


# pred_age = make_age_length_function(age_length, 10, true)
# pred_age(47.2)

# age_length = CSV.read(joinpath(surveydir, "age_length.csv"), DataFrame)

# key = unstack(age_length, :length, :age, :proportion, fill=0.0)
# heatmap(Array(key[:, 2:end])', xlabel="Length (cm)", ylabel="Age")
# hline!([10])


# μh=[15.185,24.67,32.525,39.263,43.459,46.961,48.92,51.052,53.433,55.359,56.938,57.871,
#     59.854,61.209,60.487,62.529,62.667,70.778,59.5,68.5,73.5,67.333,67.333,67.333,67.333]
# σh=[1.8679,2.8813,3.3516,3.1506,3.1858,3.3012,3.8953,4.4303,4.9431,4.9798,4.987,5.3796,
#     5.9018,5.52,5.6787,7.4417,3.7238,7.4963,4.2032,14.849,7.7782,11.59,11.59,11.59,11.59]
# min_h=[8,15,19,22,31,34,39,40,41,43,45,46,46,51,53,49,57,60,54,58,68,54,54,54,54]
# max_h=[21,37,49,54,56,60,66,66,70,69,73,76,76,74,75,79,67,82,64,79,79,75,75,75,75]

# dists = truncated.(Normal.(μh, σh), min_h, max_h)

# all_ages = DataFrame(age = 1:25)
# p_age = @chain age_length begin
#     @by(:age, :p = length(:length))
#     DataFramesMeta.@transform(:p = :p / sum(:p))
#     rightjoin(all_ages, on=:age)
#     DataFramesMeta.@transform(:p = replace(:p, missing => 0))
#     select(:p)
#     Array()
#     vec()
# end

# mm = MixtureModel(dists[1:25], p_age)
# plot(mm)
