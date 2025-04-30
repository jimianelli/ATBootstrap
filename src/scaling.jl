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
        DataFramesMeta.groupby(df, [:haul_id, :class]))
end

function get_trawl_category_means(scaling, aged_species, predict_weight)
    use_ages = in(aged_species)
    trawl_means_cat = @chain scaling begin
        DataFramesMeta.@transform(
            :category = _category.(use_ages, :species_code, :age),
            :weight = predict_weight.(:primary_length),
            :p_nasc = :sigma_bs .* :w
        )
        #=
        If a catch filter is applied, all specimens may be listed twice in the same haul, 
        with "user_defined_expansion" set to 0.0 for one listing and 1.0 for the other.
        This is done to split a haul into two parts, each applied to a different scaling
        stratum. Eliminate weighting factors == 0 here, since they don't affect any of the
        numbers and introduce NaNs if all specimens in a haul/category have :w .== 0.
        =# 
        @subset(:w .> 0)
        @by([:haul_id, :category], 
            :sigma_bs = mean(:sigma_bs),#, Weights(:w)),
            :p_nasc = sum(:p_nasc),
            :weight = mean(:weight)#, Weights(:w)),
        )
        # make proportions sum to 1.0 for each haul
        DataFrames.groupby(:haul_id)
        DataFramesMeta.@transform(:p_nasc = :p_nasc ./ sum(:p_nasc))
    end
    return trawl_means_cat
end

function _category(use_ages, species_code, age)
    if use_ages(species_code)
        return string.(species_code) .* "@" .* string.(age)
    else
        return string.(species_code) .* "@-1"
    end
end
