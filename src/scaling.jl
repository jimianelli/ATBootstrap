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

"""
Create a data frame of all the scatterer types into which backscatter should be allocated.
By default, this will be every species caught in the trawls, and further subdivided by 
year-class for pollock. 
"""
function make_all_categories(scaling, age_max, aged_species=[21740])
    event_ids = unique(scaling.event_id)
    ages = 0:age_max
    unaged_species = setdiff(event_ids, aged_species)
    with_ages = allcombinations(DataFrame,
        event_id = event_ids,
        age = ages,
        species_code = aged_species)
    without_ages = allcombinations(DataFrame,
        event_id = event_ids,
        age = [-1],
        species_code = unaged_species)
    all_categories = @chain [with_ages; without_ages] begin
        @orderby(:event_id, :species_code, :age)
        DataFramesMeta.@transform(
            :category = [string(s) * "@" * string(a) for (s, a) in zip(:species_code, :age)]
        )
    end
    return all_categories
end

function pollock_weights_at_age(scaling, length_weight, all_ages, stochastic=false)
    length_weight_pollock = @subset(length_weight, :species_code .== 21740)
    predict_weight = make_weight_function(length_weight_pollock, stochastic)
    res = @chain scaling begin
        @subset(:species_code .== 21740) # only calculate for pollock
        DataFramesMeta.@transform(
            :weight = predict_weight.(:primary_length))
        rightjoin(all_ages, on=[:event_id, :age])
        DataFramesMeta.@transform(:weight = replace(:weight, missing => 0.0))
        @by(:age, :weight = mean(:weight))
        DataFramesMeta.@transform(:species_code = 21740)
    end
    return res
end

function _category(use_ages, species_code, age)
    if use_ages(species_code)
        return string.(species_code) .* "@" .* string.(age)
    else
        return string.(species_code) .* "@-1"
    end
end

function proportion_at_category(scaling, all_categories, aged_species=[21740])
    use_ages = in(aged_species)
    scaling.category .= ""
    comp = @chain scaling begin
        DataFramesMeta.@transform(:category = _category.(use_ages, :species_code, :age))
        @by([:event_id, :category], :p_cat=sum(:w))
        DataFrames.groupby(:event_id)
        DataFramesMeta.@transform(:p_cat = :p_cat ./ sum(:p_cat))
        rightjoin(all_categories, on=[:event_id, :category])
        DataFramesMeta.@transform(:p_cat = replace(:p_cat, missing => 0.0))
        select(:event_id, :category, :p_cat)
        @orderby(:event_id, :category)
    end
    return DataFrames.unstack(comp, :category, :p_cat, fill=0.0)
end
