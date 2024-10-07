# glmm_fit_sel=exp(betas[1]+length_bins*betas[2])/(1+exp(betas[1]+length_bins*betas[2]))

function make_selectivity_function(stochastic=true)
    # LFS curve
    μ_lfs = [-2.4664253, 0.2698528]
    Σ_lfs = [0.41987457 -0.0233237;
            -0.0233237  0.001459591]
    D_lfs = MvNormal(μ_lfs, Σ_lfs)
    β_lfs = stochastic ? rand(D_lfs) : μ_lfs

    # AWT curve
    μ_awt = [-1.0558410, 0.1741619]
    Σ_awt = [0.49771768 -0.01810365
            -0.01810365  0.0008656043]
    D_awt = MvNormal(μ_awt, Σ_awt)
    β_awt = stochastic ? rand(D_awt) : μ_awt

    function selectivity(L, survey)
        if survey < 202001
            return exp(β_awt[1]+L*β_awt[2])/(1+exp(β_awt[1]+L*β_awt[2]))
        else
            return exp(β_lfs[1]+L*β_lfs[2])/(1+exp(β_lfs[1]+L*β_lfs[2]))
        end
    end

    return selectivity
end

function apply_selectivity!(scaling, selectivity_function)
    # w = Array{eltype(scaling.w)}(undef, nrow(scaling))
    # for (i, r) in enumerate(eachrow(scaling))
    #     L = r.primary_length
    #     if r.species_code == 21740 && r.class != "BT"
    #         s = selectivity_function(L, r.survey)
    #         scaling.sample_correction_scalar[i] = 1 / s
    #         w[i] = r.catch_sampling_expansion * r.user_defined_expansion *
    #             r.sample_correction_scalar * r.haul_weight
    #     end
    # end
    # scaling[:, :w] .= w

    @eachrow! scaling begin
        L = :primary_length
        if :species_code == 21740 && :class != "BT"
            s = selectivity_function(L, :survey)
            :sample_correction_scalar = 1 / s
            :w = :catch_sampling_expansion * :user_defined_expansion * :sample_correction_scalar * :haul_weight
        end
    end
end



