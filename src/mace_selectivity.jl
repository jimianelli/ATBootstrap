
#=
Questions for Kresimir:
    - Which functions to apply in EBS?
        - 
    - Which column in SKSD are these calculating?
        - sample_correction_scalar (reciprocal of p)
    - What are confidence intervals for parameters of each?
    - Citation for these?
        - 201906 report?
        - 202107 report...(not out yet)?
    - How is age calculated?
    - How is weight calculated?
    - How to re-run analysis for nearest-haul myself?
        -
=#

# (Intercept)      Length 
#  -4.0207641   0.3533825 

# 2 x 2 Matrix of class "dpoMatrix"
#             (Intercept)       Length
# (Intercept)  1.04110776 -0.085418314
# Length      -0.08541831  0.007097557


# > L50_95_range
#     2.5%    97.5% 
# 10.35141 11.98987 
# > SR_95_range
#      2.5%     97.5% 
#  4.207268 11.592823 


# glmm_fit_sel=exp(betas[1]+length_bins*betas[2])/(1+exp(betas[1]+length_bins*betas[2]))


β = [-4.0207641, 0.3533825]
f(L) = exp(β[1]+L*β[2])/(1+exp(β[1]+L*β[2]))
plot(f, 0, 75)


function make_selectivity_function(species, age, L, )
end