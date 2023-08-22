


cals = CSV.read(joinpath(@__DIR__, "../surveydata/calibration_results.csv"), DataFrame)
using StatsPlots
@df cals plot(:sa_corr_calc, marker=:o, label="Calculated")
@df cals plot!(:sa_corr_ek80, marker=:o, label="EK80 wizard")
@df cals scatter(:sa_corr_calc, :sa_corr_ek80)

@df @subset(cals, :sphere_id.=="WC38.1_33") histogram(:sa_corr_calc)
@df @subset(cals, :sphere_id.=="WC38.1_52") vline!(:sa_corr_calc)

std(cals.sa_corr_calc)
std(skipmissing(cals.sa_corr_ek80))
std(cals.sa_corr_calc) / sqrt(nrow(cals))
std(skipmissing(cals.sa_corr_ek80)) / sqrt(sum(!ismissing, cals.sa_corr_ek80))
std(cals.sa_corr_calc) / sqrt(3)
std(skipmissing(cals.sa_corr_ek80)) / sqrt(3)
10log10(std(exp10.(cals.sa_corr_calc / 10)) / sqrt(3) + 1)
10log10(std(skipmissing(exp10.(cals.sa_corr_ek80 / 10))) / sqrt(3) + 1)

const CAL_ERROR = 0.02

function simulate_cal_error(cal_error, stochastic=true) 
    if stochastic
        return exp10(cal_error*randn() / 10)
    else
        return exp10(0.0)
    end
end