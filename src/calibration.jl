


cals = CSV.read(joinpath(@__DIR__, "../surveydata/calibration_results.csv"), DataFrame)
# using StatsPlots
# @df cals plot(:sa_corr_calc, marker=:o, label="Calculated")
# @df cals plot!(:sa_corr_ek80, marker=:o, label="EK80 wizard")
# @df cals scatter(:sa_corr_calc, :sa_corr_ek80)

# std(cals.sa_corr_calc)
# std(skipmissing(cals.sa_corr_ek80))
# std(cals.sa_corr_calc) / sqrt(nrow(cals))
# std(skipmissing(cals.sa_corr_ek80)) / sqrt(sum(!ismissing, cals.sa_corr_ek80))

const CAL_ERROR = 0.05

function simulate_cal_error(cal_error, stochastic=true) 
    if stochastic
        return exp10(cal_error*randn() / 10)
    else
        return exp10(0.0)
    end
end