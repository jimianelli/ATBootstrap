using CSV
using DataFrames
using DataFramesMeta
using Turing
using StatsPlots, StatsPlots.PlotMeasures

calfile = joinpath(@__DIR__, "../../surveydata/calibration_results.csv")
cals = CSV.read(calfile, DataFrame)
cals = @chain cals begin
    @subset(
        .!ismissing.(:sphere_id),
        .!ismissing.(:pulse_length),
        .!ismissing.(:sa_corr_calc)
    )
    @orderby(:sounder, :pulse_length, :sphere_id)
    @transform(:sounder_pulse = :sounder .* " " .* string.(:pulse_length) .* " ms")
end

cal_summary = @chain cals begin
    @by([:sounder, :pulse_length, :sphere_id], 
        :n = length(:sa_corr_calc),
        :sd = std(:sa_corr_calc))
    @subset(:n .> 1)
    @transform(:se = :sd ./ sqrt.(:n))
end
sa_std = mean(cal_summary.sd[isfinite.(cal_summary.sd) .& (.! ismissing.(cal_summary.sd))])
exp10(sa_std / 10)

date_ids = indexin(cals.date, unique(cals.date))
sphere_ids = indexin(cals.sphere_id, unique(cals.sphere_id))
sounder_pulse_ids = indexin(cals.sounder_pulse, unique(cals.sounder_pulse))

sounder_pulse_labels = unique(cals.sounder_pulse)
sphere_labels = unique(cals.sphere_id)

@model function calibration_model(sa_corr_calc, sa_corr_ek80, 
        sphere_label, sounder_pulse_id)
    n = length(sa_corr_calc)
    sphere_id = indexin(cals.sphere_id, unique(cals.sphere_id))
    nsounder_pulse = length(unique(sounder_pulse_id))

    σG ~ Exponential(0.5)
    G ~ filldist(Normal(0, σG), nsounder_pulse)

    σTS_Cu ~ Exponential(0.02)
    σTS_WC ~ Exponential(0.1)
    σTS = [contains.(x, "Cu") ? σTS_Cu : σTS_WC for x in sphere_label]
    ΔTS ~ arraydist(Normal.(0, σTS))
    σmeasurement ~ Exponential(0.05)

    for t in 1:n
        i = sounder_pulse_id[t]
        j = sphere_id[t]
        μ = G[i] + ΔTS[j]
        if !ismissing(sa_corr_calc[t])
            sa_corr_calc[t] ~ Normal(μ, σmeasurement)
        end
        if !ismissing(sa_corr_ek80[t])
            sa_corr_ek80[t] ~ Normal(μ, σmeasurement)
        end
    end
end


model = calibration_model(cals.sa_corr_calc, cals.sa_corr_ek80,
    sphere_labels, sounder_pulse_ids)
chns = sample(model, NUTS(), 5000)
plot(group(chns, :G))
plot(group(chns, :ΔTS))
plot(chns[:σmeasurement])

ekpal = [range(color("red"), color("blue"), length=4);
    range(color("goldenrod"), color("green"), length=3)];
density(Array(group(chns, :G)), fill=true, alpha=0.5, 
    label=reshape(sounder_pulse_labels, 1, :), palette=ekpal, legend=:top,
    size=(800, 500), margin=10px,
    xlabel="Sₐ correction (dB)", ylabel="Probability density")

summarize(group(chns, :G))

density(Array(group(chns, :ΔTS)), fill=true, alpha=0.5,
    label=reshape(sphere_labels, 1, :),
    palette=color.(["orange4", "slategrey", "royalblue3", "mediumpurple"]),
    xlabel="Sphere deviation (dB)")
vline!([0], color=:black, label="")

corner(group(chns, :ΔTS))