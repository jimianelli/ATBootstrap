using UnderwaterAcoustics
using MonteCarloMeasurements
using StatsPlots

# Function to convert output of UnderwaterAcoustics.absorption to dB/m
# See https://github.com/org-arl/UnderwaterAcoustics.jl/issues/54
function alpha(frequency, salinity=34.0, temperature=27.0, depth=10.0, pH=8.1)
    return -20log10(absorption(frequency, 1, salinity, temperature, depth, pH)) 
end

# make it work with Particles type from MonteCarloMeasurements
register_primitive_multi(alpha) 

# Following values are from Haris et al. 2018, ICES JMS 75
# ES38-DD column of Table 3
G0_sd = [0.29, 0.06, 0.04, 0.17, 0.10, 0.18, 0.17, 0.14, 0.12, 0.22, 0.1]
G0_n = [9, 10, 7, 9, 6, 11, 9, 9, 9, 4, 4]
G0_se = G0_sd ./ sqrt.(G0_n)
G0_error = mean(G0_sd)

# ES38-DD column of table 4
ψ_sd = [0.17, 0.18, 0.28, 0.13, 0.25, 0.13, 0.2, 0.22]
ψ_n = [6, 15, 7, 4, 5, 5, 5, 10]
ψ_se = ψ_sd ./ sqrt.(ψ_n)
# ψ_error = mean(ψ_sd)
ψ_bias = 10log10.(1 .+ 0.01 * [-28, -25, -25, -25, -21, -15, -19, -21, -17, -16, 
    -19, -12, 08, -8, -5, -8, -5, -5, 
    -6, -8, -3, -5, -1, -5, -8, -10, -5])
ψ_error = abs(mean(ψ_bias))


T_true = [3.183706, 4.978313, 5.260582, 3.837085]
S_true = [32.187052, 32.437675, 32.339063, 32.363523]
T_assumed = 4.0
S_assumed = 33.0
ΔT = T_assumed .- T_true
ΔS = S_assumed .- S_true
mean(ΔT), std(ΔT)
mean(ΔS), std(ΔS)
r_true = 100                                # True range to target/volume (m)
c_true = soundspeed(mean(T_true), mean(S_true), r_true) # True sound speed (m/s)
T = T_assumed±std(ΔT)
S = S_assumed±std(ΔS)

f = 38e3                        # Frequency (Hz)
c = soundspeed(T, S, r_true)    # Assumed sound speed (m/s)
t = 2 * r_true / c_true         # Round-trip time for echo (s)
r = t * c / 2                   # Inferred/observed range (m)
α = alpha(f, S, T, r/2)         # Absorption coefficient (dB/m)
λ = c / f                       # Wavelength (m)
τ = 1.024e-3                    # Pulse length (s)
pt = 1000                       # Transmitted power (W)
ψ = exp10((-21±ψ_error) / 10)   # EBA (sr)
# ψ = exp10((-21) / 10)   # EBA (sr)
G0 = 0±G0_error                 # On-axis gain (logarithmic, dB)
g0 = exp10(G0 / 10)             # On-axis gain (linear, unitless)
g_linearity = 1.01±0.02         # Factor to correct for linearity
Δsphere = 0±0.05                # Difference of sphere TS from theoretical

# Sv = Pr + 20log10(r) + 2α*r - 10log10((pt * λ^2 * g0^2 * c * τ * ψ) / 32(π^2))
gains = 20log10(r) + 2α*r - 10log10((pt * λ^2 * g_linearity * g0^2 * c * τ * ψ) / 32(π^2)) + Δsphere
uncertainty = gains - mean(gains.particles)

cal_sd = pstd(uncertainty)
plot(uncertainty, normalize=true, label="", linewidth=0, foreground_color_legend=nothing,
    xlabel="Calibration error (dB)", ylabel="Probability density")
vline!([-cal_sd, cal_sd], label="S.D. = $(round(cal_sd, sigdigits=2)) dB", color=:black)
