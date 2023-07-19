
TSSpec(f, s) = (;f, s)

ts_lookup = Dict(                          # TS function (38 kHz, length in cm)     stdev
    "age0_pollock" =>               TSSpec(L -> 20log10(L) - 64.86,                 1.0),
    "arctic_cod" =>                 TSSpec(L -> 8.03log10(L) - 60.78,               1.0),
    "capelin" =>                    TSSpec(L -> 20log10(L) - 70.3,                  1.0),
    "chrysaora_melanaster" =>       TSSpec(L -> 10log10(π * (2L)^2) - 86.8,         1.0),
    "eulachon" =>                   TSSpec(L -> 20log10(L) - 84.5,                  1.0),
    "eulachon_new" =>               TSSpec(L -> 20log10(L) - 84.5,                  1.0),
    # "euphausiids_15_65mm_38kHz"
    "generalized_physoclist" =>     TSSpec(L -> 20log10(L) - 67.4,                  1.0),
    "generic_fish_no_swimbladder" =>TSSpec(L -> 20log10(L) - 83.2,                  1.0),
    "generic_swimbladder_fish" =>   TSSpec(L -> 20log10(L) - 67.4,                  1.0),
    "herring_75m_v2" =>             TSSpec(L -> 20log10(L) - log10(1+75/10) - 65.4, 1.0),
    "myctophids_sleucopsarus" =>    TSSpec(L -> 32.1log(log10(L)) - 64.1,           1.0),
    "pacific_hake" =>               TSSpec(L -> 20log10(L) - 68.0,                  1.0),
    "sandlance" =>                  TSSpec(L -> 56.5log10(L) - 125.1,               1.0),
    "squids" =>                     TSSpec(L -> 20log10(L) - 75.4,                  1.0),
    "standard_pollock" =>           TSSpec(L -> 20log10(L) - 66,                    0.14),
)

function make_ts_function(stochastic=false)
    function predict_ts(relationship, L)
        f, s = ts_lookup[relationship]
        err = stochastic ? s * randn() : 0.0
        f(L) + err
    end
    return predict_ts
end
