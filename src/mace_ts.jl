
TSSpec(f, s) = (;f, s)

const TS_SE_DEFAULT = 3.0

function euphausiids_15_65mm_38khz(L)
    A = -9.30429983e2
    B = 3.21027896e0
    C = 1.74003785e0
    D = 1.36133896e-8
    E = -2.26958555e-6
    F = 1.50291244e-4
    G = -4.86306872e-3
    H = 7.38748423e-2
    I = -4.08004891e-1
    J = -7.39078690e1
    Lo = 3.835e-2  
    
    #  soundspeed in m/s
    c = 1470.0
    #  ***frequency in kHz***
    nu = 38.0

    #  wavelength in m
    lam = c / (nu * 10^3)
    #  wavenumber in m^-1
    k = 1 / lam
    #  wavenumber in radians per m
    k = 2 * pi * (nu * 10^3) / c

    if L < 1.5
        TS = -105.0
    elseif L > 6.5
        TS = -73.0
    else
        #  convert input length L in cm to m
        L = L / 100

        #  Demer and Conti, 2006, Vol 63: 928-935, Eq. 10
        TS = (A * (log10(B * k * L) / (B * k * L))^C +
              D * ((k * L)^6) +
              E * ((k * L)^5) +
              F * ((k * L)^4) +
              G * ((k * L)^3) +
              H * ((k * L)^2) +
              I * (k * L) +
              J + 20.0 * log10(L/Lo))
    end
    return TS
end

age0_pollock_ts(L) = 20log10(L) - 64.86
arctic_cod_ts(L) = 8.03log10(L) - 60.78
capelin_ts(L) = 20log10(L) - 70.3
chrysaora_melanaster_ts(L) = 10log10(π * (2L)^2) - 86.8
eulachon_ts(L) = 20log10(L) - 84.5
eulachon_new_ts(L) = 20log10(L) - 84.5
generalized_physoclist_ts(L) = 20log10(L) - 67.4
generic_fish_no_swimbladder_ts(L) = 20log10(L) - 83.2
generic_swimbladder_fish_ts(L) = generalized_physoclist_ts(L)
herring_75m_v2_ts(L) = 20log10(L) - log10(1+75/10) - 65.4
myctophids_sleucopsarus_ts(L) = 32.1log(log10(L)) - 64.1
pacific_hake_ts(L) = 20log10(L) - 68.0
sandlance_ts(L) = 56.5log10(L) - 125.1
squids_ts(L) = 20log10(L) - 75.4
standard_pollock_ts(L) = 20log10(L) - 66

const ts_lookup = Dict(                   # TS function (38 kHz, length in cm)    stdev
    "age0_pollock" =>               TSSpec(age0_pollock_ts,                 TS_SE_DEFAULT),
    "arctic_cod" =>                 TSSpec(arctic_cod_ts,                   TS_SE_DEFAULT),
    "capelin" =>                    TSSpec(capelin_ts,                      TS_SE_DEFAULT),
    "chrysaora_melanaster" =>       TSSpec(chrysaora_melanaster_ts,         TS_SE_DEFAULT),
    "eulachon" =>                   TSSpec(eulachon_ts,                     TS_SE_DEFAULT),
    "eulachon_new" =>               TSSpec(eulachon_new_ts,                 TS_SE_DEFAULT),
    "euphausiids_15_65mm_38khz" =>  TSSpec(euphausiids_15_65mm_38khz,       TS_SE_DEFAULT),
    "generalized_physoclist" =>     TSSpec(generalized_physoclist_ts,       TS_SE_DEFAULT),
    "generic_fish_no_swimbladder" =>TSSpec(generic_fish_no_swimbladder_ts,  TS_SE_DEFAULT),
    "generic_swimbladder_fish" =>   TSSpec(generic_swimbladder_fish_ts,     TS_SE_DEFAULT),
    "herring_75m_v2" =>             TSSpec(herring_75m_v2_ts,               TS_SE_DEFAULT),
    "myctophids_sleucopsarus" =>    TSSpec(myctophids_sleucopsarus_ts,      TS_SE_DEFAULT),
    "pacific_hake" =>               TSSpec(pacific_hake_ts,                 TS_SE_DEFAULT),
    "sandlance" =>                  TSSpec(sandlance_ts,                    TS_SE_DEFAULT),
    "squids" =>                     TSSpec(squids_ts,                       TS_SE_DEFAULT),
    "standard_pollock" =>           TSSpec(standard_pollock_ts,             0.14),
)

function make_ts_function()#stochastic=false)
    error_dict = Dict(relationship => randn() * ts_lookup[relationship].s
       for relationship in keys(ts_lookup))
    
    function predict_ts(relationship, L, stochastic=false)
        f, s = ts_lookup[relationship]
        err = stochastic ? error_dict[relationship] : 0.0
        f(L) + err
    end
    return predict_ts
end

to_linear(x) = exp10(x / 10)
