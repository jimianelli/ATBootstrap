
TSSpec(f, s) = (;f, s)



function krill_ts(L)
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
        TS = (A * (math.log10(B * k * L) / (B * k * L))^C +
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


ts_lookup = Dict(                          # TS function (38 kHz, length in cm)     stdev
    "age0_pollock" =>               TSSpec(L -> 20log10(L) - 64.86,                 1.0),
    "arctic_cod" =>                 TSSpec(L -> 8.03log10(L) - 60.78,               1.0),
    "capelin" =>                    TSSpec(L -> 20log10(L) - 70.3,                  1.0),
    "chrysaora_melanaster" =>       TSSpec(L -> 10log10(π * (2L)^2) - 86.8,         1.0),
    "eulachon" =>                   TSSpec(L -> 20log10(L) - 84.5,                  1.0),
    "eulachon_new" =>               TSSpec(L -> 20log10(L) - 84.5,                  1.0),
    "euphausiids_15_65mm_38khz" =>  TSSpec(krill_ts,                                1.0),
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
