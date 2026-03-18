"""Target-strength relationships."""
from __future__ import annotations

import numpy as np


def TSSpec(f, s):
    return {"f": f, "s": s}


TS_SE_DEFAULT = 3.0


def euphausiids_15_65mm_38khz(L):
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

    c = 1470.0
    nu = 38.0
    lam = c / (nu * 10 ** 3)
    k = 2 * np.pi * (nu * 10 ** 3) / c

    if L < 1.5:
        TS = -105.0
    elif L > 6.5:
        TS = -73.0
    else:
        L = L / 100
        TS = (
            A * (np.log10(B * k * L) / (B * k * L)) ** C
            + D * ((k * L) ** 6)
            + E * ((k * L) ** 5)
            + F * ((k * L) ** 4)
            + G * ((k * L) ** 3)
            + H * ((k * L) ** 2)
            + I * (k * L)
            + J
            + 20.0 * np.log10(L / Lo)
        )
    return TS


def age0_pollock_ts(L):
    return 20 * np.log10(L) - 64.86


def arctic_cod_ts(L):
    return 8.03 * np.log10(L) - 60.78


def capelin_ts(L):
    return 20 * np.log10(L) - 70.3


def chrysaora_melanaster_ts(L):
    return 10 * np.log10(np.pi * (2 * L) ** 2) - 86.8


def eulachon_ts(L):
    return 20 * np.log10(L) - 84.5


def eulachon_new_ts(L):
    return 20 * np.log10(L) - 84.5


def generalized_physoclist_ts(L):
    return 20 * np.log10(L) - 67.4


def generic_fish_no_swimbladder_ts(L):
    return 20 * np.log10(L) - 83.2


def generic_swimbladder_fish_ts(L):
    return generalized_physoclist_ts(L)


def herring_75m_v2_ts(L):
    return 20 * np.log10(L) - np.log10(1 + 75 / 10) - 65.4


def myctophids_sleucopsarus_ts(L):
    return 32.1 * np.log(np.log10(L)) - 64.1


def pacific_hake_ts(L):
    return 20 * np.log10(L) - 68.0


def sandlance_ts(L):
    return 56.5 * np.log10(L) - 125.1


def squids_ts(L):
    return 20 * np.log10(L) - 75.4


def standard_pollock_ts(L):
    return 20 * np.log10(L) - 66


ts_lookup = {
    "age0_pollock": TSSpec(age0_pollock_ts, TS_SE_DEFAULT),
    "arctic_cod": TSSpec(arctic_cod_ts, TS_SE_DEFAULT),
    "capelin": TSSpec(capelin_ts, TS_SE_DEFAULT),
    "chrysaora_melanaster": TSSpec(chrysaora_melanaster_ts, TS_SE_DEFAULT),
    "eulachon": TSSpec(eulachon_ts, TS_SE_DEFAULT),
    "eulachon_new": TSSpec(eulachon_new_ts, TS_SE_DEFAULT),
    "euphausiids_15_65mm_38khz": TSSpec(euphausiids_15_65mm_38khz, TS_SE_DEFAULT),
    "generalized_physoclist": TSSpec(generalized_physoclist_ts, TS_SE_DEFAULT),
    "generic_fish_no_swimbladder": TSSpec(generic_fish_no_swimbladder_ts, TS_SE_DEFAULT),
    "generic_swimbladder_fish": TSSpec(generic_swimbladder_fish_ts, TS_SE_DEFAULT),
    "herring_75m_v2": TSSpec(herring_75m_v2_ts, TS_SE_DEFAULT),
    "myctophids_sleucopsarus": TSSpec(myctophids_sleucopsarus_ts, TS_SE_DEFAULT),
    "pacific_hake": TSSpec(pacific_hake_ts, TS_SE_DEFAULT),
    "sandlance": TSSpec(sandlance_ts, TS_SE_DEFAULT),
    "squids": TSSpec(squids_ts, TS_SE_DEFAULT),
    "standard_pollock": TSSpec(standard_pollock_ts, 0.14),
}


def make_ts_function():
    rng = np.random.default_rng()
    error_dict = {rel: rng.normal(0, ts_lookup[rel]["s"]) for rel in ts_lookup.keys()}

    def predict_ts(relationship, L, stochastic=False):
        f = ts_lookup[relationship]["f"]
        err = error_dict[relationship] if stochastic else 0.0
        return f(L) + err

    return predict_ts


def to_linear(x):
    return 10 ** (x / 10)
