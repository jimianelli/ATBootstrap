"""Nearbottom coefficients."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import truncnorm


nearbottom_intercept = 3.43


def load_species(surveydir: str):
    return pd.read_csv(f"{surveydir}/species.csv")


def make_nearbottom_dict(species: pd.DataFrame, stochastic: bool = True):
    species = species.copy()
    species["nearbottom_group"] = "Misc"
    species.loc[species["species_code"] == 21740, "nearbottom_group"] = "Pollock"
    species.loc[species["species_code"] == 21725, "nearbottom_group"] = "Arctic cod"
    species.loc[(species["species_code"] >= 10110) & (species["species_code"] <= 10120),
                "nearbottom_group"] = "Large flatfish"
    species.loc[(species["species_code"] >= 30050) & (species["species_code"] <= 30535),
                "nearbottom_group"] = "Rockfish"

    nearbottom_coefs = pd.DataFrame({
        "nearbottom_group": ["Pollock", "Arctic cod", "Large flatfish", "Rockfish", "Misc"],
        "A": [2.52, 16.39, 0.85, 93.59, 11.63],
        "lower": [2.21, 2.84, 0.16, 8.63, 3.92],
        "upper": [2.86, 62.26, 1.67, 343.43, 22.09],
        "b": [66, 67.4, 67.4, 67.4, 67.4],
    })
    nearbottom_coefs["sd_A"] = ((nearbottom_coefs["A"] - nearbottom_coefs["lower"]) +
                                (nearbottom_coefs["upper"] - nearbottom_coefs["A"])) / 2 / 2

    nearbottom_coefs = species.merge(nearbottom_coefs, on="nearbottom_group", how="left")

    if stochastic:
        a_vals = []
        for _, r in nearbottom_coefs.iterrows():
            a, lower, upper, sd = r["A"], r["lower"], r["upper"], r["sd_A"]
            if pd.isna(a):
                a_vals.append(np.nan)
                continue
            a_std = (a - lower) / sd
            b_std = (upper - a) / sd
            a_vals.append(truncnorm.rvs(-a_std, b_std, loc=a, scale=sd))
        aa = np.array(a_vals)
    else:
        aa = nearbottom_coefs["A"].to_numpy()

    aa = aa / (10 ** (-nearbottom_coefs["b"] / 10)) / 1852 ** 2 / (4 * np.pi)
    return dict(zip(nearbottom_coefs["species_code"], aa))


def apply_nearbottom_coefficient(scaling: pd.DataFrame, nearbottom_dict: dict) -> pd.DataFrame:
    scaling = scaling.copy()

    def get_coef(code):
        if code in nearbottom_dict:
            return nearbottom_dict[code]
        return nearbottom_dict.get(24166, np.nan)

    scaling["user_defined_expansion"] = scaling["species_code"].apply(get_coef)
    scaling["w"] = (
        scaling["catch_sampling_expansion"]
        * scaling["user_defined_expansion"]
        * scaling["sample_correction_scalar"]
        * scaling["haul_weight"]
    )
    return scaling
