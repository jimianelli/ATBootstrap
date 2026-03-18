"""Scaling utilities."""
from __future__ import annotations

import numpy as np
import pandas as pd


def resample_df(df: pd.DataFrame, stochastic: bool = True) -> pd.DataFrame:
    n = len(df)
    if stochastic:
        ii = np.random.randint(0, n, n)
    else:
        ii = np.arange(n)
    return df.iloc[ii].copy()


def resample_scaling(df: pd.DataFrame, stochastic: bool = True) -> pd.DataFrame:
    return (
        df.groupby(["haul_id", "class"], as_index=False, group_keys=False)
        .apply(lambda x: resample_df(x, stochastic))
        .reset_index(drop=True)
    )


def _category(use_ages, species_code, age):
    if use_ages(species_code):
        return f"{species_code}@{age}"
    return f"{species_code}@-1"


def get_trawl_category_means(scaling: pd.DataFrame, aged_species, predict_weight):
    use_ages = lambda s: s in aged_species
    df = scaling.copy()
    df["category"] = df.apply(lambda r: _category(use_ages, r["species_code"], r["age"]), axis=1)
    df["weight"] = df["primary_length"].apply(predict_weight)
    df["p_nasc"] = df["sigma_bs"] * df["w"]
    df = df[df["w"] > 0]

    trawl_means_cat = (
        df.groupby(["haul_id", "category"], as_index=False)
        .agg(sigma_bs=("sigma_bs", "mean"),
             p_nasc=("p_nasc", "sum"),
             weight=("weight", "mean"))
    )
    trawl_means_cat["p_nasc"] = trawl_means_cat.groupby("haul_id")["p_nasc"].transform(
        lambda x: x / x.sum()
    )
    return trawl_means_cat
