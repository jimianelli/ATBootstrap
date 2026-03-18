"""Selectivity functions."""
from __future__ import annotations

import numpy as np


def make_selectivity_function(stochastic: bool = True):
    mu_lfs = np.array([-2.4664253, 0.2698528])
    Sigma_lfs = np.array([[0.41987457, -0.0233237], [-0.0233237, 0.001459591]])
    mu_awt = np.array([-1.0558410, 0.1741619])
    Sigma_awt = np.array([[0.49771768, -0.01810365], [-0.01810365, 0.0008656043]])

    rng = np.random.default_rng()
    beta_lfs = rng.multivariate_normal(mu_lfs, Sigma_lfs) if stochastic else mu_lfs
    beta_awt = rng.multivariate_normal(mu_awt, Sigma_awt) if stochastic else mu_awt

    def selectivity(L, survey):
        if survey < 202001:
            return np.exp(beta_awt[0] + L * beta_awt[1]) / (1 + np.exp(beta_awt[0] + L * beta_awt[1]))
        return np.exp(beta_lfs[0] + L * beta_lfs[1]) / (1 + np.exp(beta_lfs[0] + L * beta_lfs[1]))

    return selectivity


def apply_selectivity(scaling, selectivity_function):
    scaling = scaling.copy()
    mask = (scaling["species_code"] == 21740) & (scaling["class"] != "BT")
    s = scaling.loc[mask].apply(lambda r: selectivity_function(r["primary_length"], r["survey"]), axis=1)
    scaling.loc[mask, "sample_correction_scalar"] = 1 / s
    scaling.loc[mask, "w"] = (
        scaling.loc[mask, "catch_sampling_expansion"]
        * scaling.loc[mask, "user_defined_expansion"]
        * scaling.loc[mask, "sample_correction_scalar"]
        * scaling.loc[mask, "haul_weight"]
    )
    return scaling
