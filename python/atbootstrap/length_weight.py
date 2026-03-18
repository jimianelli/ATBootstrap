"""Length-weight relationships."""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm


def make_weight_function(length_weight: pd.DataFrame, stochastic: bool = True, nmin: int = 5):
    n = len(length_weight)
    rng = np.random.default_rng()
    ii = rng.integers(0, n, n) if stochastic else np.arange(n)
    lw_boot = length_weight.iloc[ii].copy()

    X = np.log(lw_boot["fork_length"].values)
    y = np.log(lw_boot["organism_weight"].values)
    X = sm.add_constant(X)
    model = sm.OLS(y, X).fit()
    mu_loga, mu_b = model.params
    a = np.exp(mu_loga)
    b = mu_b

    Lmax = 100
    all_lengths = pd.DataFrame({"fork_length": np.arange(1, Lmax + 1)})
    binned = (
        lw_boot.merge(all_lengths, on="fork_length", how="right")
        .groupby("fork_length", as_index=False)
        .agg(mean_weight=("organism_weight", "mean"), n=("organism_weight", "size"))
    )
    binned["weight"] = np.where(
        binned["n"] >= nmin,
        binned["mean_weight"],
        a * binned["fork_length"] ** b,
    )
    weight_dict = dict(zip(binned["fork_length"], binned["weight"]))

    def weight_function(L):
        L1 = int(np.clip(np.round(L), 1, Lmax))
        return weight_dict[L1]

    return weight_function
