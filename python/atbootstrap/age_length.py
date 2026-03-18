"""Age-length relationships."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mode


def nearest_length(L: float, unique_lengths: np.ndarray) -> int:
    idx = np.argmin(np.abs(unique_lengths - L))
    return idx


def make_age_length_function(age_length: pd.DataFrame, age_max: int = 10,
                             stochastic: bool = True):
    key = age_length.copy()
    key["fork_length"] = key["fork_length"].round()
    key["age"] = np.where(key["age"] > age_max, age_max, key["age"])

    key = (
        key.groupby(["fork_length", "age"]).size().reset_index(name="n")
        .pivot_table(index="fork_length", columns="age", values="n", fill_value=0)
        .reset_index()
    )

    arr = key.iloc[:, 1:].to_numpy()
    arr = arr / arr.sum(axis=1, keepdims=True)
    fl = key["fork_length"].to_numpy()

    def sample_age(i):
        probs = arr[i]
        return np.random.choice(np.arange(1, arr.shape[1] + 1), p=probs)

    def predict_age(L):
        i = nearest_length(L, fl)
        if stochastic:
            return sample_age(i)
        return np.argmax(arr[i]) + 1

    return predict_age
