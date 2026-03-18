"""Calibration error model."""
from __future__ import annotations

import numpy as np
import pandas as pd


def load_calibration(surveydir: str):
    return pd.read_csv(f"{surveydir}/calibration_results.csv")


def simulate_cal_error(cal_error: float, stochastic: bool = True) -> float:
    if stochastic:
        return 10 ** (cal_error * np.random.randn() / 10)
    return 10 ** 0.0
