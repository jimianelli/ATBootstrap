"""Core dataclasses for ATBootstrap Python port."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd


@dataclass
class ATSurveyData:
    acoustics: pd.DataFrame
    scaling: pd.DataFrame
    age_length: pd.DataFrame
    length_weight: pd.DataFrame
    trawl_locations: pd.DataFrame
    grid: pd.DataFrame  # columns: x, y
    dA: float


@dataclass
class BootSpecs:
    selectivity: bool = True
    predict_ts: bool = True
    resample_scaling: bool = True
    nearbottom_coefs: bool = True
    age_length: bool = True
    weights_at_age: bool = True
    trawl_assignments: bool = True
    simulate_nasc: bool = True
    calibration: bool = True

    @classmethod
    def all(cls, value: bool) -> "BootSpecs":
        return cls(**{k: value for k in cls().__dict__.keys()})


@dataclass
class ScalingClassProblem:
    class_name: str
    variogram: dict
    geosetup: dict
    params: dict
    zfamily: str
    zdists: list
    cal_error: float
    age_max: int
    aged_species: Sequence[int]


@dataclass
class ATBootstrapProblem:
    class_problems: Sequence[ScalingClassProblem]
    scaling_classes: Sequence[str]
    age_max: int
    aged_species: Sequence[int]
