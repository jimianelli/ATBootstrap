"""Simple parity evaluation: Python vs Julia results for survey 202207."""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore", message="You are merging on int and float columns")

from atbootstrap import (
    preprocess_survey_data,
    read_survey_files,
    ATSurveyData,
    make_atbootstrap_problem,
    simulate,
)
from atbootstrap.preprocess import TransectRibbons


def main():
    repo_root = Path(__file__).resolve().parents[2]
    surveydir = repo_root / "surveydata" / "202207"
    results_dir = repo_root / "analyses" / "results"
    out_dir = repo_root / "docs" / "artifacts"
    out_dir.mkdir(parents=True, exist_ok=True)

    survey = "202207"
    resolution = 10.0
    km2nmi = 1 / 1.852
    dA = (resolution * km2nmi) ** 2
    log_ranges = [(250, 2190), (2752.5, 2778), (2400, 2752.49), (2780, 6100)]

    preprocess_survey_data(
        str(surveydir),
        dx=resolution,
        ebs=True,
        log_ranges=log_ranges,
        grid_method=TransectRibbons(transect_width=40.0),
    )

    files = read_survey_files(str(surveydir))
    acoustics = files["acoustics"]
    scaling = files["scaling"]

    scaling_classes = ["SS1", "SS2", "BT"]
    acoustics = acoustics[acoustics["class"].isin(scaling_classes)]
    acoustics = acoustics[acoustics["transect"] % 2 == 0]

    surveydata = ATSurveyData(
        acoustics=acoustics,
        scaling=scaling,
        age_length=files["age_length"],
        length_weight=files["length_weight"],
        trawl_locations=files["trawl_locations"],
        grid=files["surveygrid"],
        dA=dA,
    )

    atbp = make_atbootstrap_problem(
        surveydata, scaling_classes=scaling_classes
    )

    nreplicates = 100
    results_py = simulate(atbp, surveydata, nreplicates=nreplicates)

    # cache python results for comparison
    py_cache = results_dir / f"results_python_{survey}.csv"
    results_py.to_csv(py_cache, index=False)

    # Julia results baseline
    results_jl = pd.read_csv(results_dir / f"results_{survey}.csv")

    def summarize(df):
        df = df[df["species_code"] == 21740]
        return (
            df.groupby("age", as_index=False)
            .agg(n_mean=("n", "mean"), biomass_mean=("biomass", "mean"))
            .sort_values("age")
        )

    py_sum = summarize(results_py)
    jl_sum = summarize(results_jl)

    merged = py_sum.merge(jl_sum, on="age", suffixes=("_py", "_jl"))
    merged["n_rel_diff"] = (merged["n_mean_py"] - merged["n_mean_jl"]) / merged["n_mean_jl"]
    merged["biomass_rel_diff"] = (
        (merged["biomass_mean_py"] - merged["biomass_mean_jl"]) / merged["biomass_mean_jl"]
    )

    merged.to_csv(out_dir / f"python_parity_{survey}.csv", index=False)

    # Simple HTML report
    n_mape = float(np.mean(np.abs(merged["n_rel_diff"])) * 100)
    b_mape = float(np.mean(np.abs(merged["biomass_rel_diff"])) * 100)

    html = f"""<!doctype html>
<html>
<head><meta charset="utf-8"><title>Python parity evaluation {survey}</title></head>
<body>
<h1>Python parity evaluation ({survey})</h1>
<p>Comparison of Python vs Julia bootstrap mean-by-age for species 21740.</p>
<ul>
  <li>Python replicates: {nreplicates}</li>
  <li>Python cache: analyses/results/results_python_{survey}.csv</li>
  <li>Julia baseline: analyses/results/results_{survey}.csv</li>
  <li>Mean absolute percent diff (n): {n_mape:.2f}%</li>
  <li>Mean absolute percent diff (biomass): {b_mape:.2f}%</li>
</ul>
{merged.to_html(index=False)}
</body>
</html>
"""

    (repo_root / "docs" / "python-parity.html").write_text(html)

    qmd = f"""---
title: \"Python parity evaluation {survey}\"
format:
  html:
    toc: false
---

## Summary

- Python replicates: {nreplicates}
- Python cache: `analyses/results/results_python_{survey}.csv`
- Julia baseline: `analyses/results/results_{survey}.csv`
- Mean absolute percent diff (n): {n_mape:.2f}%
- Mean absolute percent diff (biomass): {b_mape:.2f}%

## Mean-by-age comparison (species 21740)

```{{r}}
library(tidyverse)
library(gt)

merged <- read_csv("docs/artifacts/python_parity_{survey}.csv", show_col_types = FALSE)
merged %>%
  gt()
```
"""

    (repo_root / "docs" / "python-parity.qmd").write_text(qmd)

    meta = {
        "survey": survey,
        "nreplicates": nreplicates,
        "n_mape": n_mape,
        "biomass_mape": b_mape,
        "note": "Full z-dist selection enabled.",
    }
    (out_dir / f"python_parity_{survey}.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
