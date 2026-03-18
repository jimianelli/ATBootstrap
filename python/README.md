# ATBootstrap Python Port

This folder will hold the Python translation of the Julia workflow in `src/`.

## Scope

- Convert core Julia modules in `src/` to Python modules with equivalent inputs/outputs.
- Maintain reproducibility of survey preprocessing, bootstrap sampling, spatial steps,
  and summary outputs.
- Preserve intermediate artifacts and tests to validate parity with Julia outputs.

## Proposed module mapping (Julia → Python)

- `src/ATBootstrap.jl` → `python/atbootstrap/__init__.py`
- `src/types.jl` → `python/atbootstrap/types.py`
- `src/preprocess_survey_data.jl` → `python/atbootstrap/preprocess.py`
- `src/calibration.jl` → `python/atbootstrap/calibration.py`
- `src/nearbottom.jl` → `python/atbootstrap/nearbottom.py`
- `src/mace_ts.jl` → `python/atbootstrap/mace_ts.py`
- `src/mace_selectivity.jl` → `python/atbootstrap/selectivity.py`
- `src/mace_age_length.jl` → `python/atbootstrap/age_length.py`
- `src/mace_length_weight.jl` → `python/atbootstrap/length_weight.py`
- `src/spatial.jl` → `python/atbootstrap/spatial.py`
- `src/scaling.jl` → `python/atbootstrap/scaling.py`
- `src/display.jl` → `python/atbootstrap/display.py`

## Dependencies (expected)

- numpy, pandas
- scipy
- matplotlib
- geopandas/shapely (if spatial ops require geometry)

## Status

Initial Python port implemented. See `docs/python-port-plan.qmd` for the porting plan
and parity checks.

## Install (editable)

From the repo root:

```
cd python
python -m pip install -e .
```

## Basic usage

```python
from atbootstrap import read_survey_files, ATSurveyData, make_atbootstrap_problem, simulate

files = read_survey_files("surveydata/202207")

surveydata = ATSurveyData(
    acoustics=files["acoustics"],
    scaling=files["scaling"],
    age_length=files["age_length"],
    length_weight=files["length_weight"],
    trawl_locations=files["trawl_locations"],
    grid=files["surveygrid"],
    dA=1.0,
)

atbp = make_atbootstrap_problem(surveydata)
res = simulate(atbp, surveydata, nreplicates=100)
```
