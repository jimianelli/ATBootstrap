"""ATBootstrap Python port."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .types import ATSurveyData, BootSpecs, ScalingClassProblem, ATBootstrapProblem
from .preprocess import get_survey_grid, preprocess_survey_data
from .spatial import (
    define_conditional_sim,
    get_lungs_params,
    parameterize_zdists,
    simulate_nasc,
    nonneg_lumult,
    trawl_assignments,
    choose_z_distribution,
    zdist_mean,
)
from .mace_ts import make_ts_function, to_linear
from .age_length import make_age_length_function
from .selectivity import make_selectivity_function, apply_selectivity
from .length_weight import make_weight_function
from .calibration import simulate_cal_error
from .scaling import resample_df, resample_scaling, get_trawl_category_means
from .nearbottom import make_nearbottom_dict, apply_nearbottom_coefficient, nearbottom_intercept


__all__ = [
    "preprocess_survey_data",
    "ATSurveyData",
    "BootSpecs",
    "ScalingClassProblem",
    "ATBootstrapProblem",
    "resample_df",
    "nonneg_lumult",
    "simulate_nasc",
    "trawl_assignments",
    "simulate",
    "simulate_class",
    "stepwise_error",
    "read_survey_files",
    "make_atbootstrap_problem",
    "make_scaling_class_problem",
]


error_labels = pd.DataFrame({
    "added_error": [
        "calibration", "simulate_nasc", "selectivity", "resample_scaling",
        "nearbottom_coefs", "trawl_assignments", "predict_ts", "age_length",
        "weights_at_age", "All",
    ],
    "error_label": [
        "Calibration", "Spatial sampling", "Selectivity", "Resample catches",
        "Nearbottom coefs", "Trawl assignment", "TS models", "Age-length",
        "Length-weight", "All",
    ],
})


def make_scaling_class_problem(surveydata: ATSurveyData, class_name: str,
                              cal_error: float = 0.1, age_max: int = 10,
                              maxlag: float = 200.0, nlags: int = 20,
                              weightfunc=lambda h: 1 / h,
                              zdist_candidates=None, aged_species=(21740,)):
    zdist_candidates = zdist_candidates or ["Gamma", "InverseGamma", "InverseGaussian", "LogNormal"]
    acoustics_sub = surveydata.acoustics[surveydata.acoustics["class"] == class_name]
    variogram, geosetup = define_conditional_sim(acoustics_sub, surveydata.grid,
                                                 maxlag=maxlag, nlags=nlags, weightfunc=weightfunc)
    params = get_lungs_params(geosetup, variogram["model"])
    optimal_dist = zdist_candidates[0]
    if len(zdist_candidates) > 1:
        optimal_dist = choose_z_distribution(zdist_candidates, acoustics_sub["nasc"].to_numpy(), params)
    zdists = parameterize_zdists(optimal_dist, params)
    return ScalingClassProblem(class_name, variogram, geosetup, params, optimal_dist,
                               zdists, cal_error, age_max, aged_species)


def make_atbootstrap_problem(surveydata: ATSurveyData, age_max: int = 10,
                             aged_species=(21740,), scaling_classes=None, cal_error: float = 0.1,
                             zdist_candidates=None, maxlag: float = 200.0,
                             nlags: int = 10, weightfunc=lambda h: 1 / h):
    scaling_classes = scaling_classes or surveydata.acoustics["class"].unique()
    class_problems = []
    for cls in scaling_classes:
        print(f"Preparing {cls}...")
        class_problems.append(
            make_scaling_class_problem(
                surveydata, cls, maxlag=maxlag, nlags=nlags,
                age_max=age_max, cal_error=cal_error,
                weightfunc=weightfunc, zdist_candidates=zdist_candidates,
                aged_species=aged_species,
            )
        )
    return ATBootstrapProblem(class_problems, scaling_classes, age_max, aged_species)


def simulate_class_iteration(scp: ScalingClassProblem, surveydata: ATSurveyData,
                             bs: BootSpecs, i: int = 1):
    scaling_sub = surveydata.scaling[surveydata.scaling["class"] == scp.class_name]
    surveygrid_coords = surveydata.grid[["x", "y"]].to_numpy()
    z0 = np.array([zdist_mean(scp.zfamily, p) for p in scp.zdists])

    selectivity_function = make_selectivity_function(bs.selectivity)
    scaling_boot = resample_scaling(scaling_sub, bs.resample_scaling)
    scaling_boot = apply_selectivity(scaling_boot, selectivity_function)

    predict_ts = make_ts_function()
    predict_age = make_age_length_function(surveydata.age_length, scp.age_max, bs.age_length)

    scaling_boot = scaling_boot.copy()
    scaling_boot["sigma_bs"] = scaling_boot.apply(
        lambda r: to_linear(predict_ts(r["ts_relationship"], r["ts_length"], bs.predict_ts)), axis=1
    )
    scaling_boot["age"] = scaling_boot["primary_length"].apply(predict_age)

    geotrawls = (
        scaling_boot.groupby(["haul_id", "class"], as_index=False)
        .first()[["haul_id", "class"]]
        .merge(surveydata.trawl_locations, on="haul_id", how="inner")
        .loc[:, ["haul_id", "x", "y"]]
    )

    ii = trawl_assignments(surveygrid_coords, geotrawls[["x", "y"]].to_numpy(),
                           bs.trawl_assignments)

    nasc = simulate_nasc(scp) if bs.simulate_nasc else nonneg_lumult(scp.params, z0)
    cal_error_sim = simulate_cal_error(scp.cal_error, bs.calibration)
    nasc_df = pd.DataFrame({
        "nasc": nasc * cal_error_sim,
        "haul_id": geotrawls["haul_id"].to_numpy()[ii],
    })

    if scp.class_name == "BT":
        from pathlib import Path
        species_path = Path(__file__).resolve().parents[2] / "surveydata" / "species.csv"
        nearbottom_dict = make_nearbottom_dict(
            pd.read_csv(species_path), bs.nearbottom_coefs
        )
        scaling_boot = apply_nearbottom_coefficient(scaling_boot, nearbottom_dict)
        scaling_boot["sigma_bs"] = scaling_boot.apply(
            lambda r: r["sigma_bs"] if r["species_code"] == 21740
            else to_linear(predict_ts(r["ts_relationship"], r["ts_length"], False)), axis=1
        )
        nasc_df["nasc"] = np.maximum(nasc_df["nasc"] - nearbottom_intercept, 0)

    predict_weight = make_weight_function(surveydata.length_weight)
    trawl_means_cat = get_trawl_category_means(scaling_boot, scp.aged_species, predict_weight)

    category_nasc = trawl_means_cat.pivot_table(index="haul_id", columns="category",
                                                values="p_nasc", fill_value=0).reset_index()
    category_sigma = trawl_means_cat.loc[:, ["haul_id", "category", "sigma_bs"]]
    category_weight = trawl_means_cat.loc[:, ["haul_id", "category", "weight"]]

    df = nasc_df.merge(category_nasc, on="haul_id", how="left")
    df = df.melt(id_vars=["haul_id", "nasc"], var_name="category", value_name="p_nasc")
    df = df.merge(category_sigma, on=["haul_id", "category"], how="left")
    df = df.merge(category_weight, on=["haul_id", "category"], how="left")

    df["n"] = df["nasc"] * df["p_nasc"] / (4 * np.pi * df["sigma_bs"]) * surveydata.dA
    df["biomass"] = df["n"] * df["weight"]

    df = (
        df.groupby("category", as_index=False)
        .agg(n=("n", "sum"), biomass=("biomass", "sum"))
    )
    df["species_code"] = df["category"].str.split("@").str[0].astype(int)
    df["age"] = df["category"].str.split("@").str[1].astype(int)
    df["i"] = i
    df["class"] = scp.class_name
    df = df.loc[:, ["class", "i", "species_code", "age", "category", "n", "biomass"]]
    df = df.sort_values(["i", "species_code", "age", "n", "biomass"])
    return df


def simulate_class(scp: ScalingClassProblem, surveydata: ATSurveyData,
                   nreplicates: int = 500, bs: BootSpecs | None = None):
    bs = bs or BootSpecs()
    print(f"Bootstrapping {scp.class_name}...")
    results = [simulate_class_iteration(scp, surveydata, bs, i + 1) for i in range(nreplicates)]
    class_results = pd.concat(results, ignore_index=True)
    class_results["class"] = scp.class_name
    return class_results


def simulate(atbp: ATBootstrapProblem, surveydata: ATSurveyData, nreplicates: int = 500,
             bs: BootSpecs | None = None, report_species=(21740,), report_ages=None):
    bs = bs or BootSpecs()
    class_results = [simulate_class(p, surveydata, nreplicates, bs) for p in atbp.class_problems]

    report_species = list(report_species) if report_species is not None else []
    if len(report_species) == 0:
        report_species = surveydata.scaling["species_code"].unique().tolist()
    if report_ages is None:
        report_ages = list(range(1, atbp.class_problems[0].age_max + 1))
    if len(report_ages) == 0:
        report_ages = [-1] + list(range(1, atbp.class_problems[0].age_max + 1))

    results = pd.concat(class_results, ignore_index=True)
    results = results[
        results["age"].isin(report_ages) & results["species_code"].isin(report_species)
    ]
    results = (
        results.groupby(["i", "species_code", "age"], as_index=False)
        .agg(n=("n", "sum"), biomass=("biomass", "sum"))
    )
    results["age"] = results["age"].replace({-1: np.nan})
    return results


def stepwise_error(atbp: ATBootstrapProblem, surveydata: ATSurveyData,
                   nreplicates: int = 500, remove: bool = False):
    error_sources = list(BootSpecs().__dict__.keys())
    colname = "eliminated_error" if remove else "added_error"
    results = []
    for i, err in enumerate(error_sources, start=1):
        prefix = "Omitting" if remove else "Adding"
        print(f"\n{prefix} {err} ({i}/{len(error_sources)})...")
        errs = {k: remove for k in error_sources}
        errs[err] = not remove
        bs = BootSpecs(**errs)
        res = simulate(atbp, surveydata, nreplicates=nreplicates, bs=bs)
        res[colname] = err
        results.append(res)
    return pd.concat(results, ignore_index=True)


def read_survey_files(surveydir: str):
    acoustics = pd.read_csv(f"{surveydir}/acoustics_projected.csv")
    trawl_locations = pd.read_csv(f"{surveydir}/trawl_locations_projected.csv")
    scaling = pd.read_csv(f"{surveydir}/scaling.csv")
    scaling["sample_correction_scalar"] = scaling["sample_correction_scalar"].astype(float)
    age_length = pd.read_csv(f"{surveydir}/age_length.csv")
    length_weight = pd.read_csv(f"{surveydir}/length_weight.csv")
    surveygrid = pd.read_csv(f"{surveydir}/surveygrid.csv")
    surveygrid = surveygrid.sample(frac=1, random_state=0).reset_index(drop=True)
    return {
        "acoustics": acoustics,
        "scaling": scaling,
        "age_length": age_length,
        "length_weight": length_weight,
        "trawl_locations": trawl_locations,
        "surveygrid": surveygrid,
    }
