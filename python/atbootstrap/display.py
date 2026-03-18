"""Plotting utilities (matplotlib-based)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .spatial import simulate_nasc
from .scaling import resample_df


def plot_class_variograms(atbp, size=(8, 6), **kwargs):
    n = len(atbp.class_problems)
    fig, axes = plt.subplots(1, n, figsize=(size[0] * n, size[1]))
    if n == 1:
        axes = [axes]
    for ax, cp in zip(axes, atbp.class_problems):
        vg_emp = cp.variogram["empirical"]
        vg_mod = cp.variogram["model"]
        x = vg_emp["lags"]
        y = vg_emp["gamma"]
        ax.plot(x, y, marker="o", label="Empirical")
        h = np.linspace(0, np.nanmax(x), 200)
        nugget, sill, rng = vg_mod["nugget"], vg_mod["sill"], vg_mod["range"]
        ax.plot(h, nugget + sill * (1 - np.exp(-h / rng)), label="Model")
        ax.set_title(cp.class_name)
        ax.set_xlabel("Lag (km)")
        ax.set_ylabel("γ")
        ax.legend()
    return fig


def plot_simulated_nasc(scp, surveydata, simdomain=None, alpha=0.3,
                        markersize=2.2, max_bubblesize=15, **kwargs):
    simdomain = simdomain or surveydata.grid
    sim_field = simulate_nasc(scp)
    df = surveydata.acoustics[surveydata.acoustics["class"] == scp.class_name]
    bubble_factor = max_bubblesize / df["nasc"].max()

    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(simdomain["x"], simdomain["y"], c=sim_field, s=markersize,
                    marker="s", **kwargs)
    ax.scatter(df["x"], df["y"], color="white", s=df["nasc"] * bubble_factor,
               alpha=alpha, edgecolors="none")
    ax.set_aspect("equal")
    ax.set_title(str(scp.class_name))
    ax.set_xlabel("Easting (km)")
    ax.set_ylabel("Northing (km)")
    fig.colorbar(sc, ax=ax)
    return fig


def plot_geosim_stats(atbp, surveydata, n=500):
    nclasses = len(atbp.class_problems)
    fig, axes = plt.subplots(3, nclasses, figsize=(6 * nclasses, 12))
    if nclasses == 1:
        axes = np.array(axes).reshape(3, 1)

    for j, prob in enumerate(atbp.class_problems):
        obs_nasc = surveydata.acoustics.loc[surveydata.acoustics["class"] == prob.class_name, "nasc"].values
        sim_nascs = [simulate_nasc(prob) for _ in range(n)]

        qq = np.linspace(0.01, 0.99, 99)
        q_obs = np.quantile(obs_nasc, qq)
        q_sims = np.array([np.quantile(s, qq) for s in sim_nascs])
        ax = axes[0, j]
        ax.plot(q_obs, q_sims.mean(axis=0), marker="o", linestyle="")
        ax.plot([0, np.quantile(obs_nasc, 0.99)], [0, np.quantile(obs_nasc, 0.99)])
        ax.set_title(prob.class_name)
        ax.set_xlabel("Observed")
        ax.set_ylabel("Simulated")

        means = np.array([s.mean() for s in sim_nascs])
        ax = axes[1, j]
        ax.hist(means, bins=30, density=True)
        ax.axvline(obs_nasc.mean(), linewidth=2)
        ax.set_xlabel("Mean NASC")
        ax.set_ylabel("Density")

        sds = np.array([s.std() for s in sim_nascs])
        ax = axes[2, j]
        ax.hist(sds, bins=30, density=True)
        ax.axvline(obs_nasc.std(), linewidth=2)
        ax.set_xlabel("Std. dev. NASC")
        ax.set_ylabel("Density")
    return fig


def plot_boot_results(results, size=(9, 4)):
    pk = results[results["species_code"] == 21740]
    ages = sorted(pk["age"].unique())
    fig, axes = plt.subplots(1, 2, figsize=size)
    axes[0].boxplot([pk[pk["age"] == a]["n"] / 1e9 for a in ages], labels=ages)
    axes[0].set_xlabel("Age class")
    axes[0].set_ylabel("Abundance (billions)")
    axes[1].boxplot([pk[pk["age"] == a]["biomass"] / 1e9 for a in ages], labels=ages)
    axes[1].set_xlabel("Age class")
    axes[1].set_ylabel("Biomass (Mt)")
    return fig


def plot_error_sources(results_totals, xlims=None):
    stds_boot = []
    for _ in range(1000):
        df = resample_df(results_totals)
        g = df.groupby("error_label").agg(
            n_cv=("n", lambda x: x.std() / x.mean()),
            biomass_cv=("biomass", lambda x: x.std() / x.mean()),
        )
        stds_boot.append(g)
    stds_boot = pd.concat(stds_boot).reset_index()
    if xlims is None:
        xmax = max(stds_boot["n_cv"].max(), stds_boot["biomass_cv"].max()) * 1.05
        xlims = (-0.005, xmax)

    fig, axes = plt.subplots(2, 1, figsize=(8, 8))
    axes[0].boxplot([stds_boot[stds_boot["error_label"] == e]["n_cv"]
                     for e in stds_boot["error_label"].unique()],
                    labels=stds_boot["error_label"].unique(), vert=False)
    axes[0].set_xlim(*xlims)
    axes[0].set_xlabel("C.V. (Numbers)")
    axes[1].boxplot([stds_boot[stds_boot["error_label"] == e]["biomass_cv"]
                     for e in stds_boot["error_label"].unique()],
                    labels=stds_boot["error_label"].unique(), vert=False)
    axes[1].set_xlim(*xlims)
    axes[1].set_xlabel("C.V. (Biomass)")
    return fig


def summarize_bootstrap(results, variable="n", species_codes=21740):
    in_spp = results[results["species_code"].isin(np.atleast_1d(species_codes))]
    df = in_spp.melt(id_vars=["i", "species_code", "age"], value_vars=["n", "biomass"],
                     var_name="variable", value_name="value")
    df = df[df["variable"] == variable]
    summary = (
        df.groupby("age", as_index=False)
        .agg(mean=("value", "mean"), std=("value", "std"))
    )
    summary["cv"] = summary["std"] / summary["mean"] * 100
    summary = summary.rename(columns={"mean": variable})
    return summary


def summarize_stepwise_bootstrap(results, variable="n", species_codes=21740):
    in_spp = results[results["species_code"].isin(np.atleast_1d(species_codes))]
    df = in_spp.melt(id_vars=["i", "species_code", "age", "added_error"],
                     value_vars=["n", "biomass"], var_name="variable", value_name="value")
    df = df[df["variable"] == variable]
    summary = (
        df.groupby(["added_error", "age"], as_index=False)
        .agg(mean=("value", "mean"), std=("value", "std"))
    )
    summary["cv"] = summary["std"] / summary["mean"] * 100
    summary = summary.rename(columns={"mean": variable})
    return summary


def plot_error_source_by_age(results_step, results, variable="n", species_codes=21740):
    stepwise_summary = summarize_stepwise_bootstrap(results_step, variable, species_codes)
    results_summary = summarize_bootstrap(results, variable, species_codes)
    fig, ax = plt.subplots(figsize=(8, 5))
    for key, grp in stepwise_summary.groupby("added_error"):
        ax.plot(grp["age"], grp["std"] / 1e9, marker="o", label=key)
    ax.plot(results_summary["age"], results_summary["std"] / 1e9,
            marker="o", linewidth=2, color="black", label="All")
    ax.set_xlabel("Age class")
    ax.set_ylabel("S.D. (Biomass, MT)")
    ax.legend()
    return fig


def merge_results(results, results_step):
    stepwise_totals = results_step.groupby(["added_error", "i"], as_index=False).agg(
        n=("n", "sum"), biomass=("biomass", "sum")
    )
    results_totals = results.groupby("i", as_index=False).agg(
        n=("n", "sum"), biomass=("biomass", "sum")
    )
    results_totals["added_error"] = "All"
    merged = pd.concat([results_totals, stepwise_totals], ignore_index=True)
    return merged
