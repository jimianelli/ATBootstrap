"""Spatial and geostatistical utilities for ATBootstrap Python port."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree, distance_matrix
from scipy.optimize import curve_fit
from scipy.stats import gamma, invgauss, invgamma, lognorm


def in_intervals(x: float, intervals: Sequence[tuple]) -> bool:
    for lo, hi in intervals:
        if lo <= x <= hi:
            return True
    return False


def _empirical_variogram(coords: np.ndarray, values: np.ndarray,
                          nlags: int = 10, maxlag: float = 200.0):
    dists = distance_matrix(coords, coords)
    diffs = values[:, None] - values[None, :]
    semi = 0.5 * (diffs ** 2)
    # upper triangle
    iu = np.triu_indices_from(dists, k=1)
    d = dists[iu]
    g = semi[iu]
    bins = np.linspace(0, maxlag, nlags + 1)
    bin_idx = np.digitize(d, bins) - 1
    centers = 0.5 * (bins[:-1] + bins[1:])
    emp = np.array([g[bin_idx == i].mean() if np.any(bin_idx == i) else np.nan
                    for i in range(nlags)])
    return centers, emp


def _exp_variogram(h, nugget, sill, rng):
    return nugget + sill * (1 - np.exp(-h / rng))


def define_conditional_sim(acoustics: pd.DataFrame, sim_domain: pd.DataFrame,
                            maxlag: float = 200.0, nlags: int = 10,
                            weightfunc: Callable[[float], float] | None = None):
    geonasc = acoustics.loc[:, ["nasc", "x", "y"]].copy()
    geonasc["x"] += 1e-3 * np.random.randn(len(geonasc))
    geonasc["y"] += 1e-3 * np.random.randn(len(geonasc))

    coords = geonasc[["x", "y"]].to_numpy()
    values = geonasc["nasc"].to_numpy()
    h, emp = _empirical_variogram(coords, values, nlags=nlags, maxlag=maxlag)

    # weighted least squares fit
    valid = np.isfinite(emp)
    weights = None
    if weightfunc is not None:
        weights = np.array([weightfunc(x) for x in h[valid]])
    p0 = [0.0, np.nanmax(emp[valid]), maxlag / 3]
    popt, _ = curve_fit(_exp_variogram, h[valid], emp[valid], p0=p0,
                        sigma=None if weights is None else 1 / weights,
                        maxfev=10000)
    nugget, sill, rng = popt
    variogram = {
        "empirical": {"lags": h, "gamma": emp},
        "model": {"nugget": nugget, "sill": sill, "range": rng},
    }
    setup = {"data": geonasc, "domain": sim_domain}
    return variogram, setup


def _covariance_from_variogram(h, nugget, sill, rng):
    # sill is partial sill; total variance = nugget + sill
    gamma = _exp_variogram(h, nugget, sill, rng)
    return (nugget + sill) - gamma


def get_lungs_params(setup: dict, theoretical_variogram: dict, variable: str = "nasc"):
    data = setup["data"]
    domain = setup["domain"]

    dcoords = data[["x", "y"]].to_numpy()
    scoords = domain[["x", "y"]].to_numpy()
    z = data[variable].to_numpy()

    nugget = theoretical_variogram["nugget"]
    sill = theoretical_variogram["sill"]
    rng = theoretical_variogram["range"]

    # covariance matrices
    Cdd = _covariance_from_variogram(distance_matrix(dcoords, dcoords), nugget, sill, rng)
    Csd = _covariance_from_variogram(distance_matrix(scoords, dcoords), nugget, sill, rng)
    Css = _covariance_from_variogram(distance_matrix(scoords, scoords), nugget, sill, rng)

    # conditional mean and covariance (ordinary kriging with constant mean)
    mu = np.mean(z)
    zc = z - mu
    Cdd_inv = np.linalg.pinv(Cdd)
    mu_x = mu + Csd @ Cdd_inv @ zc
    Ccond = Css - Csd @ Cdd_inv @ Csd.T
    # jitter for numerical stability
    jitter = 1e-10
    for _ in range(10):
        try:
            L = np.linalg.cholesky(Ccond + np.eye(Ccond.shape[0]) * jitter)
            break
        except np.linalg.LinAlgError:
            jitter *= 10
    else:
        # final fallback: shift by absolute min eigenvalue
        eigmin = np.min(np.linalg.eigvalsh(Ccond))
        shift = abs(eigmin) + 1e-6
        L = np.linalg.cholesky(Ccond + np.eye(Ccond.shape[0]) * shift)

    return {
        "data": z,
        "mu_x": mu_x,
        "L": L,
        "mu": mu,
        "dlocs": np.arange(len(z)),
        "slocs": np.arange(len(mu_x)),
    }


def _dist_params(dist_name: str, mu: float, var: float):
    if dist_name == "Gamma":
        k = var / mu
        theta = mu ** 2 / var
        return k, theta
    if dist_name == "InverseGaussian":
        return mu, mu ** 3 / var
    if dist_name == "InverseGamma":
        return mu ** 2 / var + 2, mu ** 3 / var + mu
    if dist_name == "LogNormal":
        sigma = np.sqrt(np.log(var / np.exp(2 * np.log(mu)) + 1))
        mu_l = np.log(mu) - sigma ** 2 / 2
        return mu_l, sigma
    raise ValueError("Unsupported distribution")


def parameterize_zdists(dist_name: str, lungs_params: dict, eps: float | None = None):
    eps = eps or np.cbrt(np.finfo(float).eps)
    mu_x = lungs_params["mu_x"].copy()
    mu_x[mu_x < eps] = eps
    mu_x = mu_x * lungs_params["data"].mean() / mu_x.mean()
    mu_z = np.linalg.solve(lungs_params["L"], mu_x)
    mu_z[mu_z <= eps] = eps
    vz = np.ones_like(mu_z)

    params = [_dist_params(dist_name, m, v) for m, v in zip(mu_z, vz)]
    return params


def zdist_mean(dist_name: str, params):
    if dist_name == "Gamma":
        k, theta = params
        return k * theta
    if dist_name == "InverseGaussian":
        mu, lam = params
        return mu
    if dist_name == "InverseGamma":
        a, scale = params
        return scale / (a - 1) if a > 1 else np.nan
    if dist_name == "LogNormal":
        mu_l, sigma = params
        return np.exp(mu_l + 0.5 * sigma ** 2)
    raise ValueError("Unsupported distribution")


def nonneg_lumult(params: dict, z: np.ndarray):
    # Simulated values at simulation locations only
    return params["L"] @ z


def nonneg_lusim(params: dict, zdists: list, dist_name: str):
    rng = np.random.default_rng()
    if dist_name == "Gamma":
        z = np.array([gamma.rvs(a=p[0], scale=p[1], random_state=rng) for p in zdists])
    elif dist_name == "InverseGaussian":
        z = np.array([invgauss.rvs(mu=p[0], scale=p[1], random_state=rng) for p in zdists])
    elif dist_name == "InverseGamma":
        z = np.array([invgamma.rvs(a=p[0], scale=p[1], random_state=rng) for p in zdists])
    elif dist_name == "LogNormal":
        z = np.array([lognorm.rvs(s=p[1], scale=np.exp(p[0]), random_state=rng) for p in zdists])
    else:
        raise ValueError("Unsupported distribution")
    return nonneg_lumult(params, z)


def simulate_nasc(scp):
    return nonneg_lusim(scp.params, scp.zdists, scp.zfamily)


def choose_z_distribution(candidate_dists, nasc: np.ndarray, lungs_params: dict,
                          nreplicates: int = 500, verbose: bool = False):
    bin_edges = np.concatenate([[0], 2 ** np.arange(0, 15)])
    h_nasc, _ = np.histogram(nasc, bins=bin_edges, density=True)
    fits = []
    for dist_name in candidate_dists:
        if verbose:
            print(f"Comparing with {dist_name}...")
        zdists = parameterize_zdists(dist_name, lungs_params)
        for _ in range(nreplicates):
            x = nonneg_lusim(lungs_params, zdists, dist_name)
            h_sim, _ = np.histogram(x, bins=bin_edges, density=True)
            # KL divergence
            mask = (h_nasc > 0) & (h_sim > 0)
            kld = np.sum(h_nasc[mask] * np.log(h_nasc[mask] / h_sim[mask]))
            fits.append({"distribution": dist_name, "kld": kld})
    df = pd.DataFrame(fits)
    df = df[np.isfinite(df["kld"])]
    dist_fits = df.groupby("distribution").agg(mean_kld=("kld", "mean"),
                                                se_kld=("kld", "sem")).reset_index()
    return dist_fits.loc[dist_fits["mean_kld"].idxmin(), "distribution"]


def zdists(atbp):
    return pd.DataFrame(
        [{"class": cp.class_name, "zdist": cp.zfamily} for cp in atbp.class_problems]
    )


def trawl_assignments(pixel_coords: np.ndarray, trawl_coords: np.ndarray,
                      stochastic: bool = True, nneighbors: int = 4, a: float = 2.0):
    nneighbors = min(nneighbors, len(trawl_coords))
    kdtree = cKDTree(trawl_coords)
    dists, idx = kdtree.query(pixel_coords, k=nneighbors)
    if nneighbors == 1:
        dists = dists[:, None]
        idx = idx[:, None]
    assignments = np.zeros(len(pixel_coords), dtype=int)
    rng = np.random.default_rng()
    for i in range(len(assignments)):
        if stochastic:
            weights = dists[i] ** (-a)
            weights = weights / weights.sum()
            assignments[i] = rng.choice(idx[i], p=weights)
        else:
            assignments[i] = idx[i][np.argmin(dists[i])]
    return assignments
