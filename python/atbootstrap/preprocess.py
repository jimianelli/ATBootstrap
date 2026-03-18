"""Preprocessing utilities for ATBootstrap Python port."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
import alphashape
from pyproj import Transformer

from .spatial import in_intervals


class AbstractSurveyDomain:
    pass


@dataclass
class TransectRibbons(AbstractSurveyDomain):
    transect_width: float = 20.0
    buffer: float = 0.1


@dataclass
class SurveyHull(AbstractSurveyDomain):
    k: int = 10


def transect_ribbon(transect: pd.DataFrame, transect_width: float, dx: float,
                    buffer: float = 0.1, order: str = "y") -> Polygon:
    tr1 = (
        transect.loc[:, ["x", "y", "log"]]
        .melt(id_vars=[order], value_vars=[c for c in ["x", "y"] if c != order])
    )
    tr1[order] = (tr1[order] / dx).round() * dx
    tr1 = (
        tr1.groupby([order, "variable"], as_index=False)["value"].mean()
        .pivot(index=order, columns="variable", values="value")
        .reset_index()
        .sort_values(order)
        .loc[:, ["x", "y"]]
        .drop_duplicates()
    )

    X = tr1[["x", "y"]].to_numpy().T
    v = np.diff(X, axis=1)
    v = v / np.linalg.norm(v, axis=0)
    w = transect_width / 2 * 1.852 * (1 + buffer)
    v = v * w
    v = np.column_stack([v, v[:, -1]])
    # rotation matrices
    a = np.pi / 2
    R1 = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    R2 = np.array([[np.cos(-a), -np.sin(-a)], [np.sin(-a), np.cos(-a)]])
    left = R1 @ v + X
    right = R2 @ v + X

    ribbon_bounds = [tuple(pt) for pt in np.column_stack([left, np.flip(right, axis=1)]).T]
    ribbon_bounds.append(tuple(left[:, 0]))
    return Polygon(ribbon_bounds)


def survey_domain(acoustics: pd.DataFrame, method: AbstractSurveyDomain,
                  order: str, dx: float, dy: float | None = None):
    dy = dy or dx
    if isinstance(method, TransectRibbons):
        polys = []
        for tr_id in acoustics["transect"].unique():
            tr = acoustics.loc[acoustics["transect"] == tr_id]
            polys.append(transect_ribbon(tr, method.transect_width, dx, method.buffer, order))
        return unary_union(polys)
    if isinstance(method, SurveyHull):
        unique_points = acoustics.loc[:, ["x", "y"]].drop_duplicates()
        pts = unique_points.to_numpy()
        hull = alphashape.alphashape(pts, method.k)
        return hull
    raise ValueError("Unknown survey domain method")


def grid_domain(domain, dx: float, dy: float | None = None) -> pd.DataFrame:
    dy = dy or dx
    minx, miny, maxx, maxy = domain.bounds
    xgrid = np.arange(minx, maxx + dx, dx)
    ygrid = np.arange(miny, maxy + dy, dy)
    rows = []
    for x in xgrid:
        for y in ygrid:
            p = Point(x, y)
            if domain.contains(p):
                rows.append((x, y))
    return pd.DataFrame(rows, columns=["x", "y"])


def get_survey_grid(acoustics: pd.DataFrame, method: AbstractSurveyDomain | None = None,
                    dx: float = 10.0, dy: float | None = None, order: str = "y"):
    if method is None:
        method = TransectRibbons()
    if len(acoustics) <= 3:
        raise ValueError("Not enough locations to define survey grid.")
    domain = survey_domain(acoustics, method, order, dx, dy or dx)
    grid = grid_domain(domain, dx, dy or dx)
    return grid, domain


def make_haul_id(prefix, ship, event_id):
    return prefix + "-" + ship.astype(str) + "-" + event_id.astype(str)


def merge_scaling(scaling_mace: pd.DataFrame, scaling_gap: pd.DataFrame) -> pd.DataFrame:
    ts_key = (
        scaling_mace.groupby("species_code", as_index=False)["ts_relationship"].first()
        .rename(columns={"ts_relationship": "ts_relationship"})
    )

    scaling_mace1 = scaling_mace.copy()
    scaling_mace1["haul_id"] = make_haul_id("MACE", pd.Series(1, index=scaling_mace1.index),
                                            scaling_mace1["event_id"])
    scaling_mace1 = scaling_mace1.loc[:, ["survey", "ship", "haul_id", "class", "species_code",
                                          "primary_length", "ts_length", "ts_relationship",
                                          "catch_sampling_expansion", "user_defined_expansion",
                                          "sample_correction_scalar", "haul_weight", "w"]]

    scaling_gap1 = scaling_gap.copy()
    scaling_gap1["survey"] = scaling_gap1["cruise"]
    scaling_gap1["ship"] = scaling_gap1["vessel"]
    scaling_gap1["haul_id"] = make_haul_id("GAP", scaling_gap1["vessel"], scaling_gap1["haul"])
    scaling_gap1["class"] = "BT"
    scaling_gap1["primary_length"] = scaling_gap1["length"] / 10.0
    scaling_gap1["ts_length"] = scaling_gap1["length"] / 10.0
    scaling_gap1["catch_sampling_expansion"] = 1.0
    scaling_gap1["sample_correction_scalar"] = 1.0
    scaling_gap1["user_defined_expansion"] = 1.0
    scaling_gap1["haul_weight"] = 1.0
    scaling_gap1["w"] = 1.0
    scaling_gap1 = scaling_gap1.merge(ts_key, on="species_code", how="left")
    scaling_gap1["ts_relationship"] = scaling_gap1["ts_relationship"].fillna(
        "generic_swimbladder_fish"
    )
    scaling_gap1 = scaling_gap1.loc[:, ["survey", "ship", "haul_id", "class", "species_code",
                                         "primary_length", "ts_length", "ts_relationship",
                                         "catch_sampling_expansion", "user_defined_expansion",
                                         "sample_correction_scalar", "haul_weight", "w"]]
    scaling_gap1 = scaling_gap1.dropna()

    return pd.concat([scaling_mace1, scaling_gap1], ignore_index=True)


def merge_trawl_locations(trawl_locations_mace: pd.DataFrame,
                           trawl_locations_gap: pd.DataFrame) -> pd.DataFrame:
    survey = trawl_locations_mace["survey"].unique()[0]
    trawl_locations_gap = trawl_locations_gap.copy()
    trawl_locations_gap["survey"] = survey
    tl_mace = trawl_locations_mace.copy()
    tl_mace["haul_id"] = make_haul_id("MACE", pd.Series(1, index=tl_mace.index), tl_mace["event_id"])
    tl_mace = tl_mace.loc[:, ["survey", "haul_id", "latitude", "longitude"]]
    tl_gap = trawl_locations_gap.copy()
    tl_gap["haul_id"] = make_haul_id("GAP", tl_gap["vessel"], tl_gap["event_id"])
    tl_gap = tl_gap.loc[:, ["survey", "haul_id", "latitude", "longitude"]]

    return pd.concat([tl_mace, tl_gap], ignore_index=True)


def preprocess_survey_data(
    surveydir: str,
    ebs: bool = True,
    log_ranges: Sequence[tuple] | None = None,
    dx: float = 10.0,
    dy: float | None = None,
    missingstring: Sequence[str] = (".", "NA"),
    grid_method: AbstractSurveyDomain | None = None,
    transect_order: str = "y",
):
    dy = dy or dx
    scaling_mace = pd.read_csv(f"{surveydir}/scaling_mace.csv", na_values=list(missingstring))
    trawl_locations_mace = pd.read_csv(
        f"{surveydir}/trawl_locations_mace.csv", na_values=list(missingstring)
    )

    if ebs:
        scaling_gap = pd.read_csv(f"{surveydir}/scaling_gap.csv", na_values=list(missingstring))
        scaling = merge_scaling(scaling_mace, scaling_gap)
        trawl_locations_gap = pd.read_csv(
            f"{surveydir}/trawl_locations_gap.csv", na_values=list(missingstring)
        )
        trawl_locations = merge_trawl_locations(trawl_locations_mace, trawl_locations_gap)
    else:
        scaling = scaling_mace
        trawl_locations = trawl_locations_mace

    acoustics = pd.read_csv(f"{surveydir}/acoustics.csv", na_values=list(missingstring))

    scaling["class"] = scaling["class"].str.replace("_FILTERED", "", regex=False)
    acoustics["class"] = acoustics["class"].str.replace("_FILTERED", "", regex=False)
    scaling_classes = scaling["class"].unique()

    if log_ranges is None:
        log_ranges = [tuple(acoustics["start_vessel_log"].agg(["min", "max"]).values)]

    acoustics = (
        acoustics.loc[:, ["transect", "interval", "class", "start_latitude",
                          "start_longitude", "start_vessel_log", "nasc"]]
        .rename(columns={"start_latitude": "lat", "start_longitude": "lon",
                         "start_vessel_log": "log"})
    )
    acoustics = acoustics[
        acoustics["class"].isin(scaling_classes)
        & acoustics["log"].apply(lambda x: in_intervals(x, log_ranges))
        & (acoustics["lon"].abs() < 360)
        & (acoustics["lat"].abs() < 360)
    ]
    acoustics = (
        acoustics.groupby(["transect", "interval", "class", "lat", "lon", "log"], as_index=False)
        .agg(nasc=("nasc", "sum"))
        .pivot_table(index=["transect", "interval", "lat", "lon", "log"],
                     columns="class", values="nasc", fill_value=0)
        .reset_index()
        .melt(id_vars=["transect", "interval", "lat", "lon", "log"],
              var_name="class", value_name="nasc")
    )
    acoustics["nasc"] = acoustics["nasc"].fillna(0)

    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32603", always_xy=True)
    xs, ys = transformer.transform(acoustics["lon"].values, acoustics["lat"].values)
    acoustics["x"] = xs / 1e3
    acoustics["y"] = ys / 1e3

    surveygrid, surveyhull = get_survey_grid(acoustics, method=grid_method,
                                             dx=dx, dy=dy, order=transect_order)

    acoustics["x"] = (acoustics["x"] / dx).round() * dx
    acoustics["y"] = (acoustics["y"] / dy).round() * dy
    acoustics = (
        acoustics.groupby(["transect", "class", "x", "y"], as_index=False)
        .agg(lon=("lon", "mean"), lat=("lat", "mean"), log=("log", "mean"), nasc=("nasc", "mean"))
    )

    trawl_locations = trawl_locations.copy()
    xs, ys = transformer.transform(trawl_locations["longitude"].values,
                                   trawl_locations["latitude"].values)
    trawl_locations["x"] = xs / 1e3
    trawl_locations["y"] = ys / 1e3
    trawl_locations = trawl_locations[
        trawl_locations.apply(lambda r: surveyhull.contains(Point(r["x"], r["y"])), axis=1)
    ]

    length_weight = pd.read_csv(f"{surveydir}/measurements.csv", na_values=list(missingstring))
    length_weight.columns = [c.lower() for c in length_weight.columns]
    length_weight = (
        length_weight.pivot_table(index=["specimen_id", "species_code", "event_id"],
                                  columns="measurement_type", values="measurement_value")
        .reset_index()
        .dropna()
    )

    scaling.to_csv(f"{surveydir}/scaling.csv", index=False)
    acoustics.to_csv(f"{surveydir}/acoustics_projected.csv", index=False)
    length_weight.to_csv(f"{surveydir}/length_weight.csv", index=False)
    trawl_locations.to_csv(f"{surveydir}/trawl_locations_projected.csv", index=False)
    surveygrid.to_csv(f"{surveydir}/surveygrid.csv", index=False)
