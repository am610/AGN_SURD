"""Continuous-estimator cross-check for the histogram-based SURD analysis.

This script uses KSG nearest-neighbour estimators to test total mutual
information and conditional mutual information. It does not decompose those
quantities into SURD atoms; instead, it checks whether lag-dependent dependence
and the long-lag null conclusion survive without a joint histogram.
"""

from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/surd-matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import digamma
from sklearn.neighbors import KDTree, NearestNeighbors


ROOT = Path(__file__).resolve().parent
PROCESSED = ROOT / "agn_surd_project" / "processed"
CONTINUUM_PATH = (
    ROOT / "agn_surd_project" / "agn_data" / "ngc5548_agnwatch" / "c5100.dat"
)
LINE_PATH = PROCESSED / "ngc5548_hb_velocity_bins.csv"
MJD_MIN = 47512.0
MJD_MAX = 49255.0
LAGS = np.arange(1, 201)
CONDITIONAL_LAGS = np.arange(1, 121)
K_VALUES = (3, 5, 10)
TARGETS = ("continuum", "blue_wing", "core", "red_wing")
SURD_PEAKS = {"continuum": 115, "blue_wing": 93, "core": 128, "red_wing": 198}


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / values.std(ddof=0)


def load_real_series() -> np.ndarray:
    grid = np.arange(MJD_MIN, MJD_MAX + 1.0)
    continuum = pd.read_csv(
        CONTINUUM_PATH,
        sep=r"\s+",
        header=None,
        names=["mjd", "flux", "err"],
    ).sort_values("mjd")
    line = pd.read_csv(LINE_PATH).sort_values("mjd")
    arrays = [
        np.interp(grid, continuum["mjd"], continuum["flux"]),
        np.interp(grid, line["mjd"], line["blue_wing_flux"]),
        np.interp(grid, line["mjd"], line["core_flux"]),
        np.interp(grid, line["mjd"], line["red_wing_flux"]),
    ]
    return np.vstack([zscore(values) for values in arrays])


def as_2d(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[:, None] if values.ndim == 1 else values


def deterministic_jitter(values: np.ndarray) -> np.ndarray:
    """Break exact interpolation ties without adding stochastic variation."""
    values = np.asarray(values, dtype=float).copy()
    scale = np.maximum(values.std(axis=0, ddof=0), 1.0)
    index = np.arange(len(values), dtype=float)[:, None]
    values += 1e-10 * scale * np.sin(index * 1.61803398875 + np.arange(values.shape[1]))
    return values


def radius_counts(values: np.ndarray, radii: np.ndarray) -> np.ndarray:
    tree = KDTree(values, metric="chebyshev")
    strict_radii = np.nextafter(radii, 0.0)
    return tree.query_radius(values, r=strict_radii, count_only=True) - 1


def ksg_mi_multi(x: np.ndarray, y: np.ndarray, k_values=K_VALUES) -> dict[int, float]:
    x = deterministic_jitter(as_2d(x))
    y = deterministic_jitter(as_2d(y))
    joint = np.column_stack([x, y])
    n_samples = len(joint)
    max_k = max(k_values)
    distances = NearestNeighbors(
        metric="chebyshev", n_neighbors=max_k
    ).fit(joint).kneighbors(return_distance=True)[0]

    estimates = {}
    for k in k_values:
        # kneighbors(X=None) excludes the training sample itself.
        radii = distances[:, k - 1]
        nx = radius_counts(x, radii)
        ny = radius_counts(y, radii)
        estimates[k] = float(
            digamma(k)
            + digamma(n_samples)
            - np.mean(digamma(nx + 1) + digamma(ny + 1))
        )
    return estimates


def ksg_cmi_multi(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    k_values=K_VALUES,
) -> dict[int, float]:
    x = deterministic_jitter(as_2d(x))
    y = deterministic_jitter(as_2d(y))
    z = deterministic_jitter(as_2d(z))
    xz = np.column_stack([x, z])
    yz = np.column_stack([y, z])
    joint = np.column_stack([x, y, z])
    max_k = max(k_values)
    distances = NearestNeighbors(
        metric="chebyshev", n_neighbors=max_k
    ).fit(joint).kneighbors(return_distance=True)[0]

    estimates = {}
    for k in k_values:
        radii = distances[:, k - 1]
        nxz = radius_counts(xz, radii)
        nyz = radius_counts(yz, radii)
        nz = radius_counts(z, radii)
        estimates[k] = float(
            digamma(k)
            + np.mean(
                digamma(nz + 1) - digamma(nxz + 1) - digamma(nyz + 1)
            )
        )
    return estimates


def scan_total_mi(series: np.ndarray) -> pd.DataFrame:
    rows = []
    for target_index, target_name in enumerate(TARGETS):
        predictor_indices = [index for index in range(4) if index != target_index]
        for lag in LAGS:
            estimates = ksg_mi_multi(
                series[predictor_indices, :-lag].T,
                series[target_index, lag:],
            )
            for k, estimate in estimates.items():
                rows.append(
                    {"target": target_name, "lag": int(lag), "k": k, "mi_nats": estimate}
                )
    return pd.DataFrame(rows)


def scan_conditional_mi(series: np.ndarray) -> pd.DataFrame:
    rows = []
    for lag in CONDITIONAL_LAGS:
        estimates = ksg_cmi_multi(
            series[[0, 1], :-lag].T,
            series[2, lag:],
            series[2, :-lag],
        )
        for k, estimate in estimates.items():
            rows.append({"lag": int(lag), "k": k, "cmi_nats": estimate})
    return pd.DataFrame(rows)


_SERIES: np.ndarray | None = None


def initialize_worker(series: np.ndarray) -> None:
    global _SERIES
    _SERIES = series


def run_surrogate(seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    if _SERIES is None:
        raise RuntimeError("Worker series were not initialized")
    rng = np.random.default_rng(seed)
    length = _SERIES.shape[1]
    shifted = np.vstack(
        [np.roll(row, rng.integers(220, length - 220)) for row in _SERIES]
    )
    mi = scan_total_mi(shifted)
    mi["seed"] = seed
    cmi = scan_conditional_mi(shifted)
    cmi["seed"] = seed
    return mi, cmi


def empirical_p(real_value: float, null_values: np.ndarray) -> float:
    return float((1 + np.sum(null_values >= real_value)) / (1 + len(null_values)))


def summarize_mi(real: pd.DataFrame, null: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (target, k), real_group in real.groupby(["target", "k"], sort=True):
        real_group = real_group.sort_values("lag")
        null_group = null[(null["target"] == target) & (null["k"] == k)]
        peak = real_group.loc[real_group["mi_nats"].idxmax()]
        long_real = real_group[real_group["lag"] >= 61]
        long_peak = long_real.loc[long_real["mi_nats"].idxmax()]
        null_max = null_group.groupby("seed")["mi_nats"].max().to_numpy()
        null_long_max = (
            null_group[null_group["lag"] >= 61]
            .groupby("seed")["mi_nats"]
            .max()
            .to_numpy()
        )
        surd_lag = SURD_PEAKS[target]
        at_surd = real_group[real_group["lag"] == surd_lag].iloc[0]
        null_at_surd = null_group[null_group["lag"] == surd_lag]["mi_nats"].to_numpy()
        rows.append(
            {
                "target": target,
                "k": k,
                "real_peak_lag_1_200": int(peak["lag"]),
                "real_peak_mi_nats_1_200": peak["mi_nats"],
                "global_p_1_200": empirical_p(peak["mi_nats"], null_max),
                "real_peak_lag_61_200": int(long_peak["lag"]),
                "real_peak_mi_nats_61_200": long_peak["mi_nats"],
                "global_p_61_200": empirical_p(long_peak["mi_nats"], null_long_max),
                "surd_peak_lag": surd_lag,
                "mi_nats_at_surd_peak": at_surd["mi_nats"],
                "pointwise_p_at_surd_peak": empirical_p(
                    at_surd["mi_nats"], null_at_surd
                ),
                "surrogate_count": len(null_max),
            }
        )
    return pd.DataFrame(rows)


def summarize_cmi(real: pd.DataFrame, null: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for k, real_group in real.groupby("k", sort=True):
        null_group = null[null["k"] == k]
        peak = real_group.loc[real_group["cmi_nats"].idxmax()]
        null_max = null_group.groupby("seed")["cmi_nats"].max().to_numpy()
        default_lag = 73
        at_default = real_group[real_group["lag"] == default_lag].iloc[0]
        null_at_default = null_group[null_group["lag"] == default_lag]["cmi_nats"].to_numpy()
        rows.append(
            {
                "k": k,
                "real_peak_lag_1_120": int(peak["lag"]),
                "real_peak_cmi_nats": peak["cmi_nats"],
                "global_p_1_120": empirical_p(peak["cmi_nats"], null_max),
                "cmi_nats_at_histogram_peak_73d": at_default["cmi_nats"],
                "pointwise_p_at_73d": empirical_p(
                    at_default["cmi_nats"], null_at_default
                ),
                "surrogate_count": len(null_max),
            }
        )
    return pd.DataFrame(rows)


def plot_crosscheck(
    real_mi: pd.DataFrame,
    null_mi: pd.DataFrame,
    real_cmi: pd.DataFrame,
    null_cmi: pd.DataFrame,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharex=False)
    axes = axes.flat
    k = 5

    for axis, target in zip(axes[:4], TARGETS):
        real = real_mi[(real_mi["target"] == target) & (real_mi["k"] == k)]
        null = null_mi[(null_mi["target"] == target) & (null_mi["k"] == k)]
        envelope = null.groupby("lag")["mi_nats"].agg(
            median="median",
            low=lambda values: values.quantile(0.025),
            high=lambda values: values.quantile(0.975),
        )
        axis.fill_between(
            envelope.index,
            envelope["low"],
            envelope["high"],
            color="#d8c3a5",
            alpha=0.45,
            label="95% shift-null envelope",
        )
        axis.plot(envelope.index, envelope["median"], color="#a66f2c", linestyle="--")
        axis.plot(real["lag"], real["mi_nats"], color="#174a7e", linewidth=1.8,
                  label="Real KNN MI")
        axis.axvline(SURD_PEAKS[target], color="#9b2226", linestyle=":",
                     linewidth=1.5, label="SURD maximum")
        axis.set_title(target.replace("_", " ").title())
        axis.set_xlabel("Lag (days)")
        axis.set_ylabel("Mutual information (nats)")

    real = real_cmi[real_cmi["k"] == k]
    null = null_cmi[null_cmi["k"] == k]
    envelope = null.groupby("lag")["cmi_nats"].agg(
        median="median",
        low=lambda values: values.quantile(0.025),
        high=lambda values: values.quantile(0.975),
    )
    axis = axes[4]
    axis.fill_between(envelope.index, envelope["low"], envelope["high"],
                      color="#d8c3a5", alpha=0.45)
    axis.plot(envelope.index, envelope["median"], color="#a66f2c", linestyle="--")
    axis.plot(real["lag"], real["cmi_nats"], color="#2a6f4e", linewidth=1.8)
    axis.axvline(73, color="#9b2226", linestyle=":", linewidth=1.5)
    axis.set_title("Core Conditional MI")
    axis.set_xlabel("Lag (days)")
    axis.set_ylabel("Conditional MI (nats)")

    axes[5].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[5].legend(handles, labels, loc="center", frameon=False, fontsize=11)
    fig.tight_layout()
    fig.savefig(ROOT / "overleaf_draft" / "figure10_knn_crosscheck.png", dpi=300)
    plt.close(fig)


def validate_estimators() -> None:
    rng = np.random.default_rng(42)
    n = 4000
    x = rng.normal(size=n)
    y = 0.7 * x + np.sqrt(1 - 0.7**2) * rng.normal(size=n)
    expected = -0.5 * np.log(1 - 0.7**2)
    estimated = ksg_mi_multi(x, y, (5,))[5]

    z = rng.normal(size=n)
    x_cond = z + 0.5 * rng.normal(size=n)
    y_cond = z + 0.5 * rng.normal(size=n)
    null_cmi = ksg_cmi_multi(x_cond, y_cond, z, (5,))[5]
    direct_cmi = ksg_cmi_multi(x_cond, y_cond + 0.5 * x_cond, z, (5,))[5]
    print(
        f"Estimator checks: Gaussian MI expected={expected:.3f}, estimated={estimated:.3f}; "
        f"common-driver CMI={null_cmi:.3f}; direct-link CMI={direct_cmi:.3f}"
    )
    if abs(estimated - expected) > 0.08 or abs(null_cmi) > 0.08 or direct_cmi <= 0.05:
        raise RuntimeError("KNN estimator validation failed")


def run_crosscheck(surrogates: int, workers: int) -> None:
    validate_estimators()
    series = load_real_series()
    print("Running real KNN mutual-information scans...", flush=True)
    real_mi = scan_total_mi(series)
    real_cmi = scan_conditional_mi(series)
    real_mi.to_csv(PROCESSED / "knn_mi_real_curves.csv", index=False)
    real_cmi.to_csv(PROCESSED / "knn_cmi_real_curves.csv", index=False)

    print(f"Running {surrogates} circular-shift surrogate scans...", flush=True)
    mi_frames = []
    cmi_frames = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_worker,
        initargs=(series,),
    ) as executor:
        for index, (mi, cmi) in enumerate(
            executor.map(run_surrogate, range(9000, 9000 + surrogates)), start=1
        ):
            mi_frames.append(mi)
            cmi_frames.append(cmi)
            if index % 10 == 0 or index == surrogates:
                print(f"  surrogate scans complete: {index}/{surrogates}", flush=True)
    null_mi = pd.concat(mi_frames, ignore_index=True)
    null_cmi = pd.concat(cmi_frames, ignore_index=True)
    null_mi.to_csv(PROCESSED / "knn_mi_surrogate_curves.csv", index=False)
    null_cmi.to_csv(PROCESSED / "knn_cmi_surrogate_curves.csv", index=False)

    mi_summary = summarize_mi(real_mi, null_mi)
    cmi_summary = summarize_cmi(real_cmi, null_cmi)
    mi_summary.to_csv(PROCESSED / "knn_mi_summary.csv", index=False)
    cmi_summary.to_csv(PROCESSED / "knn_cmi_summary.csv", index=False)
    plot_crosscheck(real_mi, null_mi, real_cmi, null_cmi)
    print("\nTotal MI summary")
    print(mi_summary.to_string(index=False))
    print("\nConditional MI summary")
    print(cmi_summary.to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--surrogates", type=int, default=999)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run_crosscheck(arguments.surrogates, arguments.workers)
