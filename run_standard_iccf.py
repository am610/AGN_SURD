"""Run a standard bidirectional ICCF with FR/RSS lag uncertainties."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/surd-matplotlib")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PROCESSED = ROOT / "agn_surd_project" / "processed"
LAGS = np.arange(-100.0, 201.0, 1.0)
N_FRRSS = int(os.environ.get("ICCF_N_FRRSS", "5000"))
CENTROID_THRESHOLD = 0.8
MIN_PAIRS = 11


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < MIN_PAIRS:
        return np.nan
    x = x - np.mean(x)
    y = y - np.mean(y)
    denom = np.sqrt(np.dot(x, x) * np.dot(y, y))
    return float(np.dot(x, y) / denom) if denom > 0 else np.nan


def bidirectional_iccf(
    t_cont: np.ndarray,
    f_cont: np.ndarray,
    t_line: np.ndarray,
    f_line: np.ndarray,
    lags: np.ndarray = LAGS,
) -> np.ndarray:
    """Average both interpolation directions; positive lag means line follows continuum."""
    values = np.full(len(lags), np.nan)
    for index, lag in enumerate(lags):
        directional = []

        line_times = t_cont + lag
        valid = (line_times >= t_line[0]) & (line_times <= t_line[-1])
        if np.count_nonzero(valid) >= MIN_PAIRS:
            directional.append(
                _pearson(f_cont[valid], np.interp(line_times[valid], t_line, f_line))
            )

        cont_times = t_line - lag
        valid = (cont_times >= t_cont[0]) & (cont_times <= t_cont[-1])
        if np.count_nonzero(valid) >= MIN_PAIRS:
            directional.append(
                _pearson(np.interp(cont_times[valid], t_cont, f_cont), f_line[valid])
            )

        finite = np.asarray(directional)[np.isfinite(directional)]
        if finite.size:
            values[index] = np.mean(finite)
    return values


def peak_and_centroid(
    lags: np.ndarray, ccf: np.ndarray, threshold: float = CENTROID_THRESHOLD
) -> tuple[float, float, float]:
    """Return peak lag, contiguous 0.8-rmax centroid, and rmax."""
    if not np.any(np.isfinite(ccf)):
        return np.nan, np.nan, np.nan
    peak_index = int(np.nanargmax(ccf))
    if peak_index in (0, len(ccf) - 1):
        return np.nan, np.nan, np.nan
    rmax = float(ccf[peak_index])
    if rmax < 0.2:
        return np.nan, np.nan, rmax

    cutoff = threshold * rmax
    left = peak_index
    right = peak_index
    while left > 0 and np.isfinite(ccf[left - 1]) and ccf[left - 1] >= cutoff:
        left -= 1
    while right < len(ccf) - 1 and np.isfinite(ccf[right + 1]) and ccf[right + 1] >= cutoff:
        right += 1
    weights = ccf[left : right + 1]
    centroid = np.sum(lags[left : right + 1] * weights) / np.sum(weights)
    return float(lags[peak_index]), float(centroid), rmax


def _rss_fr(
    time: np.ndarray, flux: np.ndarray, error: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    indices, counts = np.unique(rng.integers(0, len(time), len(time)), return_counts=True)
    rss_error = error[indices] / np.sqrt(counts)
    return time[indices], rng.normal(flux[indices], rss_error)


def _single_frrss(args: tuple) -> tuple[float, float, float]:
    t_cont, f_cont, e_cont, t_line, f_line, e_line, seed = args
    rng = np.random.default_rng(seed)
    tc, fc = _rss_fr(t_cont, f_cont, e_cont, rng)
    tl, fl = _rss_fr(t_line, f_line, e_line, rng)
    return peak_and_centroid(LAGS, bidirectional_iccf(tc, fc, tl, fl))


def _interval(values: np.ndarray) -> tuple[float, float, float]:
    finite = values[np.isfinite(values)]
    return tuple(np.percentile(finite, [16, 50, 84])) if finite.size else (np.nan,) * 3


def load_native_light_curves():
    continuum = pd.read_csv(
        ROOT / "agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat",
        sep=r"\s+",
        header=None,
        names=["time", "flux", "error"],
    ).sort_values("time")
    lines = pd.read_csv(PROCESSED / "ngc5548_hb_velocity_bins.csv").sort_values("mjd")
    lower = max(continuum.time.min(), lines.mjd.min())
    upper = min(continuum.time.max(), lines.mjd.max())
    continuum = continuum[continuum.time.between(lower, upper)].reset_index(drop=True)
    lines = lines[lines.mjd.between(lower, upper)].reset_index(drop=True)
    return continuum, lines, lower, upper


def main():
    continuum, lines, lower, upper = load_native_light_curves()
    components = {
        "blue_wing": ("blue_wing_flux", "blue_wing_error"),
        "core": ("core_flux", "core_error"),
        "red_wing": ("red_wing_flux", "red_wing_error"),
    }
    curve_rows = []
    distribution_rows = []
    summary_rows = []

    for component_index, (name, (flux_col, error_col)) in enumerate(components.items()):
        ccf = bidirectional_iccf(
            continuum.time.to_numpy(),
            continuum.flux.to_numpy(),
            lines.mjd.to_numpy(),
            lines[flux_col].to_numpy(),
        )
        peak, centroid, rmax = peak_and_centroid(LAGS, ccf)
        jobs = [
            (
                continuum.time.to_numpy(),
                continuum.flux.to_numpy(),
                continuum.error.to_numpy(),
                lines.mjd.to_numpy(),
                lines[flux_col].to_numpy(),
                lines[error_col].to_numpy(),
                20260716 + component_index * 100_000 + index,
            )
            for index in range(N_FRRSS)
        ]
        with ProcessPoolExecutor() as executor:
            samples = np.asarray(list(executor.map(_single_frrss, jobs, chunksize=20)))

        peak16, peak50, peak84 = _interval(samples[:, 0])
        cent16, cent50, cent84 = _interval(samples[:, 1])
        success = np.isfinite(samples[:, 1])
        summary_rows.append(
            {
                "component": name,
                "native_peak": peak,
                "native_centroid": centroid,
                "r_max": rmax,
                "peak_p16": peak16,
                "peak_median": peak50,
                "peak_p84": peak84,
                "centroid_p16": cent16,
                "centroid_median": cent50,
                "centroid_p84": cent84,
                "centroid_err_minus": cent50 - cent16,
                "centroid_err_plus": cent84 - cent50,
                "n_frrss": N_FRRSS,
                "n_success": int(success.sum()),
                "success_fraction": float(success.mean()),
                "time_min": lower,
                "time_max": upper,
            }
        )
        curve_rows.extend(
            {"component": name, "lag": lag, "r": value}
            for lag, value in zip(LAGS, ccf)
        )
        distribution_rows.extend(
            {
                "component": name,
                "realization": index,
                "peak_lag": row[0],
                "centroid_lag": row[1],
                "r_max": row[2],
            }
            for index, row in enumerate(samples)
        )
        print(
            f"{name}: centroid={cent50:.2f} -{cent50-cent16:.2f}/+{cent84-cent50:.2f} d "
            f"(native {centroid:.2f} d, rmax={rmax:.3f}, success={success.mean():.3f})"
        )

    pd.DataFrame(curve_rows).to_csv(PROCESSED / "iccf_curves.csv", index=False)
    pd.DataFrame(distribution_rows).to_csv(PROCESSED / "iccf_frrss_distributions.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(PROCESSED / "iccf_summary.csv", index=False)


if __name__ == "__main__":
    main()
