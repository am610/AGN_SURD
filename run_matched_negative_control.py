"""Matched three-predictor red-noise control and observing-window ablation."""

from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/surd-matplotlib")

import numpy as np
import pandas as pd
from scipy.stats import zscore


ROOT = Path(__file__).resolve().parent
PROCESSED = ROOT / "agn_surd_project" / "processed"
sys.path.extend([str(ROOT / "SURD"), str(ROOT / "SURD" / "utils")])
from utils import surd


LAGS = np.arange(1, 201)
NBINS = 8
N_REALIZATIONS = int(os.environ.get("ABLATION_N_REALIZATIONS", "200"))


def simulate_ou(n: int, timescale: float, rng: np.random.Generator) -> np.ndarray:
    phi = np.exp(-1.0 / timescale)
    noise_scale = np.sqrt(1.0 - phi * phi)
    values = np.empty(n + 500)
    values[0] = rng.normal()
    for index in range(1, len(values)):
        values[index] = phi * values[index - 1] + noise_scale * rng.normal()
    return zscore(values[500:])


def surd_metrics(target, predictors, lag, row_mask=None):
    future = target[lag:]
    current = [predictor[:-lag] for predictor in predictors]
    if row_mask is not None:
        future = future[row_mask]
        current = [predictor[row_mask] for predictor in current]
    data = np.column_stack([future, *current])
    if len(data) < 10:
        return np.nan, np.nan, np.nan, len(data), np.nan
    hist, _ = np.histogramdd(data, bins=NBINS)
    hist /= hist.sum()
    empty_fraction = float(np.mean(hist == 0))
    i_r, i_s, mi, _ = surd.surd(hist)
    joint_mi = mi.get((1, 2, 3), np.nan)
    synergy = sum(i_s.get(key, 0.0) for key in ((1, 2), (1, 3), (2, 3), (1, 2, 3)))
    normalized = synergy / joint_mi if joint_mi > 0 else np.nan
    return normalized, synergy, joint_mi, len(data), empty_fraction


def load_masks():
    continuum = pd.read_csv(
        ROOT / "agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat",
        sep=r"\s+",
        header=None,
        names=["time", "flux", "error"],
    )
    lines = pd.read_csv(PROCESSED / "ngc5548_hb_velocity_bins.csv")
    lower = max(continuum.time.min(), lines.mjd.min())
    upper = min(continuum.time.max(), lines.mjd.max())
    grid = np.arange(lower, upper + 1.0)
    return grid, np.isin(grid, continuum.time), np.isin(grid, lines.mjd)


GRID, CONTINUUM_MASK, LINE_MASK = load_masks()


def run_realization(index: int):
    rng = np.random.default_rng(750_000 + index)
    complete = [simulate_ou(len(GRID), scale, rng) for scale in (45.0, 65.0, 55.0, 60.0)]

    interpolated = []
    for variable_index, values in enumerate(complete):
        mask = CONTINUUM_MASK if variable_index == 0 else LINE_MASK
        interpolated.append(zscore(np.interp(GRID, GRID[mask], values[mask])))

    rows = []
    for lag in LAGS:
        native_rows = LINE_MASK[lag:] & CONTINUUM_MASK[:-lag] & LINE_MASK[:-lag]
        for condition, arrays, mask in (
            ("complete_daily", complete, None),
            ("observed_only", complete, native_rows),
            ("observed_plus_interpolation", interpolated, None),
        ):
            normalized, synergy, joint_mi, n_eff, empty = surd_metrics(
                arrays[2], [arrays[0], arrays[1], arrays[3]], lag, mask
            )
            rows.append((condition, lag, normalized, synergy, joint_mi, n_eff, empty))
    return index, rows


def main():
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(run_realization, range(N_REALIZATIONS), chunksize=2))

    records = []
    peak_records = []
    for realization, rows in results:
        frame = pd.DataFrame(
            rows,
            columns=["condition", "lag", "normalized_synergy", "synergy_bits", "joint_mi", "n_eff", "empty_fraction"],
        )
        frame.insert(0, "realization", realization)
        records.append(frame)
        for condition, subset in frame.groupby("condition"):
            valid = subset.dropna(subset=["normalized_synergy"])
            norm_peak = valid.loc[valid.normalized_synergy.idxmax()]
            abs_peak = valid.loc[valid.synergy_bits.idxmax()]
            peak_records.append(
                {
                    "realization": realization,
                    "condition": condition,
                    "normalized_peak_lag": norm_peak.lag,
                    "normalized_peak": norm_peak.normalized_synergy,
                    "normalized_peak_joint_mi": norm_peak.joint_mi,
                    "absolute_peak_lag": abs_peak.lag,
                    "absolute_peak_bits": abs_peak.synergy_bits,
                    "n_realizations": N_REALIZATIONS,
                }
            )

    all_curves = pd.concat(records, ignore_index=True)
    grouped = all_curves.groupby(["condition", "lag"], as_index=False).agg(
        normalized_median=("normalized_synergy", "median"),
        normalized_p16=("normalized_synergy", lambda x: np.nanpercentile(x, 16)),
        normalized_p84=("normalized_synergy", lambda x: np.nanpercentile(x, 84)),
        synergy_bits_median=("synergy_bits", "median"),
        joint_mi_median=("joint_mi", "median"),
        n_eff_median=("n_eff", "median"),
        empty_fraction_median=("empty_fraction", "median"),
    )
    grouped["n_realizations"] = N_REALIZATIONS
    grouped.to_csv(PROCESSED / "matched_negative_control_curves.csv", index=False)
    pd.DataFrame(peak_records).to_csv(PROCESSED / "matched_negative_control_peaks.csv", index=False)
    print(pd.DataFrame(peak_records).groupby("condition").median(numeric_only=True).to_string())


if __name__ == "__main__":
    main()
