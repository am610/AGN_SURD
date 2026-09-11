"""Run the synthetic-recovery and conditional-binning validation suites.

The outputs from this script are the source of truth for the corresponding
tables and figures in the manuscript. Run from the repository root with:

    python run_final_validations.py
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


ROOT = Path(__file__).resolve().parent
PROCESSED = ROOT / "agn_surd_project" / "processed"
OVERLEAF = ROOT / "overleaf_draft"
CONTINUUM_PATH = (
    ROOT / "agn_surd_project" / "agn_data" / "ngc5548_agnwatch" / "c5100.dat"
)
LINE_PATH = PROCESSED / "ngc5548_hb_velocity_bins.csv"

sys.path.insert(0, str(ROOT / "SURD" / "utils"))
import surd  # noqa: E402


MJD_MIN = 47512.0
MJD_MAX = 49255.0
LAGS = np.arange(1, 41)
CONDITIONAL_LAGS = np.arange(1, 121)
ATOMS = ("U1", "U2", "R12", "S12")


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / values.std(ddof=0)


def load_observing_footprint() -> dict[str, np.ndarray]:
    continuum = pd.read_csv(
        CONTINUUM_PATH,
        sep=r"\s+",
        header=None,
        names=["mjd", "flux", "err"],
    )
    continuum = continuum[
        continuum["mjd"].between(MJD_MIN, MJD_MAX)
    ].dropna()

    line = pd.read_csv(LINE_PATH)
    line = line[line["mjd"].between(MJD_MIN, MJD_MAX)].dropna()

    grid = np.arange(MJD_MIN, MJD_MAX + 1.0)
    return {
        "grid": grid,
        "continuum_mjd": continuum["mjd"].to_numpy(),
        "continuum_noise": (
            continuum["err"] / continuum["flux"].std(ddof=1)
        ).to_numpy(),
        "line_mjd": line["mjd"].to_numpy(),
        "line_noise": (
            line["core_error"] / line["core_flux"].std(ddof=1)
        ).to_numpy(),
    }


def load_real_series() -> np.ndarray:
    footprint = load_observing_footprint()
    grid = footprint["grid"]

    continuum = pd.read_csv(
        CONTINUUM_PATH,
        sep=r"\s+",
        header=None,
        names=["mjd", "flux", "err"],
    ).sort_values("mjd")
    line = pd.read_csv(LINE_PATH).sort_values("mjd")

    cont = np.interp(grid, continuum["mjd"], continuum["flux"])
    blue = np.interp(grid, line["mjd"], line["blue_wing_flux"])
    core = np.interp(grid, line["mjd"], line["core_flux"])
    return np.vstack([zscore(cont), zscore(blue), zscore(core)])


def generate_drw(
    rng: np.random.Generator,
    size: int,
    tau: float,
    sigma: float = 1.0,
) -> np.ndarray:
    values = np.empty(size, dtype=float)
    values[0] = rng.normal(0.0, sigma)
    coefficient = np.exp(-1.0 / tau)
    innovation_sigma = sigma * np.sqrt(1.0 - coefficient**2)
    for index in range(1, size):
        values[index] = (
            coefficient * values[index - 1] + rng.normal(0.0, innovation_sigma)
        )
    return values


def observe_and_interpolate(
    signal: np.ndarray,
    grid: np.ndarray,
    observed_mjd: np.ndarray,
    noise_scale: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    observed = np.interp(observed_mjd, grid, zscore(signal))
    observed += rng.normal(0.0, noise_scale)

    # Several spectra share an integer MJD. Average those measurements before
    # interpolation so every observation contributes without duplicate x values.
    sampled = pd.DataFrame({"mjd": observed_mjd, "value": observed})
    sampled = sampled.groupby("mjd", as_index=False)["value"].mean()
    return zscore(np.interp(grid, sampled["mjd"], sampled["value"]))


def surd_atoms(target: np.ndarray, first: np.ndarray, second: np.ndarray, lag: int,
               bins: int = 6) -> dict[str, float]:
    sample = np.column_stack(
        [target[lag:], first[:-lag], second[:-lag]]
    )
    histogram, _ = np.histogramdd(sample, bins=bins)
    redundancy, synergy, _, _ = surd.surd(histogram)
    return {
        "U1": float(redundancy.get((1,), 0.0)),
        "U2": float(redundancy.get((2,), 0.0)),
        "R12": float(redundancy.get((1, 2), 0.0)),
        "S12": float(synergy.get((1, 2), 0.0)),
    }


def scan_atoms(series: np.ndarray, bins: int = 6) -> dict[str, np.ndarray]:
    curves = {atom: [] for atom in ATOMS}
    for lag in LAGS:
        values = surd_atoms(series[2], series[0], series[1], int(lag), bins)
        for atom in ATOMS:
            curves[atom].append(values[atom])
    return {atom: np.asarray(values) for atom, values in curves.items()}


def synthetic_realization(seed: int) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(seed)
    footprint = load_observing_footprint()
    grid = footprint["grid"]
    margin = 60
    size = len(grid) + margin

    recovery_rows: list[dict] = []
    curve_rows: list[dict] = []

    for case in range(1, 5):
        first_full = generate_drw(rng, size, tau=50.0)
        second_full = generate_drw(rng, size, tau=30.0)
        now = slice(margin, None)

        first = first_full[now]
        independent_second = second_full[now]
        if case == 1:
            second = independent_second
            target = first_full[margin - 15 : -15]
        elif case == 2:
            second = first + rng.normal(0.0, 0.20, len(first))
            target = first_full[margin - 15 : -15]
        elif case == 3:
            second = independent_second
            first_lagged = zscore(first_full[margin - 15 : -15])
            second_lagged = zscore(second_full[margin - 15 : -15])
            target = first_lagged * second_lagged
        else:
            second = independent_second
            target = (
                0.5 * zscore(first_full[margin - 10 : -10])
                + 0.5 * zscore(second_full[margin - 20 : -20])
            )

        observed_first = observe_and_interpolate(
            first,
            grid,
            footprint["continuum_mjd"],
            footprint["continuum_noise"],
            rng,
        )
        observed_second = observe_and_interpolate(
            second,
            grid,
            footprint["line_mjd"],
            footprint["line_noise"],
            rng,
        )
        observed_target = observe_and_interpolate(
            target,
            grid,
            footprint["line_mjd"],
            footprint["line_noise"],
            rng,
        )
        curves = scan_atoms(
            np.vstack([observed_first, observed_second, observed_target]), bins=6
        )

        for lag_index, lag in enumerate(LAGS):
            row = {"seed": seed, "case": case, "lag": int(lag)}
            row.update({atom: curves[atom][lag_index] for atom in ATOMS})
            curve_rows.append(row)

        diagnostics = {
            1: [("U1", 15)],
            2: [("R12", 15)],
            3: [("S12", 15)],
            4: [("U1", 10), ("U2", 20)],
        }[case]
        for atom, true_lag in diagnostics:
            peak_index = int(np.argmax(curves[atom]))
            true_index = int(np.where(LAGS == true_lag)[0][0])
            competitors = [name for name in ATOMS if name != atom]
            recovery_rows.append(
                {
                    "seed": seed,
                    "case": case,
                    "diagnostic": atom,
                    "true_lag": true_lag,
                    "recovered_lag": int(LAGS[peak_index]),
                    "within_5d": abs(int(LAGS[peak_index]) - true_lag) <= 5,
                    "intended_atom_largest_at_true_lag": curves[atom][true_index]
                    >= max(curves[name][true_index] for name in competitors),
                    "atom_margin_at_true_lag": curves[atom][true_index]
                    - max(curves[name][true_index] for name in competitors),
                }
            )

    return recovery_rows, curve_rows


def summarize_synthetic(recovery: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["case", "diagnostic", "true_lag"]
    for key, group in recovery.groupby(keys, sort=True):
        recovered = group["recovered_lag"].to_numpy(dtype=float)
        true_lag = float(key[2])
        rows.append(
            {
                "case": int(key[0]),
                "diagnostic": key[1],
                "true_lag": true_lag,
                "n_realizations": len(group),
                "mean_recovered_lag": recovered.mean(),
                "sd_recovered_lag": recovered.std(ddof=1),
                "median_recovered_lag": np.median(recovered),
                "bias": (recovered - true_lag).mean(),
                "rmse": np.sqrt(np.mean((recovered - true_lag) ** 2)),
                "recovery_rate_within_5d": group["within_5d"].mean(),
                "intended_atom_largest_rate": group[
                    "intended_atom_largest_at_true_lag"
                ].mean(),
                "median_atom_margin": group["atom_margin_at_true_lag"].median(),
            }
        )
    return pd.DataFrame(rows)


def plot_synthetic(curves: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    titles = {
        1: "A: Single driver (15 d)",
        2: "B: Redundant proxies (15 d)",
        3: "C: Joint nonlinear drivers (15 d)",
        4: "D: Two independent drivers (10 d, 20 d)",
    }
    styles = {
        "U1": ("#1768ac", ":", r"Unique $U_1$"),
        "U2": ("#d1495b", "-.", r"Unique $U_2$"),
        "R12": ("#5b5f62", "--", r"Redundancy $R_{12}$"),
        "S12": ("#7a3e9d", "-", r"Synergy $S_{12}$"),
    }
    intended = {1: ["U1"], 2: ["R12"], 3: ["S12"], 4: ["U1", "U2"]}

    for case, axis in zip(range(1, 5), axes.flat):
        case_data = curves[curves["case"] == case]
        for atom in ATOMS:
            pivot = case_data.pivot(index="seed", columns="lag", values=atom)
            median = pivot.median(axis=0)
            p16 = pivot.quantile(0.16, axis=0)
            p84 = pivot.quantile(0.84, axis=0)
            color, linestyle, label = styles[atom]
            axis.plot(median.index, median, color=color, linestyle=linestyle,
                      linewidth=2.0, label=label)
            if atom in intended[case]:
                axis.fill_between(median.index, p16, p84, color=color, alpha=0.14)
        for true_lag in ({1: [15], 2: [15], 3: [15], 4: [10, 20]}[case]):
            axis.axvline(true_lag, color="#1b4332", linewidth=1.2, alpha=0.8)
        axis.set_title(titles[case])
        axis.set_ylabel("Information (bits)")
        axis.legend(fontsize=8)
    for axis in axes[1]:
        axis.set_xlabel("Lag (days)")
    fig.tight_layout()
    fig.savefig(OVERLEAF / "figure5_synthetic_benchmarks.png", dpi=300)
    plt.close(fig)


def histogram_edges(data: np.ndarray, bins: int, method: str) -> int | list[np.ndarray]:
    if method == "equal_width":
        return bins
    edges = []
    quantiles = np.linspace(0.0, 1.0, bins + 1)
    for column in data.T:
        column_edges = np.unique(np.quantile(column, quantiles))
        if len(column_edges) != bins + 1:
            return bins
        column_edges[0] = np.nextafter(column_edges[0], -np.inf)
        column_edges[-1] = np.nextafter(column_edges[-1], np.inf)
        edges.append(column_edges)
    return edges


def conditional_values(
    series: np.ndarray,
    lag: int,
    bins: int,
    method: str,
) -> dict[str, float]:
    target = series[2, lag:]
    first = series[0, :-lag]
    second = series[1, :-lag]
    history = series[2, :-lag]
    data = np.column_stack([target, first, second, history])
    histogram, _ = np.histogramdd(data, bins=histogram_edges(data, bins, method))

    total = histogram.sum()
    probability = histogram / total
    conditional_synergy = 0.0
    conditional_leak = 0.0
    for history_bin in range(histogram.shape[-1]):
        history_probability = probability[:, :, :, history_bin].sum()
        if history_probability <= 0.0:
            continue
        conditioned = probability[:, :, :, history_bin] / history_probability
        _, synergy, _, leak = surd.surd(conditioned.copy())
        conditional_synergy += history_probability * synergy.get((1, 2), 0.0)
        conditional_leak += history_probability * leak

    nonempty = histogram[histogram > 0]
    return {
        "conditional_synergy": float(conditional_synergy),
        "conditional_leak": float(conditional_leak),
        "empty_fraction": float(np.mean(histogram == 0)),
        "nonempty_bins": int(len(nonempty)),
        "mean_nonempty_occupancy": float(nonempty.mean()),
        "singleton_fraction_nonempty": float(np.mean(nonempty == 1)),
    }


def conditional_scan(series: np.ndarray, bins: int, method: str) -> pd.DataFrame:
    rows = []
    for lag in CONDITIONAL_LAGS:
        row = {"lag": int(lag), "bins": bins, "method": method}
        row.update(conditional_values(series, int(lag), bins, method))
        rows.append(row)
    return pd.DataFrame(rows)


_SURROGATE_SERIES: np.ndarray | None = None


def initialize_surrogate_worker(series: np.ndarray) -> None:
    global _SURROGATE_SERIES
    _SURROGATE_SERIES = series


def conditional_surrogate_max(task: tuple[str, int, int]) -> tuple[str, int, int, float]:
    method, bins, seed = task
    if _SURROGATE_SERIES is None:
        raise RuntimeError("Surrogate worker was not initialized")
    rng = np.random.default_rng(seed)
    length = _SURROGATE_SERIES.shape[1]
    shifted = np.vstack(
        [
            np.roll(row, rng.integers(120, length - 120))
            for row in _SURROGATE_SERIES
        ]
    )
    maximum = conditional_scan(shifted, bins, method)["conditional_synergy"].max()
    return method, bins, seed, float(maximum)


def plot_binning(scans: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, scans["bins"].nunique()))
    for axis, method in zip(axes, ("equal_width", "quantile")):
        method_data = scans[scans["method"] == method]
        for color, bins in zip(colors, sorted(method_data["bins"].unique())):
            subset = method_data[method_data["bins"] == bins]
            axis.plot(subset["lag"], subset["conditional_synergy"], color=color,
                      linewidth=1.6, label=f"{bins} bins")
        axis.set_title("Equal-width bins" if method == "equal_width" else "Quantile bins")
        axis.set_xlabel("Lag (days)")
        axis.set_ylabel(r"Conditioned synergy $S_{12}$ (bits)")
        axis.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(OVERLEAF / "figure9_conditional_binning_sensitivity.png", dpi=300)
    plt.close(fig)


def run_validations(realizations: int, surrogates: int, workers: int) -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    OVERLEAF.mkdir(parents=True, exist_ok=True)

    print(f"Running {realizations} synthetic realizations...")
    recovery_rows: list[dict] = []
    curve_rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for index, (recovery, curves) in enumerate(
            executor.map(synthetic_realization, range(1000, 1000 + realizations)), start=1
        ):
            recovery_rows.extend(recovery)
            curve_rows.extend(curves)
            if index % 10 == 0 or index == realizations:
                print(f"  synthetic realizations complete: {index}/{realizations}", flush=True)

    recovery = pd.DataFrame(recovery_rows)
    curves = pd.DataFrame(curve_rows)
    summary = summarize_synthetic(recovery)
    recovery.to_csv(PROCESSED / "synthetic_validation_realizations.csv", index=False)
    curves.to_csv(PROCESSED / "synthetic_validation_curves.csv", index=False)
    summary.to_csv(PROCESSED / "synthetic_validation_summary.csv", index=False)
    plot_synthetic(curves)
    print(summary.to_string(index=False))

    print("Running conditional binning scans...")
    real_series = load_real_series()
    scan_frames = []
    for method in ("equal_width", "quantile"):
        for bins in range(3, 9):
            scan = conditional_scan(real_series, bins, method)
            scan_frames.append(scan)
            peak = scan.loc[scan["conditional_synergy"].idxmax()]
            print(
                f"  {method}, {bins} bins: peak={int(peak['lag'])} d, "
                f"S={peak['conditional_synergy']:.4f}, "
                f"empty={100 * peak['empty_fraction']:.1f}%",
                flush=True,
            )
    scans = pd.concat(scan_frames, ignore_index=True)
    scans.to_csv(PROCESSED / "conditional_binning_sensitivity_curves.csv", index=False)

    print(f"Running {surrogates} auxiliary surrogates per binning setting...")
    tasks = [
        (method, bins, 7000 + method_index * 10000 + bins * 100 + seed)
        for method_index, method in enumerate(("equal_width", "quantile"))
        for bins in range(3, 9)
        for seed in range(surrogates)
    ]
    surrogate_rows = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_surrogate_worker,
        initargs=(real_series,),
    ) as executor:
        for index, result in enumerate(executor.map(conditional_surrogate_max, tasks), start=1):
            method, bins, seed, maximum = result
            surrogate_rows.append(
                {"method": method, "bins": bins, "seed": seed, "null_max": maximum}
            )
            if index % 50 == 0 or index == len(tasks):
                print(f"  surrogate scans complete: {index}/{len(tasks)}", flush=True)
    surrogate_results = pd.DataFrame(surrogate_rows)
    surrogate_results.to_csv(
        PROCESSED / "conditional_binning_sensitivity_surrogates.csv", index=False
    )

    summary_rows = []
    for (method, bins), scan in scans.groupby(["method", "bins"], sort=True):
        peak = scan.loc[scan["conditional_synergy"].idxmax()]
        null = surrogate_results[
            (surrogate_results["method"] == method)
            & (surrogate_results["bins"] == bins)
        ]["null_max"]
        reference = scans[(scans["method"] == method) & (scans["bins"] == 6)]
        correlation = np.corrcoef(
            scan.sort_values("lag")["conditional_synergy"],
            reference.sort_values("lag")["conditional_synergy"],
        )[0, 1]
        p_global = (1 + np.sum(null >= peak["conditional_synergy"])) / (1 + len(null))
        summary_rows.append(
            {
                "method": method,
                "bins": bins,
                "peak_lag": int(peak["lag"]),
                "peak_synergy": peak["conditional_synergy"],
                "empty_fraction_at_peak": peak["empty_fraction"],
                "nonempty_bins_at_peak": int(peak["nonempty_bins"]),
                "mean_nonempty_occupancy_at_peak": peak["mean_nonempty_occupancy"],
                "correlation_with_6_bin_curve": correlation,
                "surrogate_count": len(null),
                "median_null_max": null.median(),
                "global_p_value": p_global,
            }
        )
    binning_summary = pd.DataFrame(summary_rows)
    binning_summary.to_csv(
        PROCESSED / "conditional_binning_sensitivity_summary.csv", index=False
    )
    plot_binning(scans)
    print(binning_summary.to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--realizations", type=int, default=100)
    parser.add_argument("--surrogates", type=int, default=49)
    parser.add_argument(
        "--workers", type=int, default=min(8, os.cpu_count() or 1)
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run_validations(arguments.realizations, arguments.surrogates, arguments.workers)
