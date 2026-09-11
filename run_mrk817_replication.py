"""Independent-object KNN lag scan using the public Mrk 817 STORM 2 data.

The analysis keeps the native irregular observation times. For each integer
lag, continuum and line observations are paired one-to-one when their actual
time separation lies within the stated tolerance. Circular shifts preserve
each line's serial structure while breaking its alignment with the continuum.
"""

from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/surd-matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_knn_estimator_crosscheck import ksg_mi_multi


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "agn_surd_project" / "agn_data"
PROCESSED = ROOT / "agn_surd_project" / "processed"
FIGURE = ROOT / "overleaf_draft" / "figure11_mrk817_replication.png"
CONTINUUM = DATA / "mrk817_storm2_continuum-1180_clean.csv"
LINE_FILES = {
    "Lyalpha": DATA / "mrk817_storm2_lya-1235-1248_clean.csv",
    "N V": DATA / "mrk817_storm2_nv-1268-1290_clean.csv",
    "Si IV": DATA / "mrk817_storm2_siiv-1424-1460_clean.csv",
    "C IV": DATA / "mrk817_storm2_civ-1590-1638_clean.csv",
    "He II": DATA / "mrk817_storm2_heii-1680-1700_clean.csv",
}
LITERATURE_REST_LAGS = {
    "Lyalpha": 10.4,
    "N V": 15.5,
    "Si IV": 8.2,
    "C IV": 11.8,
    "He II": 9.0,
}
REDSHIFT = 0.03145
LAGS = np.arange(0, 41)
K_VALUES = (3, 5, 10)
PAIR_TOLERANCE_DAYS = 1.1


def load_series() -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    continuum = pd.read_csv(CONTINUUM)
    times = continuum["hjd"].to_numpy(dtype=float)
    flux = continuum["flux"].to_numpy(dtype=float)
    lines = {}
    for name, path in LINE_FILES.items():
        frame = pd.read_csv(path)
        line_times = frame["hjd"].to_numpy(dtype=float)
        if not np.array_equal(times, line_times):
            raise ValueError(f"Timestamp grid differs for {path}")
        lines[name] = frame["flux"].to_numpy(dtype=float)
    return times, flux, lines


def lag_pairs(times: np.ndarray, lag: float, tolerance: float) -> np.ndarray:
    """Greedily select unique source/target pairs nearest the requested lag."""
    candidates = []
    for source, time in enumerate(times):
        insertion = int(np.searchsorted(times, time + lag))
        for target in (insertion - 1, insertion):
            if 0 <= target < len(times):
                residual = abs((times[target] - time) - lag)
                if residual <= tolerance:
                    candidates.append((residual, source, target))

    selected = []
    used_sources: set[int] = set()
    used_targets: set[int] = set()
    for residual, source, target in sorted(candidates):
        if source not in used_sources and target not in used_targets:
            selected.append((source, target, times[target] - times[source]))
            used_sources.add(source)
            used_targets.add(target)
    return np.asarray(selected, dtype=float)


def build_pair_map(times: np.ndarray, tolerance: float) -> dict[int, np.ndarray]:
    pairs = {int(lag): lag_pairs(times, lag, tolerance) for lag in LAGS}
    minimum = min(len(value) for value in pairs.values())
    if minimum <= max(K_VALUES) + 2:
        raise ValueError(f"Too few pairs at some lags: minimum={minimum}")
    return pairs


def scan(
    continuum: np.ndarray,
    lines: dict[str, np.ndarray],
    pairs: dict[int, np.ndarray],
) -> pd.DataFrame:
    rows = []
    for line_name, line_flux in lines.items():
        for lag, lag_map in pairs.items():
            source = lag_map[:, 0].astype(int)
            target = lag_map[:, 1].astype(int)
            estimates = ksg_mi_multi(
                continuum[source], line_flux[target], k_values=K_VALUES
            )
            median_separation = float(np.median(lag_map[:, 2]))
            for k, estimate in estimates.items():
                rows.append(
                    {
                        "line": line_name,
                        "lag_days": lag,
                        "k": k,
                        "mi_nats": estimate,
                        "pair_count": len(lag_map),
                        "median_actual_lag_days": median_separation,
                    }
                )
    return pd.DataFrame(rows)


_TIMES: np.ndarray | None = None
_CONTINUUM_FLUX: np.ndarray | None = None
_LINE_FLUX: dict[str, np.ndarray] | None = None
_PAIR_MAP: dict[int, np.ndarray] | None = None


def initialize_worker(
    times: np.ndarray,
    continuum: np.ndarray,
    lines: dict[str, np.ndarray],
    pairs: dict[int, np.ndarray],
) -> None:
    global _TIMES, _CONTINUUM_FLUX, _LINE_FLUX, _PAIR_MAP
    _TIMES = times
    _CONTINUUM_FLUX = continuum
    _LINE_FLUX = lines
    _PAIR_MAP = pairs


def run_surrogate(seed: int) -> pd.DataFrame:
    if _CONTINUUM_FLUX is None or _LINE_FLUX is None or _PAIR_MAP is None:
        raise RuntimeError("Worker data were not initialized")
    rng = np.random.default_rng(seed)
    length = len(_CONTINUUM_FLUX)
    shifted = {
        name: np.roll(values, int(rng.integers(25, length - 25)))
        for name, values in _LINE_FLUX.items()
    }
    result = scan(_CONTINUUM_FLUX, shifted, _PAIR_MAP)
    result["seed"] = seed
    return result


def empirical_p(value: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= value)) / (1 + len(null)))


def summarize(real: pd.DataFrame, null: pd.DataFrame) -> pd.DataFrame:
    rows = []
    family_null = null.groupby(["seed", "k"])["mi_nats"].max()
    literature_rows = []
    for line, rest_lag in LITERATURE_REST_LAGS.items():
        tested_lag = int(np.rint(rest_lag * (1 + REDSHIFT)))
        literature_rows.append(
            null[(null["line"] == line) & (null["lag_days"] == tested_lag)]
        )
    literature_null = pd.concat(literature_rows).groupby(["seed", "k"])[
        "mi_nats"
    ].max()
    for (line, k), group in real.groupby(["line", "k"], sort=True):
        group = group.sort_values("lag_days")
        null_group = null[(null["line"] == line) & (null["k"] == k)]
        peak = group.loc[group["mi_nats"].idxmax()]
        line_null_max = null_group.groupby("seed")["mi_nats"].max().to_numpy()
        all_line_null_max = family_null.xs(k, level="k").to_numpy()
        literature_rest = LITERATURE_REST_LAGS[line]
        literature_observed = literature_rest * (1 + REDSHIFT)
        literature_lag = int(np.rint(literature_observed))
        at_literature = group[group["lag_days"] == literature_lag].iloc[0]
        null_at_literature = null_group[
            null_group["lag_days"] == literature_lag
        ]["mi_nats"].to_numpy()
        literature_family_max = literature_null.xs(k, level="k").to_numpy()
        rows.append(
            {
                "line": line,
                "k": k,
                "peak_lag_days": int(peak["lag_days"]),
                "peak_mi_nats": peak["mi_nats"],
                "peak_pair_count": int(peak["pair_count"]),
                "linewise_global_p_0_40d": empirical_p(
                    peak["mi_nats"], line_null_max
                ),
                "five_line_familywise_p_0_40d": empirical_p(
                    peak["mi_nats"], all_line_null_max
                ),
                "literature_rest_lag_days": literature_rest,
                "literature_observed_lag_days": literature_observed,
                "tested_literature_lag_days": literature_lag,
                "mi_at_literature_lag_nats": at_literature["mi_nats"],
                "pointwise_p_at_literature_lag": empirical_p(
                    at_literature["mi_nats"], null_at_literature
                ),
                "five_line_familywise_p_at_literature_lags": empirical_p(
                    at_literature["mi_nats"], literature_family_max
                ),
                "surrogate_count": len(line_null_max),
            }
        )
    return pd.DataFrame(rows)


def tolerance_sensitivity(
    times: np.ndarray,
    continuum: np.ndarray,
    lines: dict[str, np.ndarray],
) -> pd.DataFrame:
    frames = []
    for tolerance in (0.75, 1.1, 1.5):
        result = scan(continuum, lines, build_pair_map(times, tolerance))
        peaks = result.loc[result.groupby(["line", "k"])["mi_nats"].idxmax()].copy()
        peaks["pair_tolerance_days"] = tolerance
        frames.append(peaks)
    return pd.concat(frames, ignore_index=True)


def plot_results(real: pd.DataFrame, null: pd.DataFrame) -> None:
    order = list(LINE_FILES)
    colors = {3: "#2A6F97", 5: "#C7522A", 10: "#648767"}
    fig, axes = plt.subplots(2, 3, figsize=(10.8, 6.4), sharex=True)
    for axis, line in zip(axes.flat, order):
        for k in K_VALUES:
            subset = real[(real["line"] == line) & (real["k"] == k)]
            null_subset = null[(null["line"] == line) & (null["k"] == k)]
            envelope = null_subset.groupby("lag_days")["mi_nats"].quantile(0.95)
            axis.plot(
                subset["lag_days"], subset["mi_nats"], color=colors[k],
                lw=1.6, label=f"k={k}"
            )
            if k == 5:
                axis.plot(
                    envelope.index, envelope.values, color="#555555", lw=1.0,
                    ls="--", label="95% pointwise null"
                )
        observed_lag = LITERATURE_REST_LAGS[line] * (1 + REDSHIFT)
        axis.axvline(observed_lag, color="#9C2F2F", lw=1.1, ls=":")
        axis.set_title(line)
        axis.grid(alpha=0.18)
    axes[1, 2].axis("off")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    axes[1, 2].legend(handles, labels, loc="center", frameon=False)
    fig.supxlabel("Observed-frame lag (days)")
    fig.supylabel("KNN mutual information (nats)")
    fig.suptitle("Independent-object check: Mrk 817 STORM 2", y=0.99)
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--surrogates", type=int, default=999)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    PROCESSED.mkdir(parents=True, exist_ok=True)
    times, continuum, lines = load_series()
    pairs = build_pair_map(times, PAIR_TOLERANCE_DAYS)
    real = scan(continuum, lines, pairs)
    seeds = list(range(41001, 41001 + args.surrogates))
    workers = args.workers or min(os.cpu_count() or 1, 8)
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_worker,
        initargs=(times, continuum, lines, pairs),
    ) as executor:
        null = pd.concat(executor.map(run_surrogate, seeds), ignore_index=True)

    summary = summarize(real, null)
    sensitivity = tolerance_sensitivity(times, continuum, lines)
    real.to_csv(PROCESSED / "mrk817_knn_real_curves.csv", index=False)
    null.to_csv(PROCESSED / "mrk817_knn_surrogate_curves.csv", index=False)
    summary.to_csv(PROCESSED / "mrk817_knn_summary.csv", index=False)
    sensitivity.to_csv(
        PROCESSED / "mrk817_knn_pairing_sensitivity.csv", index=False
    )
    plot_results(real, null)
    print(summary.to_string(index=False))
    print("\nPairing sensitivity:")
    print(
        sensitivity[
            ["line", "k", "pair_tolerance_days", "lag_days", "mi_nats", "pair_count"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
