"""Final closure checks for estimator dependence and blocked validation.

This module records two checks that are intentionally separate from the primary
significance analysis.  The first applies temporal neighbour exclusion to the
continuous conditional information estimator.  The second replaces the two
season diagnostic with four contiguous native date blocks.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import digamma
from sklearn.linear_model import Ridge
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors

from .adopted_data import ROOT, load_adopted
from .decompose_all_atoms import extract_aligned_grid
from .incremental_information import cmi_histogram, cmi_ksg

HERE = Path(__file__).resolve().parent


def cmi_ksg_temporal(y, x, z, times, k=5, exclusion_days=0.0):
    """Estimate KSG CMI while excluding temporally adjacent neighbours."""
    y = np.asarray(y, dtype=float)[:, None]
    x = np.asarray(x, dtype=float)[:, None]
    z = np.asarray(z, dtype=float)
    if z.ndim == 1:
        z = z[:, None]
    times = np.asarray(times, dtype=float)

    def standardize(values):
        spread = np.std(values, axis=0, ddof=1)
        spread = np.where(np.isfinite(spread) & (spread > 0), spread, 1.0)
        return (values - np.mean(values, axis=0)) / spread

    y = standardize(y)
    x = standardize(x)
    z = standardize(z)
    idx = np.arange(len(y), dtype=float)[:, None]
    y = y + 1e-10 * np.sin(idx * 2.718 + 2.0)
    x = x + 1e-10 * np.sin(idx * 1.618 + 1.0)
    z = z + 1e-10 * np.sin(idx * 3.141 + 3.0)

    joint = np.column_stack([y, x, z])
    yz = np.column_stack([y, z])
    xz = np.column_stack([x, z])
    joint_distances = pairwise_distances(joint, metric="chebyshev")
    temporal = np.abs(times[:, None] - times[None, :]) <= exclusion_days

    eligible = ~temporal
    np.fill_diagonal(eligible, False)
    radii = np.empty(len(joint), dtype=float)
    for i in range(len(joint)):
        candidates = joint_distances[i, eligible[i]]
        radii[i] = np.partition(candidates, k - 1)[k - 1]

    def counts(points):
        distances = pairwise_distances(points, metric="chebyshev")
        keep = distances <= radii[:, None]
        keep &= eligible
        return keep.sum(axis=1, dtype=float)

    nyz = counts(yz)
    nxz = counts(xz)
    nz = counts(z)
    value = digamma(k) + np.mean(digamma(nz + 1) - digamma(nxz + 1) - digamma(nyz + 1))
    return float(max(0.0, value / np.log(2.0)))


def run_temporal_exclusion(grid, target, candidate, conditions, lags, exclusions):
    rows = []
    for lag in lags:
        y = grid[target].iloc[lag:].to_numpy()
        x = grid[candidate].iloc[:-lag].to_numpy()
        z = np.column_stack([grid[c].iloc[:-lag].to_numpy() for c in conditions])
        times = grid.time.iloc[:-lag].to_numpy()
        valid = np.isfinite(y) & np.isfinite(x) & np.isfinite(z).all(axis=1)
        if valid.sum() < 32:
            continue
        for exclusion in exclusions:
            rows.append({
                "target": target,
                "candidate": candidate,
                "lag_days": int(lag),
                "exclusion_days": float(exclusion),
                "tuple_count": int(valid.sum()),
                "cmi_ksg_k5_bits": cmi_ksg_temporal(y[valid], x[valid], z[valid], times[valid], k=5, exclusion_days=exclusion),
            })
    return pd.DataFrame(rows)


def blocked_prediction(grid, observed_times, target, candidate, conditions, lag=15, blocks=4):
    """Evaluate fixed contiguous native date blocks with one held out at a time."""
    labels = np.full(len(grid), -1, dtype=int)
    unique = np.unique(np.asarray(observed_times, dtype=float))
    edges = np.linspace(0, len(unique), blocks + 1, dtype=int)
    for block in range(blocks):
        dates = unique[edges[block]:edges[block + 1]]
        labels[np.isin(grid.time.to_numpy(), dates)] = block
    y = grid[target].iloc[lag:].to_numpy()
    x = grid[candidate].iloc[:-lag].to_numpy()
    z = np.column_stack([grid[c].iloc[:-lag].to_numpy() for c in conditions])
    target_block = labels[lag:]
    source_block = labels[:-lag]
    valid = np.isfinite(y) & np.isfinite(x) & np.isfinite(z).all(axis=1)
    valid &= (target_block >= 0) & (target_block == source_block)
    y, x, z, fold = y[valid], x[valid], z[valid], target_block[valid]
    rows = []
    for heldout in range(blocks):
        test = fold == heldout
        train = fold != heldout
        if test.sum() < 10 or train.sum() < 30:
            continue
        base = Ridge(alpha=1.0).fit(z[train], y[train])
        aug = Ridge(alpha=1.0).fit(np.column_stack([z[train], x[train]]), y[train])
        mse_base = float(np.mean((y[test] - base.predict(z[test])) ** 2))
        mse_aug = float(np.mean((y[test] - aug.predict(np.column_stack([z[test], x[test]]))) ** 2))
        variance = float(np.var(y[test], ddof=1))
        rows.append({
            "heldout_block": heldout,
            "test_samples": int(test.sum()),
            "mse_baseline": mse_base,
            "mse_augmented": mse_aug,
            "delta_r2": float((mse_base - mse_aug) / variance) if variance > 0 else np.nan,
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=str(ROOT / "agn_surd_project/processed/reconciled/final_closure"))
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    sampling = json.loads((HERE / "sampling_config.json").read_text())
    config = json.loads((HERE / "information_config.json").read_text())
    lines, cont = load_adopted(sampling["continuum"])
    grid = extract_aligned_grid(lines, cont, sampling, gap_limited=True)
    exclusion = run_temporal_exclusion(grid, "core", "blue", ["continuum", "core"], [1, 5, 15, 30, 60], [0, 1, 3, 5])
    exclusion.to_csv(out / "temporal_neighbor_exclusion.csv", index=False)
    blocked = blocked_prediction(grid, lines.jd_offset.to_numpy(), "core", "blue", ["continuum", "core"], lag=15, blocks=4)
    blocked.to_csv(out / "blocked_prediction.csv", index=False)
    atom = pd.read_csv(ROOT / "agn_surd_project/processed/synthetic_validation_summary.csv")
    atom.to_csv(out / "intended_atom_identification.csv", index=False)
    manifest = {
        "status": "Final closure checks completed",
        "temporal_exclusion_days": [0, 1, 3, 5],
        "blocked_validation_blocks": 4,
        "atom_source": "synthetic_validation_summary.csv",
        "notes": "These are estimator and validation sensitivity checks. They do not turn observational information into causal evidence.",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(exclusion.groupby("exclusion_days")["cmi_ksg_k5_bits"].agg(["max", "mean"]).to_string())
    print(blocked.to_string(index=False))


if __name__ == "__main__":
    main()
