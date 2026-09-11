"""Incremental predictive information tests: conditioning on continuum driver and target history.

Tests whether candidate velocity wings add genuine incremental predictive information
beyond continuum and target autocorrelation, using both discrete quantile CMI,
continuous KSG nearest-neighbor CMI, circular-shift null controls preserving the (Y, Z)
relationship, and out-of-sample held-out season prediction.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import digamma
from sklearn.linear_model import Ridge
from sklearn.neighbors import NearestNeighbors

from .adopted_data import ROOT, load_adopted
from .information import marginal_edges
from .sampling import interpolate
from .decompose_all_atoms import extract_aligned_grid

HERE = Path(__file__).resolve().parent


# -----------------------------------------------------------------------------
# 1. Discrete Quantile CMI Estimator
# -----------------------------------------------------------------------------

def cmi_histogram(y, x, z, n_bins=3):
    """Compute discrete CMI I(Y; X | Z) from joint marginal-quantile histogram.

    Parameters
    ----------
    y : 1D array, target Y
    x : 1D array, candidate predictor X
    z : 1D or 2D array, conditioning variables Z (e.g. continuum and target past)
    n_bins : int, number of quantile bins per marginal
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    z = np.asarray(z, dtype=float)
    if z.ndim == 1:
        z = z[:, None]

    n_samples = len(y)
    if n_samples < 32:
        return None

    # Marginal quantile edges
    def edges(vals):
        interior = np.unique(np.quantile(vals, np.arange(1, n_bins) / n_bins))
        return np.r_[-np.inf, interior, np.inf]

    all_data = [y, x] + [z[:, j] for j in range(z.shape[1])]
    all_edges = [edges(col) for col in all_data]

    hist, _ = np.histogramdd(np.column_stack(all_data), bins=all_edges)
    total_cells = hist.size
    empty_fraction = float(np.mean(hist == 0))

    counts = hist.astype(float)
    if counts.sum() == 0:
        return None
    p = counts / counts.sum()

    # p(y, x, z) has shape (n_bins_y, n_bins_x, n_bins_z1, ...)
    # Summing over (y, x) gives p(z)
    z_axes = tuple(range(2, p.ndim))
    pz = p.sum(axis=(0, 1))

    # Summing over x gives p(y, z)
    pyz = p.sum(axis=1)

    # Summing over y gives p(x, z)
    pxz = p.sum(axis=0)

    # Expand pz to shape of pyz and pxz for proper broadcasting
    # In general CMI: I(Y; X | Z) = sum p(y, x, z) * log2( p(y,x,z)*p(z) / (p(y,z)*p(x,z)) )
    mask = p > 0
    ratio = np.zeros_like(p)

    # Broadcast pz, pyz, pxz across full tensor
    pz_full = np.broadcast_to(pz[None, None, ...], p.shape)
    pyz_full = np.broadcast_to(pyz[:, None, ...], p.shape)
    pxz_full = np.broadcast_to(pxz[None, :, ...], p.shape)

    denom = pyz_full * pxz_full
    valid_cells = mask & (denom > 0)
    ratio[valid_cells] = (p[valid_cells] * pz_full[valid_cells]) / denom[valid_cells]
    ratio[ratio <= 0] = 1.0

    cmi_bits = float(np.sum(p[valid_cells] * np.log2(ratio[valid_cells])))
    return {
        'cmi_bits': max(0.0, cmi_bits),
        'empty_fraction': empty_fraction,
        'total_cells': total_cells,
        'samples_per_cell': n_samples / total_cells
    }


# -----------------------------------------------------------------------------
# 2. Continuous KSG Nearest-Neighbor CMI Estimator
# -----------------------------------------------------------------------------

def radius_counts(points, radii):
    """Count neighbor points within radius r using Chebyshev metric (L_inf)."""
    tree = NearestNeighbors(metric="chebyshev").fit(points)
    # radius_neighbors returns neighbors strictly within radius (or <=)
    counts = tree.radius_neighbors(points, radius=radii, return_distance=False)
    # Exclude the point itself
    return np.array([len(c) - 1 for c in counts], dtype=float)


def cmi_ksg(y, x, z, k=5):
    """Continuous KSG conditional mutual information estimator I(Y; X | Z)."""
    y = np.asarray(y, dtype=float)[:, None] if y.ndim == 1 else y
    x = np.asarray(x, dtype=float)[:, None] if x.ndim == 1 else x
    z = np.asarray(z, dtype=float)[:, None] if z.ndim == 1 else z

    # The Chebyshev metric is scale sensitive at finite sample size.  Standardize
    # every marginal before the neighbour search so flux units do not determine
    # which dimension sets the joint radius.
    def standardize(values):
        center = np.mean(values, axis=0)
        spread = np.std(values, axis=0, ddof=1)
        spread = np.where(np.isfinite(spread) & (spread > 0), spread, 1.0)
        return (values - center) / spread

    y = standardize(y)
    x = standardize(x)
    z = standardize(z)

    # Slight deterministic jitter to break exact interpolation ties
    n_samples = len(y)
    scale = np.maximum(np.std(x, axis=0), 1.0)
    idx = np.arange(n_samples, dtype=float)[:, None]
    x_jit = x + 1e-10 * scale * np.sin(idx * 1.618 + 1.0)
    y_jit = y + 1e-10 * np.maximum(np.std(y, axis=0), 1.0) * np.sin(idx * 2.718 + 2.0)
    z_jit = z + 1e-10 * np.maximum(np.std(z, axis=0), 1.0) * np.sin(idx * 3.141 + 3.0)

    joint = np.column_stack([y_jit, x_jit, z_jit])
    yz = np.column_stack([y_jit, z_jit])
    xz = np.column_stack([x_jit, z_jit])

    # Find Chebyshev distance to k-th neighbor in the joint (y, x, z) space
    distances = NearestNeighbors(metric="chebyshev", n_neighbors=k + 1).fit(joint).kneighbors(return_distance=True)[0]
    radii = distances[:, k]  # k-th neighbor distance

    nyz = radius_counts(yz, radii)
    nxz = radius_counts(xz, radii)
    nz = radius_counts(z_jit, radii)

    # KSG CMI formula: I(Y; X | Z) = psi(k) + < psi(N_z + 1) - psi(N_xz + 1) - psi(N_yz + 1) >
    val_nats = digamma(k) + np.mean(digamma(nz + 1) - digamma(nxz + 1) - digamma(nyz + 1))
    val_bits = float(val_nats / np.log(2.0))
    return max(0.0, val_bits)


# -----------------------------------------------------------------------------
# 3. Incremental Information Scanner
# -----------------------------------------------------------------------------

def scan_incremental_pair(grid, target='core', candidate='blue', conditions=['continuum', 'core'],
                          lags=range(1, 121), n_bins=3, k_values=[3, 5]):
    """Scan candidate predictor X against target Y conditioning on Z over lags."""
    t_vals = grid['time'].to_numpy()
    rows = []

    for lag in lags:
        # Target Y at future time t
        y = grid[target].to_numpy()

        # Candidate predictor X at past time t - lag
        query = t_vals - lag
        x = grid[candidate].to_numpy()  # aligned series at t - lag requires shifting or indexing
        # Because grid has uniform 1.0-day steps, query = t - lag is simply indexing with shift
        if lag >= len(grid):
            continue

        # Target future slice: grid[target] from index `lag` onward
        y_fut = grid[target].iloc[lag:].to_numpy()
        # Candidate past slice: grid[candidate] up to len - lag
        x_past = grid[candidate].iloc[:-lag].to_numpy()
        # Condition past slices: grid[c] up to len - lag
        z_past = np.column_stack([grid[c].iloc[:-lag].to_numpy() for c in conditions])

        # Mask NaNs (e.g. from gap limits)
        valid = np.isfinite(y_fut) & np.isfinite(x_past)
        for j in range(z_past.shape[1]):
            valid &= np.isfinite(z_past[:, j])

        n_tuples = int(np.sum(valid))
        if n_tuples < 32:
            continue

        y_v = y_fut[valid]
        x_v = x_past[valid]
        z_v = z_past[valid]

        # Discrete Quantile CMI
        hist_res = cmi_histogram(y_v, x_v, z_v, n_bins=n_bins)
        cmi_hist = hist_res['cmi_bits'] if hist_res else np.nan
        empty_frac = hist_res['empty_fraction'] if hist_res else np.nan

        # Continuous KSG CMI
        cmi_k3 = cmi_ksg(y_v, x_v, z_v, k=3)
        cmi_k5 = cmi_ksg(y_v, x_v, z_v, k=5)

        rows.append({
            'target': target,
            'candidate': candidate,
            'conditions': ';'.join(conditions),
            'lag_days': lag,
            'tuple_count': n_tuples,
            'cmi_hist_3bin_bits': cmi_hist,
            'hist_empty_fraction': empty_frac,
            'cmi_ksg_k3_bits': cmi_k3,
            'cmi_ksg_k5_bits': cmi_k5,
        })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# 4. Circular-Shift Null Calibration for Incremental Information
# -----------------------------------------------------------------------------

def run_null_incremental_worker(job):
    """Run one null realization: circular shift candidate X relative to (Y, Z)."""
    grid, target, candidate, conditions, lag_list, n_bins, shift, seed = job
    shifted_grid = grid.copy()
    # Circularly roll candidate column
    vals = shifted_grid[candidate].to_numpy()
    valid_mask = np.isfinite(vals)
    valid_vals = vals[valid_mask]
    rolled = np.roll(valid_vals, shift)
    new_vals = vals.copy()
    new_vals[valid_mask] = rolled
    shifted_grid[candidate] = new_vals

    df = scan_incremental_pair(shifted_grid, target=target, candidate=candidate,
                               conditions=conditions, lags=lag_list, n_bins=n_bins, k_values=[5])
    if df.empty:
        return {'max_cmi_hist': 0.0, 'max_cmi_k5': 0.0}
    return {
        'max_cmi_hist': float(df['cmi_hist_3bin_bits'].max()),
        'max_cmi_k5': float(df['cmi_ksg_k5_bits'].max())
    }


def calibrate_incremental_null(grid, target='core', candidate='blue', conditions=['continuum', 'core'],
                               lags=range(1, 61), n_surrogates=50, workers=2):
    min_shift = 30
    max_shift = len(grid) - 30
    shifts = [int(s) for s in np.linspace(min_shift, max_shift, n_surrogates)]
    seeds = [200000 + i for i in range(n_surrogates)]

    jobs = [(grid, target, candidate, conditions, list(lags), 3, shifts[i], seeds[i])
            for i in range(n_surrogates)]

    if workers == 1:
        results = [run_null_incremental_worker(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(run_null_incremental_worker, jobs))

    null_hist = np.array([r['max_cmi_hist'] for r in results])
    null_k5 = np.array([r['max_cmi_k5'] for r in results])
    return null_hist, null_k5


# -----------------------------------------------------------------------------
# 5. Out-of-Sample Held-Out Season Prediction Test
# -----------------------------------------------------------------------------

def observing_season_labels(times, observed_times, boundary_gap_days=60):
    """Assign times to seasons defined only by gaps in native observations."""
    times = np.asarray(times, dtype=float)
    observed_times = np.asarray(observed_times, dtype=float)
    breaks = np.where(np.diff(observed_times) > boundary_gap_days)[0]
    edges = (observed_times[breaks] + observed_times[breaks + 1]) / 2.0
    labels = np.searchsorted(edges, times).astype(int)
    labels[(times < observed_times[0]) | (times > observed_times[-1])] = -1
    return labels


def evaluate_heldout_seasons(grid, observed_times, boundary_gap_days=60,
                             target='core', candidate='blue',
                             conditions=['continuum', 'core'], lag=15):
    """Test out-of-sample predictive improvement by holding out each observing season."""
    t_vals = grid['time'].to_numpy()
    season_labels = observing_season_labels(
        t_vals, observed_times, boundary_gap_days=boundary_gap_days
    )

    # Align data at lag
    y_fut = grid[target].iloc[lag:].to_numpy()
    x_cand = grid[candidate].iloc[:-lag].to_numpy()
    z_base = np.column_stack([grid[c].iloc[:-lag].to_numpy() for c in conditions])
    target_seasons = season_labels[lag:]
    predictor_seasons = season_labels[:-lag]

    valid = np.isfinite(y_fut) & np.isfinite(x_cand)
    for j in range(z_base.shape[1]):
        valid &= np.isfinite(z_base[:, j])
    valid &= (target_seasons >= 0) & (target_seasons == predictor_seasons)

    y_fut = y_fut[valid]
    x_cand = x_cand[valid]
    z_base = z_base[valid]
    seasons = target_seasons[valid]

    unique_seasons = np.unique(seasons)
    results = []

    for s_test in unique_seasons:
        train_idx = (seasons != s_test)
        test_idx = (seasons == s_test)

        if np.sum(test_idx) < 10 or np.sum(train_idx) < 30:
            continue

        # Baseline: Predict Y from Z alone
        model_base = Ridge(alpha=1.0)
        model_base.fit(z_base[train_idx], y_fut[train_idx])
        pred_base = model_base.predict(z_base[test_idx])
        mse_base = float(np.mean((y_fut[test_idx] - pred_base) ** 2))

        # Augmented: Predict Y from Z AND candidate X
        x_aug_train = np.column_stack([z_base[train_idx], x_cand[train_idx]])
        x_aug_test = np.column_stack([z_base[test_idx], x_cand[test_idx]])

        model_aug = Ridge(alpha=1.0)
        model_aug.fit(x_aug_train, y_fut[train_idx])
        pred_aug = model_aug.predict(x_aug_test)
        mse_aug = float(np.mean((y_fut[test_idx] - pred_aug) ** 2))

        # Variance of test set
        var_y = float(np.var(y_fut[test_idx], ddof=1))
        r2_base = 1.0 - (mse_base / var_y) if var_y > 0 else np.nan
        r2_aug = 1.0 - (mse_aug / var_y) if var_y > 0 else np.nan

        results.append({
            'heldout_season': int(s_test),
            'test_samples': int(np.sum(test_idx)),
            'mse_baseline': mse_base,
            'mse_augmented': mse_aug,
            'delta_mse': mse_aug - mse_base,  # negative means augmented improved prediction
            'r2_baseline': r2_base,
            'r2_augmented': r2_aug,
            'delta_r2': r2_aug - r2_base      # positive means augmented improved prediction
        })

    return pd.DataFrame(results)


def summarize_prediction_folds(pred_df, rng_seed=24017, n_bootstrap=5000):
    """Summarize season results without treating interpolated days as peers."""
    if pred_df.empty:
        return {
            'n_heldout_seasons': 0,
            'mean_delta_r2': np.nan,
            'weighted_delta_r2': np.nan,
            'weighted_delta_r2_ci_low': np.nan,
            'weighted_delta_r2_ci_high': np.nan,
            'improved_season_fraction': np.nan,
        }
    values = pred_df['delta_r2'].to_numpy(dtype=float)
    weights = pred_df['test_samples'].to_numpy(dtype=float)
    weighted = float(np.average(values, weights=weights))
    rng = np.random.default_rng(rng_seed)
    draws = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        idx = rng.integers(0, len(values), len(values))
        draws[i] = np.average(values[idx], weights=weights[idx])
    return {
        'n_heldout_seasons': int(len(values)),
        'mean_delta_r2': float(np.mean(values)),
        'weighted_delta_r2': weighted,
        'weighted_delta_r2_ci_low': float(np.percentile(draws, 2.5)),
        'weighted_delta_r2_ci_high': float(np.percentile(draws, 97.5)),
        'improved_season_fraction': float(np.mean(values > 0)),
    }


# -----------------------------------------------------------------------------
# 6. Main Orchestrator
# -----------------------------------------------------------------------------

def run_all_incremental_tests(out_dir=None, n_surrogates=50, workers=2):
    sampling = json.loads((HERE / 'sampling_config.json').read_text())
    config = json.loads((HERE / 'information_config.json').read_text())

    lines, cont = load_adopted(sampling['continuum'])
    out = Path(out_dir) if out_dir else (
        ROOT / 'agn_surd_project/processed/reconciled/incremental_information'
    )
    out.mkdir(parents=True, exist_ok=True)

    grid = extract_aligned_grid(lines, cont, sampling, gap_limited=True)

    # Test cases:
    # 1. Target = Core, Candidate = Blue Wing, Condition = [Continuum, Core]
    # 2. Target = Core, Candidate = Red Wing, Condition = [Continuum, Core]
    # 3. Target = Red Wing, Candidate = Blue Wing, Condition = [Continuum, Red] (Outflow test: does Blue lead Red?)
    # 4. Target = Blue Wing, Candidate = Red Wing, Condition = [Continuum, Blue] (Inflow test: does Red lead Blue?)
    test_cases = [
        ('core', 'blue', ['continuum', 'core']),
        ('core', 'red', ['continuum', 'core']),
        ('red', 'blue', ['continuum', 'red']),
        ('blue', 'red', ['continuum', 'blue'])
    ]

    all_scans = []
    summary_records = []

    print("==================================================================")
    print("Running Incremental Information Scans and Null Calibrations...")
    print("==================================================================")

    for target, cand, conds in test_cases:
        case_name = f"{target}_pred_by_{cand}_cond_on_{'_'.join(conds)}"
        print(f"\n--- Case: Target={target}, Candidate={cand}, Cond={conds} ---")

        # 1. Real Scan
        df_scan = scan_incremental_pair(grid, target=target, candidate=cand,
                                        conditions=conds, lags=range(1, 101), n_bins=3)
        df_scan.to_csv(out / f'{case_name}_scan.csv', index=False)
        all_scans.append(df_scan)

        # Max observed values
        obs_max_hist = float(df_scan['cmi_hist_3bin_bits'].max())
        obs_peak_hist_lag = int(df_scan.loc[df_scan['cmi_hist_3bin_bits'].idxmax(), 'lag_days'])
        obs_max_k5 = float(df_scan['cmi_ksg_k5_bits'].max())
        obs_peak_k5_lag = int(df_scan.loc[df_scan['cmi_ksg_k5_bits'].idxmax(), 'lag_days'])

        # 2. Null Calibration
        print(f"Running {n_surrogates} circular-shift null surrogates (preserving target + continuum)...")
        null_hist, null_k5 = calibrate_incremental_null(grid, target=target, candidate=cand,
                                                        conditions=conds, lags=range(1, 101),
                                                        n_surrogates=n_surrogates, workers=workers)

        p_hist = float((1 + np.sum(null_hist >= obs_max_hist)) / (len(null_hist) + 1))
        p_k5 = float((1 + np.sum(null_k5 >= obs_max_k5)) / (len(null_k5) + 1))

        # 3. Held-out Season Prediction
        prediction_lag = 15
        pred_df = evaluate_heldout_seasons(
            grid,
            lines.jd_offset.to_numpy(),
            boundary_gap_days=sampling['season_boundary_gap_days'],
            target=target,
            candidate=cand,
            conditions=conds,
            lag=prediction_lag,
        )
        pred_df.to_csv(out / f'{case_name}_heldout_prediction.csv', index=False)
        pred_summary = summarize_prediction_folds(pred_df)

        summary_records.append({
            'target': target,
            'candidate_wing': cand,
            'conditioning_vars': ';'.join(conds),
            'obs_peak_hist_lag': obs_peak_hist_lag,
            'obs_max_cmi_hist_bits': obs_max_hist,
            'null_hist_95pct': float(np.percentile(null_hist, 95)),
            'global_p_val_hist': p_hist,
            'obs_peak_ksg_k5_lag': obs_peak_k5_lag,
            'obs_max_cmi_ksg_k5_bits': obs_max_k5,
            'null_ksg_k5_95pct': float(np.percentile(null_k5, 95)),
            'global_p_val_ksg_k5': p_k5,
            'prediction_lag_days': prediction_lag,
            'n_heldout_seasons': pred_summary['n_heldout_seasons'],
            'mean_out_of_sample_delta_r2': pred_summary['mean_delta_r2'],
            'weighted_out_of_sample_delta_r2': pred_summary['weighted_delta_r2'],
            'weighted_delta_r2_ci_low': pred_summary['weighted_delta_r2_ci_low'],
            'weighted_delta_r2_ci_high': pred_summary['weighted_delta_r2_ci_high'],
            'improved_season_fraction': pred_summary['improved_season_fraction'],
        })

    summary_table = pd.DataFrame(summary_records)
    summary_path = out / 'incremental_information_summary.csv'
    summary_table.to_csv(summary_path, index=False)

    print("\n==================================================================")
    print("INCREMENTAL INFORMATION SUMMARY:")
    print("==================================================================")
    print(summary_table[['target', 'candidate_wing', 'obs_max_cmi_hist_bits', 'global_p_val_hist',
                         'obs_max_cmi_ksg_k5_bits', 'global_p_val_ksg_k5', 'mean_out_of_sample_delta_r2']].to_string(index=False))

    # Plot 4-panel comparison figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    for idx, (target, cand, conds) in enumerate(test_cases):
        ax = axes[idx // 2, idx % 2]
        df_scan = all_scans[idx]
        ax.plot(df_scan.lag_days, df_scan.cmi_hist_3bin_bits, label='Discrete CMI (3-bin quantile)', color='tab:blue', lw=1.8)
        ax.plot(df_scan.lag_days, df_scan.cmi_ksg_k5_bits, label='Continuous KSG CMI (k=5)', color='tab:red', lw=1.8)
        ax.axhline(0, color='gray', ls='--', alpha=0.5)
        ax.set_title(f"Target: {target.title()} | Candidate: {cand.title()} Wing\nConditioned on: Continuum + {target.title()}(t)", fontsize=11)
        ax.set_ylabel('Incremental CMI (bits)', fontsize=11)
        ax.grid(True, ls='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=9)

    axes[1, 0].set_xlabel('Lag $\\tau$ (days)', fontsize=12)
    axes[1, 1].set_xlabel('Lag $\\tau$ (days)', fontsize=12)
    fig.suptitle('Incremental Predictive Information $I(Y(t+\\tau); X_{\\mathrm{wing}}(t) \\mid F_{5100}(t), Y(t))$', fontsize=14)
    plt.tight_layout()
    plot_file = out / 'incremental_information_curves.png'
    fig.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved plot: {plot_file}")

    manifest = {
        'status': 'Task 4 Incremental Predictive Information Test Completed',
        'n_surrogates': n_surrogates,
        'estimator_comparison': ['3-bin quantile histogram', 'continuous KSG (k=5)'],
        'heldout_prediction': 'Native-gap-defined observing-season cross-validation with no cross-season tuples',
        'results': summary_records
    }
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"All Task 4 outputs saved to {out}")


def main():
    parser = argparse.ArgumentParser(description='Run incremental predictive information tests.')
    parser.add_argument('--surrogates', type=int, default=50, help='Number of null surrogates (default: 50)')
    parser.add_argument('--workers', type=int, default=2, help='Number of parallel workers (default: 2)')
    parser.add_argument('--output-dir', type=str, default=None, help='Output directory')
    args = parser.parse_args()

    run_all_incremental_tests(out_dir=args.output_dir, n_surrogates=args.surrogates, workers=args.workers)


if __name__ == '__main__':
    main()
