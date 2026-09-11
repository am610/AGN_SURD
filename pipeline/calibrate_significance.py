"""Calibrate statistical significance of observed information and lag peaks against declared null distributions."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd

from .adopted_data import ROOT, load_adopted
from .information import scan, marginal_edges
from .continuum_baseline import scan_continuum, decompose_1d, peaks_continuum
from .sampling import tuples
from .benchmark_information import simulate

HERE = Path(__file__).resolve().parent


def circular_shift_series(values, shift):
    """Circularly roll an array by shift positions."""
    return np.roll(values, shift)


def run_null_realization_continuum(lines, cont, sampling, config, method, shift, seed):
    """Run a circular-shift null realization for continuum-only baseline."""
    rng = np.random.default_rng(seed)
    shifted_cont = cont.copy()
    # Apply circular shift to continuum fluxes
    shifted_cont['continuum'] = circular_shift_series(cont['continuum'].to_numpy(), shift)

    # Scan over lags
    targets = config['targets']
    edges = {name: marginal_edges(lines[name], config['bins']) for name in targets}
    edges['continuum'] = marginal_edges(shifted_cont['continuum'], config['bins'])

    rows = []
    for lag in range(sampling['lag_min_days'], sampling['lag_max_days'] + 1, sampling['lag_step_days']):
        frame = tuples(lines, shifted_cont, lag, method, sampling)
        for target in targets:
            if len(frame) >= config['minimum_tuples']:
                cols = ['target_' + target, 'continuum']
                hist, _ = np.histogramdd(frame[cols].to_numpy(), bins=[edges[target], edges['continuum']])
                res = decompose_1d(hist)
                if res is not None:
                    mi, leak, _ = res
                    rows.append(dict(target=target, lag_days=lag, mi_bits=mi, leak=leak))

    df = pd.DataFrame(rows)
    max_stats = {}
    curves = {}
    for target in targets:
        sub = df[df.target == target]
        if not sub.empty:
            max_stats[target] = {
                'max_mi': float(sub['mi_bits'].max()),
                'min_leak': float(sub['leak'].min())
            }
            curves[target] = sub.set_index('lag_days')[['mi_bits', 'leak']].to_dict()
    return max_stats, curves


def run_worker_continuum(job):
    lines, cont, sampling, config, method, shift, seed = job
    return run_null_realization_continuum(lines, cont, sampling, config, method, shift, seed)


def benjamini_hochberg(values):
    """Return Benjamini Hochberg adjusted values in original order."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def calibrate_continuum_baseline(method='gap_limited', n_surrogates=100, workers=2, out_dir=None):
    sampling = json.loads((HERE / 'sampling_config.json').read_text())
    config = json.loads((HERE / 'information_config.json').read_text())

    lines, cont = load_adopted(sampling['continuum'])
    out = Path(out_dir) if out_dir else (
        ROOT / f'agn_surd_project/processed/reconciled/calibration/continuum_{method}'
    )
    out.mkdir(parents=True, exist_ok=True)

    print(f"=== Calibrating Continuum Baseline Significance ({method}, N={n_surrogates}) ===")

    # 1. Observed scan
    print("Computing observed scan...")
    obs_curves = scan_continuum(lines, cont, sampling, config)
    obs_sub = obs_curves[(obs_curves.method == method) & (obs_curves.status == 'ok')]
    obs_peaks = peaks_continuum(obs_sub)
    print("Observed peaks:")
    print(obs_peaks[['target', 'peak_mi_lag_days', 'peak_mi_value', 'min_leak_lag_days', 'min_leak_value']].to_string(index=False))

    # 2. Surrogate jobs
    # Use deterministic shifts between 20 and len(cont) - 20 to avoid near-zero trivial alignment
    min_shift = 20
    max_shift = len(cont) - 20
    shifts = [int(s) for s in np.linspace(min_shift, max_shift, n_surrogates)]
    seeds = [100000 + i for i in range(n_surrogates)]

    jobs = [(lines, cont, sampling, config, method, shifts[i], seeds[i]) for i in range(n_surrogates)]

    print(f"Running {n_surrogates} circular-shift null realizations with {workers} workers...")
    if workers == 1:
        results = [run_worker_continuum(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(run_worker_continuum, jobs))

    # 3. Analyze null distribution
    summary = []
    pointwise_records = []
    targets = config['targets']

    for target in targets:
        obs_row = obs_peaks[obs_peaks.target == target].iloc[0]
        obs_max_mi = float(obs_row.peak_mi_value)
        obs_min_leak = float(obs_row.min_leak_value)

        null_max_mis = np.array([r[0][target]['max_mi'] for r in results if target in r[0]])
        null_min_leaks = np.array([r[0][target]['min_leak'] for r in results if target in r[0]])

        # Global p-value for max MI (greater or equal)
        p_global_mi = float((1 + np.sum(null_max_mis >= obs_max_mi)) / (len(null_max_mis) + 1))
        # Global p-value for min leak (less or equal)
        p_global_leak = float((1 + np.sum(null_min_leaks <= obs_min_leak)) / (len(null_min_leaks) + 1))

        summary.append({
            'method': method,
            'target': target,
            'observed_peak_lag': int(obs_row.peak_mi_lag_days),
            'observed_max_mi_bits': obs_max_mi,
            'null_max_mi_median': float(np.median(null_max_mis)),
            'null_max_mi_95pct': float(np.percentile(null_max_mis, 95)),
            'global_p_value_mi': p_global_mi,
            'observed_min_leak_lag': int(obs_row.min_leak_lag_days),
            'observed_min_leak': obs_min_leak,
            'null_min_leak_median': float(np.median(null_min_leaks)),
            'null_min_leak_5pct': float(np.percentile(null_min_leaks, 5)),
            'global_p_value_leak': p_global_leak,
            'n_surrogates': len(null_max_mis)
        })

    summary_df = pd.DataFrame(summary)
    summary_df['target_family_bh_q_value_mi'] = benjamini_hochberg(
        summary_df['global_p_value_mi'].to_numpy()
    )
    summary_df['target_family_bh_q_value_leak'] = benjamini_hochberg(
        summary_df['global_p_value_leak'].to_numpy()
    )
    print("\nCalibration Summary:")
    print(summary_df[['target', 'observed_peak_lag', 'observed_max_mi_bits', 'null_max_mi_95pct', 'global_p_value_mi', 'global_p_value_leak']].to_string(index=False))

    summary_df.to_csv(out / 'significance_summary.csv', index=False)
    manifest = {
        'status': 'Continuum baseline significance calibration',
        'sampling_method': method,
        'n_surrogates': n_surrogates,
        'null_type': 'circular_shift_continuum',
        'scan_family': 'Each p value controls all scanned lags for one target; BH q values control the three target family for the selected sampling method',
        'results': summary_df.to_dict(orient='records')
    }
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"Outputs saved to {out}")
    return summary_df


def main():
    parser = argparse.ArgumentParser(description='Calibrate information significance against null distributions.')
    parser.add_argument('--method', type=str, default='gap_limited',
                        choices=['daily', 'gap_limited', 'within_season', 'native'],
                        help='Sampling method to calibrate (default: gap_limited)')
    parser.add_argument('--surrogates', type=int, default=100,
                        help='Number of surrogate realizations (default: 100)')
    parser.add_argument('--workers', type=int, default=2,
                        help='Number of parallel workers (default: 2)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory')
    args = parser.parse_args()

    calibrate_continuum_baseline(method=args.method, n_surrogates=args.surrogates,
                                 workers=args.workers, out_dir=args.output_dir)


if __name__ == '__main__':
    main()
