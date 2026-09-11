"""Task 3: Calibrate false positives and detection power across a parameter grid.

Injects known delayed line responses into synthetic DRW signals sampled at actual
observational dates with propagated errors, measuring:
- Detection power (true recovery within tolerance and above significance threshold)
- False positive rate under the null (A = 0)
- Lag recovery bias and root-mean-square error
across response amplitudes and delays.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .adopted_data import ROOT, load_adopted
from .information import marginal_edges
from .continuum_baseline import decompose_1d
from .benchmark_information import ou
from .sampling import tuples

HERE = Path(__file__).resolve().parent


def simulate_delayed_line(cont, lines, amplitude, delay_days, width_days=5.0,
                          timescale_days=50.0, seed=123):
    """Simulate continuum driver and delayed line response on adopted dates."""
    rng = np.random.default_rng(seed)
    burn = 300.0
    start = min(lines.jd_offset.min(), cont.jd_offset.min()) - burn
    stop = max(lines.jd_offset.max(), cont.jd_offset.max())
    step = 0.5
    grid = np.arange(start, stop + step, step)

    # Continuum driver
    driver = ou(rng, grid, timescale_days)
    sim_cont = cont.copy()
    cont_scale = float(cont.continuum.std(ddof=1))
    sim_cont['continuum'] = (cont.continuum.mean()
                             + cont_scale * np.interp(cont.jd_offset, grid, driver)
                             + rng.normal(size=len(cont)) * cont.continuum_error.to_numpy())

    # Line response
    resp_lags = np.arange(delay_days - width_days / 2.0, delay_days + width_days / 2.0 + step / 2.0, step)
    resp = np.mean([np.interp(lines.jd_offset.to_numpy() - d, grid, driver) for d in resp_lags], axis=0)

    sim_lines = lines.copy()
    for name in ['blue', 'core', 'red']:
        indep = np.interp(lines.jd_offset, grid, ou(rng, grid, timescale_days))
        # When amplitude = 0, signal is purely independent line variability
        signal = amplitude * resp + (1.0 - amplitude) * indep if amplitude > 0 else indep
        sim_lines[name] = (lines[name].mean()
                           + lines[name].std(ddof=1) * signal
                           + rng.normal(size=len(lines)) * lines[name + '_error'].to_numpy())
    sim_lines['profile_total'] = sim_lines[['blue', 'core', 'red']].sum(axis=1)
    return sim_lines, sim_cont


def run_single_grid_realization(lines, cont, sampling_config, info_config,
                                amplitude, delay_days, seed, target='core',
                                lags=range(0, 61), method='gap_limited'):
    """Evaluate lag recovery for one synthetic realization."""
    sim_lines, sim_cont = simulate_delayed_line(cont, lines, amplitude, delay_days, seed=seed)

    # 1-predictor scan: I(Target; Continuum)
    edges_y = marginal_edges(sim_lines[target], info_config['bins'])
    edges_x = marginal_edges(sim_cont.continuum, info_config['bins'])

    mi_list = []
    for lag in lags:
        frame = tuples(sim_lines, sim_cont, lag, method, sampling_config)
        if len(frame) >= info_config['minimum_tuples']:
            cols = ['target_' + target, 'continuum']
            hist, _ = np.histogramdd(frame[cols].to_numpy(), bins=[edges_y, edges_x])
            res = decompose_1d(hist)
            if res is not None:
                mi_list.append((lag, res[0]))

    if not mi_list:
        return {'peak_lag': np.nan, 'peak_mi': 0.0}

    best_lag, best_mi = max(mi_list, key=lambda x: x[1])
    return {'peak_lag': best_lag, 'peak_mi': best_mi}


def worker_grid(job):
    lines, cont, sampling_config, info_config, amp, delay, seed = job
    res = run_single_grid_realization(lines, cont, sampling_config, info_config, amp, delay, seed)
    return {
        'amplitude': amp,
        'true_delay': delay,
        'seed': seed,
        'recovered_lag': res['peak_lag'],
        'peak_mi': res['peak_mi']
    }


def run_detection_power_grid(n_realizations=25, workers=2, out_dir=None):
    sampling = json.loads((HERE / 'sampling_config.json').read_text())
    config = json.loads((HERE / 'information_config.json').read_text())

    lines, cont = load_adopted(sampling['continuum'])
    out = Path(out_dir) if out_dir else (
        ROOT / 'agn_surd_project/processed/reconciled/detection_power'
    )
    out.mkdir(parents=True, exist_ok=True)

    # Grid definition:
    # Amplitudes: 0.0 (Null control), 0.2, 0.4, 0.6, 0.8, 1.0
    # True Delays: 10, 20, 30 days
    amplitudes = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    delays = [10, 20, 30]

    print("==================================================================")
    print("Running Task 3 Detection Power and False Positive Calibration Grid")
    print(f"Amplitudes: {amplitudes} | Delays: {delays} | M = {n_realizations} reps")
    print("==================================================================")

    jobs = []
    base_seed = 300000
    seed_idx = 0

    # 1. Null control jobs (Amplitude = 0)
    for rep in range(n_realizations):
        jobs.append((lines, cont, sampling, config, 0.0, 0, base_seed + seed_idx))
        seed_idx += 1

    # 2. Injected delay jobs
    for amp in amplitudes[1:]:
        for delay in delays:
            for rep in range(n_realizations):
                jobs.append((lines, cont, sampling, config, amp, delay, base_seed + seed_idx))
                seed_idx += 1

    print(f"Total realizations to run: {len(jobs)}...")
    if workers == 1:
        results = [worker_grid(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(worker_grid, jobs))

    res_df = pd.DataFrame(results)
    res_df.to_csv(out / 'grid_raw_realizations.csv', index=False)

    # Calculate detection threshold from the null distribution (95th percentile of max MI under amp=0)
    null_runs = res_df[res_df.amplitude == 0.0]
    detection_threshold = float(np.percentile(null_runs.peak_mi, 95))
    print(f"\nEmpirical 95% Null Detection Threshold (A=0): {detection_threshold:.4f} bits")

    # Evaluate False Positive Rate (fraction of null runs exceeding threshold)
    fpr = float(np.mean(null_runs.peak_mi >= detection_threshold))

    # Calculate Detection Power & Accuracy per (amp, delay)
    summary_rows = []
    # Add null summary
    summary_rows.append({
        'amplitude': 0.0,
        'true_delay': 0,
        'detection_power': fpr,  # FPR for null
        'recovery_within_3d': np.nan,
        'mean_recovered_lag': float(null_runs.recovered_lag.mean()),
        'lag_bias': np.nan,
        'lag_rmse': np.nan,
        'mean_peak_mi': float(null_runs.peak_mi.mean()),
        'n_realizations': len(null_runs)
    })

    tolerance = 3  # days

    for amp in amplitudes[1:]:
        for delay in delays:
            sub = res_df[(res_df.amplitude == amp) & (res_df.true_delay == delay)]
            lags = sub.recovered_lag.to_numpy()
            mis = sub.peak_mi.to_numpy()

            # True detection: within tolerance of true lag AND above detection threshold
            detected = (np.abs(lags - delay) <= tolerance) & (mis >= detection_threshold)
            within_tol = (np.abs(lags - delay) <= tolerance)

            bias = float(np.mean(lags - delay))
            rmse = float(np.sqrt(np.mean((lags - delay) ** 2)))

            summary_rows.append({
                'amplitude': amp,
                'true_delay': delay,
                'detection_power': float(np.mean(detected)),
                'recovery_within_3d': float(np.mean(within_tol)),
                'mean_recovered_lag': float(np.mean(lags)),
                'lag_bias': bias,
                'lag_rmse': rmse,
                'mean_peak_mi': float(np.mean(mis)),
                'n_realizations': len(sub)
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = out / 'detection_power_grid_summary.csv'
    summary_df.to_csv(summary_path, index=False)

    print("\n==================================================================")
    print("DETECTION POWER SUMMARY TABLE:")
    print("==================================================================")
    print(summary_df[['amplitude', 'true_delay', 'detection_power', 'recovery_within_3d', 'lag_bias', 'lag_rmse', 'mean_peak_mi']].to_string(index=False))

    # Plot Detection Power curves
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    colors = {10: 'tab:blue', 20: 'tab:green', 30: 'tab:red'}

    for delay in delays:
        sub = summary_df[summary_df.true_delay == delay]
        ax1.plot(sub.amplitude, sub.detection_power, 'o-', label=f'$\\tau_{{\\mathrm{{true}}}} = {delay}$ d', color=colors[delay], lw=2)
        ax2.plot(sub.amplitude, sub.lag_rmse, 's--', label=f'$\\tau_{{\\mathrm{{true}}}} = {delay}$ d', color=colors[delay], lw=1.8)

    ax1.axhline(0.05, color='gray', ls=':', label='Nominal False Positive Rate (5%)')
    ax1.set_xlabel('Response Amplitude $A$', fontsize=12)
    ax1.set_ylabel('Detection Power $P(\\text{Detect})$', fontsize=12)
    ax1.set_title('Detection Power vs Response Amplitude (95% Significance)', fontsize=12)
    ax1.set_ylim(-0.05, 1.05)
    ax1.grid(True, ls='--', alpha=0.5)
    ax1.legend(loc='lower right', fontsize=10)

    ax2.set_xlabel('Response Amplitude $A$', fontsize=12)
    ax2.set_ylabel('Lag Recovery RMSE (days)', fontsize=12)
    ax2.set_title('Lag Recovery Root-Mean-Square Error', fontsize=12)
    ax2.grid(True, ls='--', alpha=0.5)
    ax2.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    plot_file = out / 'detection_power_curves.png'
    fig.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved plot: {plot_file}")

    manifest = {
        'status': 'Task 3 Detection Power and False Positive Calibration Grid Completed',
        'detection_threshold_bits': detection_threshold,
        'false_positive_rate': fpr,
        'tolerance_days': tolerance,
        'grid_amplitudes': amplitudes,
        'grid_delays': delays,
        'realizations_per_cell': n_realizations
    }
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"All Task 3 outputs saved to {out}")


def main():
    parser = argparse.ArgumentParser(description='Run detection power calibration grid.')
    parser.add_argument('--reps', type=int, default=25, help='Realizations per grid point (default: 25)')
    parser.add_argument('--workers', type=int, default=2, help='Number of parallel workers (default: 2)')
    parser.add_argument('--output-dir', type=str, default=None, help='Output directory')
    args = parser.parse_args()

    run_detection_power_grid(n_realizations=args.reps, workers=args.workers, out_dir=args.output_dir)


if __name__ == '__main__':
    main()
