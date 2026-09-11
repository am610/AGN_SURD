"""Continuum-only lag recovery baseline for adopted observed data and simulations."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from .adopted_data import ROOT, load_adopted
from .information import marginal_edges
from .sampling import tuples
from .benchmark_information import simulate

HERE = Path(__file__).resolve().parent


def decompose_1d(histogram):
    counts = np.asarray(histogram, dtype=float)
    if counts.sum() <= 0 or not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError('Invalid histogram')
    prob = counts / counts.sum()
    py = prob.sum(axis=1)
    px = prob.sum(axis=0)
    entropy_y = float(-np.sum(py[py > 0] * np.log2(py[py > 0])))
    if entropy_y <= 1e-12:
        return None
    mask = prob > 0
    outer = np.outer(py, px)
    mi = float(np.sum(prob[mask] * np.log2(prob[mask] / outer[mask])))
    leak = float(max(0.0, 1.0 - mi / entropy_y))
    return mi, leak, entropy_y


def scan_continuum(lines, cont, sampling, config):
    targets = config['targets']
    edges = {name: marginal_edges(lines[name], config['bins']) for name in targets}
    edges['continuum'] = marginal_edges(cont.continuum, config['bins'])
    rows = []
    for method in sampling['methods']:
        for lag in range(sampling['lag_min_days'], sampling['lag_max_days'] + 1, sampling['lag_step_days']):
            frame = tuples(lines, cont, lag, method, sampling)
            for target in targets:
                row = dict(method=method, lag_days=lag, target=target,
                           predictor='continuum', tuple_count=len(frame),
                           line_observed_dates=len(lines), continuum_observed_dates=len(cont),
                           status='insufficient_tuples')
                if len(frame) >= config['minimum_tuples']:
                    cols = ['target_' + target, 'continuum']
                    hist, _ = np.histogramdd(frame[cols].to_numpy(), bins=[edges[target], edges['continuum']])
                    result = decompose_1d(hist)
                    row.update(histogram_cells=hist.size, occupied_cells=np.count_nonzero(hist),
                               empty_fraction=float(np.mean(hist == 0)),
                               singleton_tuple_fraction=float(np.sum(hist == 1) / len(frame)))
                    if result is None:
                        row['status'] = 'constant_target'
                    else:
                        mi, leak, entropy_y = result
                        row.update(status='ok', mutual_info_bits=mi,
                                   information_leak_fraction=leak,
                                   target_entropy_bits=entropy_y)
                rows.append(row)
    return pd.DataFrame(rows)


def peaks_continuum(curves):
    rows = []
    for (method, target), group in curves.groupby(['method', 'target']):
        valid = group[group.status == 'ok'].dropna(subset=['mutual_info_bits'])
        if valid.empty:
            continue
        max_mi_row = valid.loc[valid['mutual_info_bits'].idxmax()]
        min_leak_row = valid.loc[valid['information_leak_fraction'].idxmin()]
        rows.append(dict(method=method, target=target,
                         peak_mi_lag_days=int(max_mi_row.lag_days),
                         peak_mi_value=float(max_mi_row.mutual_info_bits),
                         min_leak_lag_days=int(min_leak_row.lag_days),
                         min_leak_value=float(min_leak_row.information_leak_fraction),
                         peak_tuple_count=int(max_mi_row.tuple_count),
                         scanned_lags=len(valid)))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description='Run continuum-only lag recovery baseline.')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Custom output directory (default: agn_surd_project/processed/reconciled/continuum_baseline)')
    args = parser.parse_args()

    sampling = json.loads((HERE / 'sampling_config.json').read_text())
    config = json.loads((HERE / 'information_config.json').read_text())

    out = Path(args.output_dir) if args.output_dir else (
        ROOT / 'agn_surd_project/processed/reconciled/continuum_baseline'
    )
    out.mkdir(parents=True, exist_ok=True)

    print("Loading adopted data...")
    lines, cont = load_adopted(sampling['continuum'])

    print("Running observed continuum-only scan...")
    observed = scan_continuum(lines, cont, sampling, config)
    obs_peaks = peaks_continuum(observed)
    observed.to_csv(out / 'observed_continuum_curves.csv', index=False)
    obs_peaks.to_csv(out / 'observed_continuum_peaks.csv', index=False)
    print("Observed continuum-only peaks:")
    print(obs_peaks.to_string(index=False))

    sim_peak_rows = []
    for case in config['simulation_cases']:
        for seed in config['pilot_seeds']:
            print(f"Running simulation: {case} seed {seed}...")
            sim_lines, sim_cont = simulate(lines, cont, config, case, seed)
            sim_curves = scan_continuum(sim_lines, sim_cont, sampling, config)
            sim_peaks = peaks_continuum(sim_curves)
            sim_peaks['case'] = case
            sim_peaks['seed'] = seed
            sim_peak_rows.append(sim_peaks)
            sim_curves.to_csv(out / f'{case}_{seed}_continuum_curves.csv', index=False)

    all_sim_peaks = pd.concat(sim_peak_rows, ignore_index=True)
    all_sim_peaks.to_csv(out / 'simulated_continuum_peaks.csv', index=False)
    print("\nSummary of simulated continuum-only peaks:")
    print(all_sim_peaks[['case', 'seed', 'method', 'target', 'peak_mi_lag_days', 'peak_mi_value']].to_string(index=False))

    manifest = {
        'status': 'Continuum-only lag recovery baseline',
        'sampling_config': sampling,
        'information_config': config,
        'description': '1-predictor baseline evaluating I(Line(t); Continuum(t - tau)) over 0-200 days'
    }
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"\nOutputs saved to {out}")


if __name__ == '__main__':
    main()
