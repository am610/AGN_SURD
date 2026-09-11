"""Multi-lag and asymmetric predictor SURD scanning (Equation 3.18) for AGN kinematic models."""
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .adopted_data import ROOT, load_adopted
from .information import marginal_edges, decompose
from .sampling import interpolate
from .decompose_all_atoms import extract_aligned_grid

HERE = Path(__file__).resolve().parent

KINEMATIC_CONFIGS = {
    'symmetric_virial': {
        'description': 'Wings respond 5 days before core',
        'response_delays': {'continuum': 0, 'blue': 5, 'red': 5, 'core': 10, 'profile_total': 10}
    },
    'blue_leading_outflow': {
        'description': 'Blue responds 5 days before red',
        'response_delays': {'continuum': 0, 'blue': 5, 'core': 8, 'red': 10, 'profile_total': 8}
    },
    'red_leading_inflow': {
        'description': 'Red responds 5 days before blue',
        'response_delays': {'continuum': 0, 'red': 5, 'core': 8, 'blue': 10, 'profile_total': 8}
    },
    'common_lag_control': {
        'description': 'All predictors at the identical lag tau_0',
        'response_delays': {'continuum': 0, 'blue': 0, 'core': 0, 'red': 0, 'profile_total': 0}
    }
}


def relative_offsets(target, predictors, response_delays):
    """Return predictor lags relative to the response phase of the target."""
    target_delay = response_delays[target]
    return {name: target_delay - response_delays[name] for name in predictors}


def scan_asymmetric_lags(lines, cont, sampling_config, info_config, target='core',
                         predictors=['continuum', 'blue', 'red'], offsets={'continuum': 0, 'blue': 0, 'red': 0},
                         method='gap_limited', max_base_lag=60):
    gap_limited = (method == 'gap_limited')
    grid = extract_aligned_grid(lines, cont, sampling_config, gap_limited=gap_limited)
    t_vals = grid['time'].to_numpy()

    # Discretization edges
    bins = info_config['bins']
    edges = {}
    edges['continuum'] = marginal_edges(cont.continuum.to_numpy(), bins)
    for col in ['blue', 'core', 'red', 'profile_total']:
        edges[col] = marginal_edges(lines[col].to_numpy(), bins)

    lags = list(range(0, max_base_lag + 1))
    rows = []

    for tau_0 in lags:
        y_series = grid[target].to_numpy()

        preds = []
        for p in predictors:
            tau_p = tau_0 + offsets.get(p, 0)
            query = t_vals - tau_p
            if p == 'continuum':
                val = interpolate(cont.jd_offset.to_numpy(), cont.continuum.to_numpy(), query,
                                  sampling_config['maximum_interpolation_gap_days'] if gap_limited else np.inf)
            else:
                val = interpolate(lines.jd_offset.to_numpy(), lines[p].to_numpy(), query,
                                  sampling_config['maximum_interpolation_gap_days'] if gap_limited else np.inf)
            preds.append(val)

        valid = np.isfinite(y_series)
        for p_val in preds:
            valid &= np.isfinite(p_val)

        n_tuples = int(np.sum(valid))
        row = {
            'target': target,
            'predictors': ';'.join(predictors),
            'base_lag_days': tau_0,
            'predictor_lags': ';'.join(f"{p}:{tau_0 + offsets.get(p, 0)}" for p in predictors),
            'tuple_count': n_tuples,
            'status': 'insufficient_tuples'
        }

        if n_tuples >= info_config['minimum_tuples']:
            data_matrix = np.column_stack([y_series[valid]] + [p[valid] for p in preds])
            hist_bins = [edges[target]] + [edges[p] for p in predictors]
            hist, _ = np.histogramdd(data_matrix, bins=hist_bins)
            result = decompose(hist)

            if result is not None:
                red, syn, joint, leak, entropy = result
                u = sum(red.get((i,), 0.0) for i in range(1, len(predictors) + 1))
                r = sum(v for k, v in red.items() if len(k) > 1)
                s = sum(syn.values())
                denom = joint if joint > 1e-12 else np.nan

                row.update({
                    'status': 'ok',
                    'joint_mi_bits': joint,
                    'target_entropy_bits': entropy,
                    'leak_fraction': leak,
                    'total_unique_bits': u,
                    'total_redundancy_bits': r,
                    'total_synergy_bits': s,
                    'norm_unique': u / denom,
                    'norm_redundancy': r / denom,
                    'norm_synergy': s / denom,
                })
        rows.append(row)

    return pd.DataFrame(rows)


def run_kinematic_comparison(out_dir=None):
    sampling = json.loads((HERE / 'sampling_config.json').read_text())
    config = json.loads((HERE / 'information_config.json').read_text())

    lines, cont = load_adopted(sampling['continuum'])
    out = Path(out_dir) if out_dir else (
        ROOT / 'agn_surd_project/processed/reconciled/asymmetric_lags'
    )
    out.mkdir(parents=True, exist_ok=True)

    print("=== Scanning Asymmetric / Multi-Lag Configurations ===")

    all_frames = []
    variable_sets = {
        'continuum_wings': ['continuum', 'blue', 'core', 'red'],
        'total_wings': ['profile_total', 'blue', 'core', 'red'],
    }

    for mode, variables in variable_sets.items():
        for target in variables:
            preds = [name for name in variables if name != target]
            for name, cfg in KINEMATIC_CONFIGS.items():
                print(f"Running {mode}, target {target}, model {name}...")
                offsets = relative_offsets(target, preds, cfg['response_delays'])
                df = scan_asymmetric_lags(
                    lines, cont, sampling, config, target=target,
                    predictors=preds, offsets=offsets,
                    method='gap_limited', max_base_lag=60,
                )
                df['mode'] = mode
                df['kinematic_model'] = name
                df['configuration_description'] = cfg['description']
                all_frames.append(df)

    res_df = pd.concat(all_frames, ignore_index=True)
    csv_path = out / 'all_targets_kinematic_asymmetric_scans.csv'
    res_df.to_csv(csv_path, index=False)
    print(f"Saved results to {csv_path}")

    # The compact paper figure retains the core target comparison.  The CSV
    # contains every target and both requested four-variable sets.
    plot_df = res_df[(res_df['mode'] == 'continuum_wings') & (res_df['target'] == 'core')]
    legacy_csv_path = out / 'core_target_kinematic_asymmetric_scans.csv'
    plot_df.to_csv(legacy_csv_path, index=False)
    print(f"Saved corrected core comparison to {legacy_csv_path}")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    colors = {'common_lag_control': 'black', 'symmetric_virial': 'tab:blue',
              'blue_leading_outflow': 'tab:cyan', 'red_leading_inflow': 'tab:red'}

    for name, group in plot_df.groupby('kinematic_model'):
        sub = group[group.status == 'ok']
        lbl = name.replace('_', ' ').title()
        ax1.plot(sub.base_lag_days, sub.norm_synergy, label=lbl, color=colors.get(name, 'gray'), lw=1.8)
        ax2.plot(sub.base_lag_days, sub.leak_fraction, label=lbl, color=colors.get(name, 'gray'), lw=1.8)

    ax1.set_xlabel('Base Lag $\\tau_0$ (days)', fontsize=12)
    ax1.set_ylabel('Normalized Synergy $S / I_{\\mathrm{joint}}$', fontsize=12)
    ax1.set_title('Normalized Synergy for Core Target under Multi-Lag Scenarios', fontsize=12)
    ax1.grid(True, ls='--', alpha=0.5)
    ax1.legend(loc='upper right', fontsize=10)

    ax2.set_xlabel('Base Lag $\\tau_0$ (days)', fontsize=12)
    ax2.set_ylabel('Information Leakage $H(Y|X)/H(Y)$', fontsize=12)
    ax2.set_title('Information Leakage for Core Target under Multi-Lag Scenarios', fontsize=12)
    ax2.grid(True, ls='--', alpha=0.5)
    ax2.legend(loc='lower right', fontsize=10)

    plt.tight_layout()
    plot_path = out / 'core_target_kinematic_asymmetric_comparison.png'
    fig.savefig(plot_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved plot: {plot_path}")


def main():
    parser = argparse.ArgumentParser(description='Run asymmetric predictor lag scans.')
    parser.add_argument('--output-dir', type=str, default=None, help='Output directory')
    args = parser.parse_args()
    run_kinematic_comparison(out_dir=args.output_dir)


if __name__ == '__main__':
    main()
