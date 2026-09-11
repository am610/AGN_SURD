"""Run full 4-target round-robin SURD decomposition exporting and plotting all 25 atoms and leakage curves."""
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .adopted_data import ROOT, load_adopted
from .information import marginal_edges, decompose
from .sampling import interpolate

HERE = Path(__file__).resolve().parent


def extract_aligned_grid(lines, cont, config, gap_limited=True):
    lt, ct = lines.jd_offset.to_numpy(), cont.jd_offset.to_numpy()
    step = config.get('grid_step_days', 1.0)
    t = np.arange(np.ceil(lt[0]), np.floor(lt[-1]) + 1, step)
    gap = config['maximum_interpolation_gap_days'] if gap_limited else np.inf

    grid_df = pd.DataFrame({'time': t})
    grid_df['continuum'] = interpolate(ct, cont.continuum.to_numpy(), t, gap)
    for col in ['blue', 'core', 'red', 'profile_total']:
        grid_df[col] = interpolate(lt, lines[col].to_numpy(), t, gap)
    return grid_df


def run_four_target_scan(lines, cont, sampling_config, info_config, mode='continuum_wings', method='gap_limited'):
    if mode == 'continuum_wings':
        var_names = ['continuum', 'blue', 'core', 'red']
    elif mode == 'total_wings':
        var_names = ['profile_total', 'blue', 'core', 'red']
    else:
        raise ValueError(f"Unknown mode: {mode}")

    gap_limited = (method == 'gap_limited')
    grid = extract_aligned_grid(lines, cont, sampling_config, gap_limited=gap_limited)

    # Precompute marginal quantile bin edges
    bins = info_config['bins']
    edges = {}
    for v in var_names:
        if v == 'continuum':
            edges[v] = marginal_edges(cont.continuum.to_numpy(), bins)
        else:
            edges[v] = marginal_edges(lines[v].to_numpy(), bins)

    lag_min = sampling_config['lag_min_days']
    lag_max = sampling_config['lag_max_days']
    lag_step = sampling_config['lag_step_days']
    lags = list(range(lag_min, lag_max + 1, lag_step))

    rows = []
    t_vals = grid['time'].to_numpy()

    for target_idx, target_name in enumerate(var_names):
        pred_names = [v for v in var_names if v != target_name]
        print(f"Scanning Target: {target_name} with predictors {pred_names} ({method})...")

        for lag in lags:
            # Query times for predictors
            query = t_vals - lag

            # Target at time t
            y_series = grid[target_name].to_numpy()

            # Predictors at time t - lag
            preds = []
            for p_name in pred_names:
                if p_name == 'continuum':
                    p_val = interpolate(cont.jd_offset.to_numpy(), cont.continuum.to_numpy(), query,
                                        sampling_config['maximum_interpolation_gap_days'] if gap_limited else np.inf)
                else:
                    p_val = interpolate(lines.jd_offset.to_numpy(), lines[p_name].to_numpy(), query,
                                        sampling_config['maximum_interpolation_gap_days'] if gap_limited else np.inf)
                preds.append(p_val)

            # Combined mask of non-NaNs
            valid = np.isfinite(y_series)
            for p in preds:
                valid &= np.isfinite(p)

            n_tuples = int(np.sum(valid))
            row = {
                'mode': mode,
                'method': method,
                'target': target_name,
                'predictors': ';'.join(pred_names),
                'lag_days': lag,
                'tuple_count': n_tuples,
                'status': 'insufficient_tuples'
            }

            if n_tuples >= info_config['minimum_tuples']:
                data_matrix = np.column_stack([y_series[valid]] + [p[valid] for p in preds])
                hist_bins = [edges[target_name]] + [edges[p] for p in pred_names]
                hist, _ = np.histogramdd(data_matrix, bins=hist_bins)
                result = decompose(hist)

                if result is not None:
                    red, syn, joint, leak, entropy = result
                    u1 = float(red.get((1,), 0.0))
                    u2 = float(red.get((2,), 0.0))
                    u3 = float(red.get((3,), 0.0))
                    r12 = float(red.get((1, 2), 0.0))
                    r13 = float(red.get((1, 3), 0.0))
                    r23 = float(red.get((2, 3), 0.0))
                    r123 = float(red.get((1, 2, 3), 0.0))
                    s12 = float(syn.get((1, 2), 0.0))
                    s13 = float(syn.get((1, 3), 0.0))
                    s23 = float(syn.get((2, 3), 0.0))
                    s123 = float(syn.get((1, 2, 3), 0.0))

                    total_u = u1 + u2 + u3
                    total_r = r12 + r13 + r23 + r123
                    total_s = s12 + s13 + s23 + s123

                    p1, p2, p3 = pred_names
                    denom = joint if joint > 1e-12 else np.nan

                    row.update({
                        'status': 'ok',
                        'joint_mi_bits': joint,
                        'target_entropy_bits': entropy,
                        'leak_fraction': leak,
                        'total_unique_bits': total_u,
                        'total_redundancy_bits': total_r,
                        'total_synergy_bits': total_s,
                        # Normalized summary fractions
                        'norm_unique': total_u / denom,
                        'norm_redundancy': total_r / denom,
                        'norm_synergy': total_s / denom,
                        # Individual atoms in bits
                        f'U_{p1}': u1, f'U_{p2}': u2, f'U_{p3}': u3,
                        f'R_{p1}_{p2}': r12, f'R_{p1}_{p3}': r13, f'R_{p2}_{p3}': r23,
                        f'R_{p1}_{p2}_{p3}': r123,
                        f'S_{p1}_{p2}': s12, f'S_{p1}_{p3}': s13, f'S_{p2}_{p3}': s23,
                        f'S_{p1}_{p2}_{p3}': s123,
                        # Normalized individual atoms
                        f'norm_U_{p1}': u1 / denom, f'norm_U_{p2}': u2 / denom, f'norm_U_{p3}': u3 / denom,
                        f'norm_R_{p1}_{p2}': r12 / denom, f'norm_R_{p1}_{p3}': r13 / denom,
                        f'norm_R_{p2}_{p3}': r23 / denom, f'norm_R_{p1}_{p2}_{p3}': r123 / denom,
                        f'norm_S_{p1}_{p2}': s12 / denom, f'norm_S_{p1}_{p3}': s13 / denom,
                        f'norm_S_{p2}_{p3}': s23 / denom, f'norm_S_{p1}_{p2}_{p3}': s123 / denom,
                    })
            rows.append(row)

    return pd.DataFrame(rows)


def plot_four_targets(df, out_path, mode_title):
    targets = df['target'].unique()
    fig, axes = plt.subplots(len(targets), 2, figsize=(14, 3.2 * len(targets)), sharex=True)
    if len(targets) == 1:
        axes = np.array([axes])

    for i, target in enumerate(targets):
        sub = df[(df.target == target) & (df.status == 'ok')]
        ax_decomp = axes[i, 0]
        ax_leak = axes[i, 1]

        # Left panel: Normalized Decomposition (U, R, S)
        ax_decomp.plot(sub.lag_days, sub.norm_unique, label='Normalized Unique (U)', color='tab:blue', lw=1.8)
        ax_decomp.plot(sub.lag_days, sub.norm_redundancy, label='Normalized Redundancy (R)', color='tab:green', lw=1.8)
        ax_decomp.plot(sub.lag_days, sub.norm_synergy, label='Normalized Synergy (S)', color='tab:red', lw=1.8)
        ax_decomp.set_ylabel(f'Target: {target}\nInformation Fraction', fontsize=11)
        ax_decomp.set_ylim(-0.05, 1.05)
        ax_decomp.grid(True, ls='--', alpha=0.5)
        if i == 0:
            ax_decomp.set_title('Normalized Information Decomposition (U, R, S)', fontsize=12)
            ax_decomp.legend(loc='upper right', framealpha=0.9, fontsize=9)

        # Right panel: Normalized Leakage & Joint MI
        ax_leak.plot(sub.lag_days, sub.leak_fraction, label='Normalized Leakage $H(Y|X)/H(Y)$', color='tab:purple', lw=1.8)
        ax_leak_tw = ax_leak.twinx()
        ax_leak_tw.plot(sub.lag_days, sub.joint_mi_bits, label='Joint MI (bits)', color='tab:gray', ls=':', lw=1.5)
        ax_leak.set_ylabel('Leakage Fraction', fontsize=11, color='tab:purple')
        ax_leak_tw.set_ylabel('Joint MI (bits)', fontsize=11, color='tab:gray')
        ax_leak.set_ylim(-0.05, 1.05)
        ax_leak.grid(True, ls='--', alpha=0.5)
        if i == 0:
            ax_leak.set_title('Information Leakage & Joint Mutual Information', fontsize=12)
            lines_1, labels_1 = ax_leak.get_legend_handles_labels()
            lines_2, labels_2 = ax_leak_tw.get_legend_handles_labels()
            ax_leak.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right', framealpha=0.9, fontsize=9)

    axes[-1, 0].set_xlabel('Lag $\\tau$ (days)', fontsize=12)
    axes[-1, 1].set_xlabel('Lag $\\tau$ (days)', fontsize=12)
    fig.suptitle(f'Four-Target SURD Information Decomposition ({mode_title})', fontsize=14, y=0.995)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved plot: {out_path}")


def main():
    parser = argparse.ArgumentParser(description='Run 4-target all-atom SURD decomposition.')
    parser.add_argument('--mode', type=str, default='all', choices=['continuum_wings', 'total_wings', 'all'],
                        help='Variable combination mode')
    parser.add_argument('--method', type=str, default='gap_limited', choices=['gap_limited', 'daily'],
                        help='Sampling method (default: gap_limited)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory')
    args = parser.parse_args()

    sampling = json.loads((HERE / 'sampling_config.json').read_text())
    config = json.loads((HERE / 'information_config.json').read_text())

    lines, cont = load_adopted(sampling['continuum'])
    out = Path(args.output_dir) if args.output_dir else (
        ROOT / 'agn_surd_project/processed/reconciled/four_targets'
    )
    out.mkdir(parents=True, exist_ok=True)

    modes = ['continuum_wings', 'total_wings'] if args.mode == 'all' else [args.mode]

    for mode in modes:
        title = "F5100 + Wings" if mode == 'continuum_wings' else "Hb_total + Wings"
        print(f"\n==================================================")
        print(f"Running Four-Target SURD Decomposition: {title}")
        print(f"==================================================")
        df = run_four_target_scan(lines, cont, sampling, config, mode=mode, method=args.method)
        csv_file = out / f'{mode}_{args.method}_all_atoms.csv'
        df.to_csv(csv_file, index=False)
        print(f"Saved table: {csv_file}")

        plot_file = out / f'{mode}_{args.method}_plots.png'
        plot_four_targets(df, plot_file, title)

    print(f"\nAll 4-target decompositions and plots completed successfully in {out}!")


if __name__ == '__main__':
    main()
