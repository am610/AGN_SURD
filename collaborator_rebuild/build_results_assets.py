"""Freeze descriptive tables and publication figures for the results section."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'collaborator_rebuild'
OUT = BASE / 'results_section'
KEYS = ['epoch', 'subset', 'target', 'scenario']
SYMBOLS = {'continuum': 'F', 'blue': 'B', 'core': 'C', 'red': 'R'}
SCENARIOS = {'common': 'Common', 'blue_earlier': 'Blue earlier', 'red_earlier': 'Red earlier'}


def short_component(label):
    for name, symbol in SYMBOLS.items():
        label = label.replace(name, symbol)
    return label


def write_table(name, rows, columns, header):
    lines = [r'\begin{tabular}{' + columns + '}', r'\toprule', header + r' \\', r'\midrule']
    lines += [' & '.join(map(str, row)) + r' \\' for row in rows]
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / name).write_text('\n'.join(lines) + '\n')


def main():
    OUT.mkdir(exist_ok=True)
    (OUT / 'figures').mkdir(exist_ok=True)
    paths = {
        'campaigns': BASE / 'published_results/campaign_review/campaign_summary.csv',
        'fixed': BASE / 'published_results/fixed_support/curve_comparison.csv',
        'atoms': BASE / 'published_results/fixed_support/individual_components.csv',
        'support': BASE / 'published_results/fixed_support/curve_support.csv',
        'previous': BASE / 'published_results/fixed_support/support_change_summary.csv',
        'dates': BASE / 'published_results/fixed_support/fixed_target_dates.csv',
        'inventory': BASE / 'published_results/bins_2/native_inventory.csv',
        'preparation': BASE / 'historical_prepared/manifest.json',
        'configuration': BASE / 'published_results/fixed_support/configuration.json',
        'scan_manifest': BASE / 'published_results/fixed_support/manifest.json',
        'calibration': BASE / 'calibration_review/summary.csv',
        'calibration_manifest': BASE / 'calibration_review/manifest.json',
        'order_audit': BASE / 'order_review/curve_summary.csv',
        'order_manifest': BASE / 'order_review/manifest.json',
    }
    campaign, fixed, atoms, support, inventory = [pd.read_csv(paths[name]) for name in
                                                 ['campaigns', 'fixed', 'atoms', 'support', 'inventory']]
    assert len(atoms) == 6240 and len(fixed) == 4 and len(support) == 120
    assert atoms.groupby(KEYS + ['bins']).lag_days.nunique().eq(30).all()
    assert atoms.groupby(KEYS + ['bins']).tuple_count.nunique().eq(1).all()
    assert atoms.groupby(KEYS + ['bins']).target_dates_sha256.nunique().eq(1).all()
    assert np.allclose(atoms.groupby(KEYS + ['bins', 'lag_days']).fraction.sum(), 1, atol=1e-10, rtol=0)
    ranges = support.groupby('epoch').fixed_tuple_count.agg(['min', 'max'])
    rows = []
    for row in campaign.itertuples():
        year = row.epoch
        counts = inventory[inventory.epoch == year].set_index('series').observations
        rows.append([year, int(counts['total']), int(counts['continuum']),
                     row.supported_configurations, f'{row.median_total_variation:.3f}',
                     f'{ranges.loc[year, "min"]} to {ranges.loc[year, "max"]}'])
    write_table('campaign_rows.tex', rows, 'rrrrrl',
                r'Campaign & $N_{\rm line}$ & $N_F$ & $N_{\rm usable}$ & Median $D_{\rm TV}$ & Fixed dates')
    rows = []
    fixed = fixed.sort_values(['epoch', 'scenario'])
    peaks = pd.read_csv(BASE / 'published_results/fixed_support/component_peak_review.csv')
    for row in fixed.itertuples():
        selected = peaks[(peaks.epoch == row.epoch) & (peaks.scenario == row.scenario)]
        disjoint = int((~selected.peak_sets_overlap).sum())
        rows.append([row.epoch, SCENARIOS[row.scenario], row.fixed_tuple_count,
                     f'{row.median_total_variation:.3f}', f'{row.median_absolute_leakage_change:.3f}',
                     f'{disjoint}/26'])
    paths['peaks'] = BASE / 'published_results/fixed_support/component_peak_review.csv'
    write_table('fixed_rows.tex', rows, 'rlrrrr',
                r'Campaign & Scenario & Dates & Median $D_{\rm TV}$ & Median $|\Delta\ell|$ & Disjoint peaks')
    extrema = []
    for key, group in atoms.groupby(KEYS + ['bins', 'kind']):
        maximum = group.fraction.max()
        selected = group[np.isclose(group.fraction, maximum, atol=1e-12, rtol=0)]
        for component, part in selected.groupby('component'):
            extrema.append(dict(zip(KEYS + ['bins', 'kind'], key), component=component,
                                lag_days=','.join(map(str, sorted(part.lag_days.astype(int)))),
                                maximum_fraction=float(maximum)))
    extrema = pd.DataFrame(extrema)
    extrema.to_csv(OUT / 'descriptive_extrema.csv', index=False)
    representative = extrema[(extrema.epoch == 1993) & (extrema.scenario == 'common')]
    representative = representative.assign(order=representative.kind.map({'U': 0, 'R': 1, 'S': 2})).sort_values(['bins', 'order'])
    write_table('extrema_rows.tex', [[r.bins, r.kind, short_component(r.component), r.lag_days,
                                    f'{r.maximum_fraction:.3f}'] for r in representative.itertuples()],
                'rllrr', r'Bins & Type & Component & Base lag (days) & Maximum fraction')
    calibration = pd.read_csv(paths['calibration']).query('model == "independent"')
    rows = []
    for count in sorted(calibration.sample_size.unique()):
        selected = calibration[calibration.sample_size == count].set_index('bins')
        rows.append([count, f'{selected.loc[2, "median_joint_mi_bits"]:.3f}',
                     f'{selected.loc[3, "median_joint_mi_bits"]:.3f}',
                     f'{selected.loc[2, "median_leakage"]:.3f}',
                     f'{selected.loc[3, "median_leakage"]:.3f}'])
    write_table('calibration_rows.tex', rows, 'rrrrr',
                r'Tuples & $I$, two states & $I$, three states & $\ell$, two states & $\ell$, three states')
    order_audit = pd.read_csv(paths['order_audit'])
    rows = [[r.epoch, SCENARIOS[r.scenario], r.bins, f'{r.affected_configurations}/30',
             f'{r.maximum_total_variation:.3f}' if r.affected_configurations else r'$<10^{-8}$']
            for r in order_audit.itertuples()]
    write_table('order_rows.tex', rows, 'rlrrr', r'Campaign & Scenario & Bins & Affected lags & Maximum $D_{\rm TV}$')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10,
                         'axes.labelsize': 9, 'legend.fontsize': 7, 'xtick.labelsize': 8, 'ytick.labelsize': 8})
    figures = []
    for key, group in atoms.groupby(KEYS):
        year, subset, target, scenario = key
        fig, axes = plt.subplots(4, 1, figsize=(6.8, 8), sharex=True)
        for ax, kind, title in zip(axes[:3], ['U', 'R', 'S'], ['Unique', 'Redundant', 'Synergistic']):
            labels = sorted(group[group.kind == kind].component.unique())
            colors = plt.get_cmap('tab20')(np.linspace(0, .95, len(labels)))
            for color, component in zip(colors, labels):
                for bins, style in [(2, ':'), (3, '-')]:
                    curve = group[(group.component == component) & (group.bins == bins)].sort_values('lag_days')
                    ax.plot(curve.lag_days, curve.fraction, color=color, linestyle=style, linewidth=1.25,
                            label=short_component(component) if bins == 3 else None)
            ax.set_title(title, loc='left')
            ax.set_ylabel('Fraction of joint MI')
            ax.set_ylim(-.025, 1.025)
            ax.legend(ncol=3 if kind != 'U' else 4, loc='upper right', columnspacing=.8, handlelength=1.6)
            ax.grid(alpha=.15)
        for bins, style in [(2, ':'), (3, '-')]:
            curve = group[group.bins == bins].drop_duplicates('lag_days').sort_values('lag_days')
            axes[3].plot(curve.lag_days, curve.normalized_leakage, color='black', linestyle=style,
                         linewidth=1.25, label=f'{bins} bins')
        axes[3].set_title('Information leakage', loc='left')
        axes[3].set_ylabel('Fraction of target entropy')
        axes[3].set_ylim(-.025, 1.025)
        axes[3].set_xlabel('Base lag in observed days')
        axes[3].set_xlim(1, 30)
        axes[3].set_xticks([1, 5, 10, 15, 20, 25, 30])
        axes[3].grid(alpha=.15)
        axes[3].legend(loc='upper right')
        fig.suptitle(f'{year}: {SCENARIOS[scenario].lower()}, continuum target, {int(group.tuple_count.iloc[0])} fixed dates', y=.995)
        fig.tight_layout(rect=(0, 0, 1, .98), h_pad=.7)
        path = OUT / 'figures' / f'{year}_{scenario}.png'
        fig.savefig(path, dpi=240)
        plt.close(fig)
        figures.append(path)
    generated = [OUT / name for name in ['campaign_rows.tex', 'fixed_rows.tex', 'extrema_rows.tex',
                                       'calibration_rows.tex', 'order_rows.tex', 'descriptive_extrema.csv']] + figures
    manifest = dict(build_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()},
                    output_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in generated},
                    purpose='Frozen descriptive results section; initial calibration and ordering audit complete; realistic observation calibration unperformed')
    (OUT / 'asset_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Created {len(generated)} assets with verified source tables.')


if __name__ == '__main__':
    main()
