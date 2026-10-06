"""Run individual SURD components within explicitly defined historical epochs."""

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline.adopted_data import load_adopted
from pipeline.sampling import interpolate
from collaborator_rebuild.official_decomposition import decompose_counts
from collaborator_rebuild.lag_design import campaign_labels, CAMPAIGNS

BASE = Path(__file__).resolve().parent
OUT = BASE / 'historical_results'
CONFIG = {
    'epoch_rule': 'Published intervals from Peterson et al. 2004 Table 2',
    'campaign_jd_offsets': CAMPAIGNS,
    'lag_days': list(range(31)),
    'maximum_interpolation_gap_days': 10,
    'bins': 3,
    'minimum_native_targets': 32,
    'target_sampling': 'Native dates of the target; predictors interpolated within the same epoch',
    'time_frame': 'Observed days',
    'total_definition': 'Broad profile flux within adopted velocity windows; blue plus core plus red',
    'velocity_km_s': {'blue': [-6000, -2000], 'core': [-2000, 2000], 'red': [2000, 6000]},
    'unequal_offsets_days': {
        'common': {'total': 0, 'continuum': 0, 'blue': 0, 'core': 0, 'red': 0},
        'blue_earlier': {'total': 0, 'continuum': 0, 'blue': 5, 'core': 0, 'red': 0},
        'red_earlier': {'total': 0, 'continuum': 0, 'blue': 0, 'core': 0, 'red': 5},
    },
    'unequal_scope': 'Exploratory one delay per predictor; not the expanded vector of thesis equation 3.18',
    'interpretation': 'Descriptive SURD estimates; no significance or causal detection claim',
}


def epoch(times):
    return campaign_labels(times)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Freeze settings before computing information curves.
    (OUT / 'configuration.json').write_text(json.dumps(CONFIG, indent=2) + '\n')
    lines, continuum = load_adopted()
    sources = {'continuum': continuum.rename(columns={'continuum': 'value', 'continuum_error': 'error'})}
    for name, column in [('total', 'profile_total'), ('blue', 'blue'), ('core', 'core'), ('red', 'red')]:
        sources[name] = lines.rename(columns={column: 'value', column + '_error': 'error'})
    sources = {n: f[['jd_offset', 'value', 'error']].assign(epoch=epoch(f.jd_offset)) for n, f in sources.items()}
    atoms, summaries, inventory = [], [], []
    for year in sorted(sources['total'].epoch.unique()):
        series = {n: f[f.epoch == year].copy() for n, f in sources.items()}
        native = OUT / str(year) / 'native'
        native.mkdir(parents=True, exist_ok=True)
        for name, frame in series.items():
            frame.to_csv(native / f'{name}.csv', index=False)
            inventory.append(dict(epoch=int(year), series=name, observations=len(frame)))
        edges = {n: np.r_[-np.inf, np.unique(np.quantile(f.value, [1/3, 2/3])), np.inf]
                 for n, f in series.items() if len(f)}
        for subset, names in [('total_wings', ['total', 'blue', 'core', 'red']),
                              ('continuum_wings', ['continuum', 'blue', 'core', 'red'])]:
            if any(len(series[n]) < 2 for n in names):
                continue
            for target in names:
                t = series[target].jd_offset.to_numpy()
                y = series[target].value.to_numpy()
                for scenario, offsets in CONFIG['unequal_offsets_days'].items():
                    for lag in CONFIG['lag_days']:
                        predictors = [interpolate(series[n].jd_offset.to_numpy(), series[n].value.to_numpy(),
                                                  t - lag - offsets[n], 10) for n in names]
                        data = np.column_stack([y] + predictors)
                        valid = np.isfinite(data).all(axis=1)
                        row = dict(epoch=int(year), subset=subset, target=target, scenario=scenario,
                                   lag_days=lag, tuple_count=int(valid.sum()), status='insufficient_support')
                        if valid.sum() >= CONFIG['minimum_native_targets']:
                            hist, _ = np.histogramdd(data[valid], bins=[edges[n] for n in [target] + names])
                            result = decompose_counts(hist, names)
                            row.update(status='ok', joint_mi_bits=result['joint_mi_bits'],
                                       normalized_leakage=result['normalized_leakage'],
                                       occupied_cells=int(np.count_nonzero(hist)), histogram_cells=hist.size)
                            for atom in result['atoms']:
                                atoms.append(dict(**row, component=atom['kind'] + '(' + ','.join(atom['predictors']) + ')',
                                                  kind=atom['kind'], bits=atom['bits'], fraction=atom['fraction']))
                        summaries.append(row)
    summary = pd.DataFrame(summaries)
    components = pd.DataFrame(atoms)
    summary.to_csv(OUT / 'support_and_leakage.csv', index=False)
    components.to_csv(OUT / 'individual_components.csv', index=False)
    pd.DataFrame(inventory).to_csv(OUT / 'native_inventory.csv', index=False)
    if not components.empty:
        for key, group in components.groupby(['epoch', 'subset', 'target', 'scenario']):
            fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
            for ax, kind in zip(axs[:3], ['U', 'R', 'S']):
                for label, curve in group[group.kind == kind].groupby('component'):
                    ax.plot(curve.lag_days, curve.fraction, label=label)
                ax.set_ylabel(kind + ' / joint MI')
                ax.legend(fontsize=7, ncol=3)
                ax.set_ylim(-0.02, 1.02)
            leak = group.drop_duplicates('lag_days')
            axs[3].plot(leak.lag_days, leak.normalized_leakage, color='black')
            axs[3].set_ylabel('Normalized leakage')
            axs[3].set_ylim(-0.02, 1.02)
            axs[3].set_xlabel('Base lag in observed days')
            fig.suptitle(' '.join(map(str, key)) + ' descriptive SURD')
            fig.tight_layout()
            path = OUT / str(key[0]) / ('_'.join(key[1:]) + '.png')
            fig.savefig(path, dpi=140)
            plt.close(fig)
    manifest = {'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'config_sha256': hashlib.sha256((OUT / 'configuration.json').read_bytes()).hexdigest(),
                'native_inventory_rows': len(inventory), 'scan_rows': len(summary),
                'usable_scan_rows': int((summary.status == 'ok').sum()),
                'atom_rows': len(components)}
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
