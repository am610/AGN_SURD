"""Per epoch SURD scans on a regular daily grid, common and uncommon predictor lags.

Each campaign is analysed independently. Every one of the five series is
linearly interpolated onto the campaign's daily grid with no extrapolation and
no bridging of gaps longer than the declared limit. The target is Y(t); the
predictors are X_i(t - L - delta_i) for the four members of a subset, including
the target's own history. delta_i = 0 for every predictor is the common lag
analysis. A nonzero delta_i gives that predictor its own lag, which generalises
the lag set of thesis equation 3.18 to one lag per predictor variable.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from collaborator_rebuild.official_decomposition import HASHES, decompose_counts
from collaborator_rebuild.lag_design import CAMPAIGNS
from pipeline.sampling import interpolate

BASE = Path(__file__).resolve().parent
NATIVE = BASE / 'published_results' / 'bins_2'
OUT = BASE / 'epoch_results'
SUBSETS = {'total_wings': ['total', 'blue', 'core', 'red'],
           'continuum_wings': ['continuum', 'blue', 'core', 'red']}
# Days added to the base lag of a predictor. The earlier responding wing is
# sampled further in the past so that it carries the same driving event.
CONFIGS = {
    'common': {'total': 0, 'continuum': 0, 'blue': 0, 'core': 0, 'red': 0},
    'blue_leading': {'total': 0, 'continuum': 0, 'blue': 5, 'core': 0, 'red': 0},
    'red_leading': {'total': 0, 'continuum': 0, 'blue': 0, 'core': 0, 'red': 5},
}
CONFIG = {
    'epochs': CAMPAIGNS,
    'native_input': 'published_results/bins_2/<epoch>/native/<series>.csv',
    'grid': 'Daily integer grid of JD less 2400000 spanning the epoch; one value per day',
    'interpolation': 'Linear; no extrapolation; no bridging of a bracketing gap above the limit',
    'maximum_interpolation_gap_days': 15,
    'lag_days': list(range(1, 31)),
    'bins': [2, 3],
    'bin_edges': 'Marginal quantiles of each full native epoch series, shared by target and predictor roles',
    'minimum_tuples': 32,
    'predictor_lag': 'base lag + delta of that predictor, in observed days',
    'predictor_offsets_days': CONFIGS,
    'subsets': SUBSETS,
    'targets': 'Every member of each subset in turn; all four members are past predictors',
    'normalization': 'Atoms divided by joint mutual information; leakage divided by target entropy',
    'estimator': 'Unchanged authors SURD routine at the pinned source hashes',
    'interpretation': 'Descriptive SURD estimates; no significance or causal detection claim',
}


def load_epoch(year):
    return {name: pd.read_csv(NATIVE / str(year) / 'native' / f'{name}.csv')
            for name in ['continuum', 'total', 'blue', 'core', 'red']}


def epoch_grid(series):
    lo = max(f.jd_offset.min() for f in series.values())
    hi = min(f.jd_offset.max() for f in series.values())
    return np.arange(np.ceil(lo), np.floor(hi) + 1.0)


def lag_tuples(series, names, target, grid, lag, offsets, maximum_gap):
    """Rows of (Y(t), X_1, ..., X_4) with every entry finite."""
    def sample(name, times):
        frame = series[name]
        return interpolate(frame.jd_offset.to_numpy(), frame.value.to_numpy(), times, maximum_gap)
    columns = [sample(target, grid)]
    columns += [sample(n, grid - lag - offsets[n]) for n in names]
    data = np.column_stack(columns)
    valid = np.isfinite(data).all(axis=1)
    return data[valid], grid[valid]


def quantile_edges(series, bins):
    return {n: np.r_[-np.inf, np.unique(np.quantile(f.value, np.arange(1, bins) / bins)), np.inf]
            for n, f in series.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    (OUT / 'configuration.json').write_text(json.dumps(CONFIG, indent=2) + '\n')
    atoms, summaries, inventory = [], [], []
    gap = CONFIG['maximum_interpolation_gap_days']
    for year in CAMPAIGNS:
        series = load_epoch(year)
        grid = epoch_grid(series)
        for name, frame in series.items():
            inventory.append(dict(epoch=year, series=name, native_dates=len(frame),
                                  first_jd=float(frame.jd_offset.min()), last_jd=float(frame.jd_offset.max())))
        edges = {b: quantile_edges(series, b) for b in CONFIG['bins']}
        for subset, names in SUBSETS.items():
            for target in names:
                for config, offsets in CONFIGS.items():
                    for lag in CONFIG['lag_days']:
                        data, dates = lag_tuples(series, names, target, grid, lag, offsets, gap)
                        for bins in CONFIG['bins']:
                            row = dict(epoch=year, subset=subset, target=target, config=config,
                                       bins=bins, lag_days=lag, tuple_count=len(data), status='insufficient_support')
                            if len(data) >= CONFIG['minimum_tuples']:
                                counts, _ = np.histogramdd(data, bins=[edges[bins][n] for n in [target] + names])
                                row.update(occupied_cells=int(np.count_nonzero(counts)), histogram_cells=int(counts.size))
                                if np.count_nonzero(counts.sum(axis=tuple(range(1, counts.ndim)))) < 2:
                                    row['status'] = 'degenerate_target'
                                else:
                                    result = decompose_counts(counts, names)
                                    row.update(status='ok' if result['normalization_defined'] else 'undefined_normalization',
                                               joint_mi_bits=result['joint_mi_bits'],
                                               normalized_leakage=result['normalized_leakage'])
                                    if row['status'] == 'ok':
                                        for atom in result['atoms']:
                                            atoms.append(dict(epoch=year, subset=subset, target=target, config=config,
                                                              bins=bins, lag_days=lag, kind=atom['kind'],
                                                              component=atom['kind'] + '(' + ','.join(atom['predictors']) + ')',
                                                              bits=atom['bits'], fraction=atom['fraction']))
                            summaries.append(row)
    summary = pd.DataFrame(summaries)
    components = pd.DataFrame(atoms)
    summary.to_csv(OUT / 'support_and_leakage.csv', index=False)
    components.to_csv(OUT / 'individual_components.csv.gz', index=False)
    pd.DataFrame(inventory).to_csv(OUT / 'native_inventory.csv', index=False)
    inputs = sorted(NATIVE.glob('*/native/*.csv'))
    manifest = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    configuration_sha256=hashlib.sha256((OUT / 'configuration.json').read_bytes()).hexdigest(),
                    official_source_sha256=HASHES,
                    input_sha256={str(p.relative_to(NATIVE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                    scan_rows=len(summary), usable_scan_rows=int((summary.status == 'ok').sum()),
                    atom_rows=len(components), status_counts=summary.status.value_counts().to_dict())
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({k: v for k, v in manifest.items() if k != 'input_sha256'}, indent=2))
    if not args.no_plots:
        from collaborator_rebuild.plot_epoch import plot_all
        plot_all()


if __name__ == '__main__':
    main()
