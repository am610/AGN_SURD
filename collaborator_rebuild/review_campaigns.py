"""Review saved published scans after verifying identical native support."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from collaborator_rebuild.compare_bins import KEYS

BASE = Path(__file__).resolve().parent / 'published_results'
CURVES = KEYS[:-1]


def verify_inputs(base):
    configs = [json.loads((base / f'bins_{b}' / 'configuration.json').read_text()) for b in (2, 3)]
    manifests = [json.loads((base / f'bins_{b}' / 'manifest.json').read_text()) for b in (2, 3)]
    if manifests[0]['source_sha256'] != manifests[1]['source_sha256']:
        raise ValueError('Scans used different runner sources')
    for b, manifest in zip((2, 3), manifests):
        actual = hashlib.sha256((base / f'bins_{b}' / 'configuration.json').read_bytes()).hexdigest()
        if manifest['config_sha256'] != actual:
            raise ValueError('Saved configuration does not match scan manifest')
    if [c.pop('bins') for c in configs] != [2, 3] or configs[0] != configs[1]:
        raise ValueError('Scan settings differ beyond bin count')
    files = [sorted(p.relative_to(base / f'bins_{b}') for p in
                    (base / f'bins_{b}').glob('*/native/*.csv')) for b in (2, 3)]
    if not files[0] or files[0] != files[1]:
        raise ValueError('Native input inventories differ or are empty')
    for path in files[0]:
        if (base / 'bins_2' / path).read_bytes() != (base / 'bins_3' / path).read_bytes():
            raise ValueError(f'Native inputs differ: {path}')
    support = [pd.read_csv(base / f'bins_{b}' / 'support_and_leakage.csv') for b in (2, 3)]
    columns = KEYS + ['tuple_count', 'status']
    pd.testing.assert_frame_equal(*[f[columns].sort_values(KEYS).reset_index(drop=True) for f in support])
    if any(f.duplicated(KEYS).any() for f in support):
        raise ValueError('Repeated scan configuration')
    return configs[0], support, len(files[0])


def compare_atoms(frames, support):
    atom_keys = KEYS + ['component', 'kind']
    paired = frames[0].merge(frames[1], on=atom_keys, suffixes=('_2', '_3'),
                             how='outer', indicator=True, validate='one_to_one')
    if not paired['_merge'].eq('both').all():
        raise ValueError('Atom inventories differ')
    paired = paired[paired.lag_days > 0].copy()
    if not np.isfinite(paired[['fraction_2', 'fraction_3']]).all().all():
        raise ValueError('Undefined atom fractions')
    expected = support.query('lag_days > 0 and status == "ok"')[KEYS]
    counts = paired.groupby(KEYS).size().rename('atom_count').reset_index()
    check = expected.merge(counts, on=KEYS, how='outer', indicator=True, validate='one_to_one')
    if not check['_merge'].eq('both').all() or not check.atom_count.eq(26).all():
        raise ValueError('Expected 26 atoms per supported configuration')
    for b in (2, 3):
        sums = paired.groupby(KEYS)[f'fraction_{b}'].sum()
        if not np.allclose(sums, 1, atol=1e-10, rtol=0):
            raise ValueError('Atom fractions do not sum to one')
    paired['absolute_change'] = abs(paired.fraction_2 - paired.fraction_3)
    changes = (paired.groupby(KEYS).absolute_change.sum() / 2).rename('total_variation').reset_index()
    leak = paired.drop_duplicates(KEYS)[KEYS + ['normalized_leakage_2', 'normalized_leakage_3']]
    changes = changes.merge(leak, on=KEYS, validate='one_to_one')
    changes['absolute_leakage_change'] = abs(changes.normalized_leakage_2 - changes.normalized_leakage_3)
    return paired, changes


def peak_review(paired):
    rows = []
    for key, group in paired.groupby(CURVES + ['component', 'kind']):
        row = dict(zip(CURVES + ['component', 'kind'], key))
        row['supported_lags'] = len(group)
        row['median_absolute_fraction_change'] = group.absolute_change.median()
        peaks = []
        for b in (2, 3):
            peak = group[f'fraction_{b}'].max()
            lags = group.loc[np.isclose(group[f'fraction_{b}'], peak, atol=1e-12, rtol=0), 'lag_days'].astype(int).tolist()
            row[f'maximum_fraction_{b}'] = peak
            row[f'peak_lags_{b}'] = json.dumps(sorted(lags))
            peaks.append(set(lags))
        row['peak_sets_overlap'] = bool(peaks[0] & peaks[1])
        row['minimum_peak_separation_days'] = min(abs(a - b) for a in peaks[0] for b in peaks[1])
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    config, supports, native_count = verify_inputs(BASE)
    frames = [pd.read_csv(BASE / f'bins_{b}' / 'individual_components.csv') for b in (2, 3)]
    paired, changes = compare_atoms(frames, supports[0])
    peaks = peak_review(paired)
    out = BASE / 'campaign_review'
    out.mkdir(exist_ok=True)
    changes.to_csv(out / 'configuration_changes.csv', index=False)
    paired.to_csv(out / 'individual_atom_changes.csv', index=False)
    peaks.to_csv(out / 'component_peak_review.csv', index=False)
    rows = []
    for year, group in supports[0].query('lag_days > 0').groupby('epoch'):
        matched = changes[changes.epoch == year]
        ok = group[group.status == 'ok']
        row = dict(epoch=int(year), configurations=len(group), supported_configurations=len(ok),
                   supported_curves=len(ok[CURVES].drop_duplicates()),
                   minimum_tuples=int(group.tuple_count.min()), maximum_tuples=int(group.tuple_count.max()),
                   median_total_variation=matched.total_variation.median(),
                   maximum_total_variation=matched.total_variation.max(),
                   fraction_tv_above_point_two=(matched.total_variation > .2).mean() if len(matched) else np.nan,
                   median_absolute_leakage_change=matched.absolute_leakage_change.median())
        for b, support in zip((2, 3), supports):
            selected = support.query('epoch == @year and lag_days > 0 and status == "ok"')
            row[f'median_occupied_cells_{b}'] = selected.occupied_cells.median()
            row[f'median_tuples_per_occupied_cell_{b}'] = (selected.tuple_count / selected.occupied_cells).median()
        eligible = peaks[(peaks.epoch == year) & (peaks.supported_lags >= 2)]
        row['component_curves_with_multiple_lags'] = len(eligible)
        row['fraction_disjoint_peak_sets'] = (~eligible.peak_sets_overlap).mean() if len(eligible) else np.nan
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(out / 'campaign_summary.csv', index=False)
    fig, axes = plt.subplots(len(rows), 1, figsize=(11, 13), squeeze=False)
    for ax, year in zip(axes[:, 0], summary.epoch):
        grid = changes[changes.epoch == year].pivot(index=['subset', 'target', 'scenario'], columns='lag_days', values='total_variation')
        inventory = supports[0].query('epoch == @year')[['subset', 'target', 'scenario']].drop_duplicates()
        index = pd.MultiIndex.from_frame(inventory)
        grid = grid.reindex(index=index, columns=[v for v in config['lag_days'] if v > 0])
        cmap = plt.get_cmap('viridis').copy()
        cmap.set_bad('lightgray')
        im = ax.imshow(grid.to_numpy(), vmin=0, vmax=1, aspect='auto', cmap=cmap,
                       extent=(.5, 30.5, len(grid) - .5, -.5))
        ax.set_yticks(range(len(grid)), [' / '.join(key) for key in grid.index], fontsize=5)
        ax.set_title(f'{year}: bin sensitivity on identical native support')
        ax.set_xlabel('Base lag in observed days')
    fig.colorbar(im, ax=axes[:, 0].tolist(), label='Total variation across 26 normalized components', fraction=.025)
    fig.suptitle('Gray cells fail the 32 tuple numerical guard; descriptive comparison only')
    fig.subplots_adjust(left=.27, right=.86, top=.94, bottom=.05, hspace=.65)
    fig.savefig(out / 'campaign_sensitivity.png', dpi=180)
    plt.close(fig)
    inputs = [BASE / f'bins_{b}' / filename for b in (2, 3) for filename in
              ('configuration.json', 'manifest.json', 'individual_components.csv', 'support_and_leakage.csv')]
    inputs += sorted((BASE / 'bins_2').glob('*/native/*.csv'))
    manifest = dict(native_files_verified=native_count, matched_positive_lag_configurations=len(changes),
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    input_sha256={str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                    interpretation='Descriptive diagnostics; no significance test or physical delay estimate')
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
