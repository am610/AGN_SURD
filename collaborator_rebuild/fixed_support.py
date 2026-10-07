"""Compare published SURD curves on a fixed intersection of target dates."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from collaborator_rebuild.official_decomposition import HASHES, decompose_counts
from collaborator_rebuild.compare_bins import KEYS
from collaborator_rebuild.review_campaigns import BASE, CURVES, compare_atoms, peak_review, verify_inputs
from pipeline.sampling import interpolate

SUBSETS = {'total_wings': ['total', 'blue', 'core', 'red'],
           'continuum_wings': ['continuum', 'blue', 'core', 'red']}
OUT = BASE / 'fixed_support'


def fixed_lag_tuples(series, names, target, lags, offsets, maximum_gap):
    """Return full lag matrices and one intersection mask for native targets."""
    lags = tuple(lags)
    if not lags or len(set(lags)) != len(lags) or any(lag <= 0 for lag in lags):
        raise ValueError('Require distinct positive lags')
    times = series[target].jd_offset.to_numpy()
    y = series[target].value.to_numpy()
    matrices = []
    for lag in lags:
        predictors = [interpolate(series[n].jd_offset.to_numpy(), series[n].value.to_numpy(),
                                  times - lag - offsets[n], maximum_gap) for n in names]
        matrices.append(np.column_stack([y] + predictors))
    valid = np.stack([np.isfinite(data).all(axis=1) for data in matrices])
    return matrices, valid.all(axis=0), valid.sum(axis=1)


def calculate_curve(matrices, keep, names, target, edges, minimum_tuples):
    """Apply the guard to the intersection, never to a changing lag subset."""
    if keep.sum() < minimum_tuples:
        return []
    results = []
    for data in matrices:
        counts, _ = np.histogramdd(data[keep], bins=[edges[n] for n in [target] + names])
        if int(counts.sum()) != int(keep.sum()):
            raise ValueError('Histogram lost fixed support tuples')
        occupancy = dict(occupied_cells=int(np.count_nonzero(counts)),
                         singleton_cells=int(np.count_nonzero(counts == 1)),
                         histogram_cells=int(counts.size))
        if np.count_nonzero(counts.sum(axis=tuple(range(1, counts.ndim)))) < 2:
            results.append(dict(**occupancy, status='degenerate_target'))
            continue
        result = decompose_counts(counts, names)
        results.append(dict(**occupancy, **result,
                            status='ok' if result['normalization_defined'] else 'undefined_normalization'))
    return results


def plot_curve(key, group, out):
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    for ax, kind in zip(axes[:3], ['U', 'R', 'S']):
        for index, (component, curve) in enumerate(group[group.kind == kind].groupby('component')):
            color = plt.get_cmap('tab20')(index)
            for bins, style in [(2, ':'), (3, '-')]:
                selected = curve[curve.bins == bins]
                ax.plot(selected.lag_days, selected.fraction, linestyle=style, color=color,
                        label=component if bins == 3 else None)
        ax.set_ylabel(kind + ' / joint MI')
        ax.set_ylim(-.02, 1.02)
        ax.legend(fontsize=7, ncol=3)
    for bins, style in [(2, ':'), (3, '-')]:
        selected = group[group.bins == bins].drop_duplicates('lag_days')
        axes[3].plot(selected.lag_days, selected.normalized_leakage, linestyle=style,
                     label=f'{bins} bins', color='black')
    axes[3].legend()
    axes[3].set_ylim(-.02, 1.02)
    axes[3].set_ylabel('Normalized leakage')
    axes[3].set_xlabel('Base lag in observed days')
    fig.suptitle(' / '.join(map(str, key)) + '\nIdentical target dates across all lags; descriptive SURD')
    fig.tight_layout()
    fig.savefig(out / ('_'.join(map(str, key)) + '.png'), dpi=150)
    plt.close(fig)


def main():
    source_config, _, native_count = verify_inputs(BASE)
    lags = [lag for lag in source_config['lag_days'] if lag > 0]
    config = dict(source_config, bins=[2, 3], lag_days=lags,
                  support_rule='Intersect finite tuple masks across every declared positive lag within each curve',
                  bin_edge_rule='Full native campaign marginal quantiles, unchanged from saved scans',
                  zero_lag_rule='Excluded from the positive lag curve and support intersection',
                  uncertainty_policy='Flux errors retained in native inputs; no uncertainty propagation or null calibration',
                  interpretation='Descriptive estimation sensitivity; no physical delay or causal detection claim')
    OUT.mkdir(exist_ok=True)
    (OUT / 'configuration.json').write_text(json.dumps(config, indent=2) + '\n')
    audits, dates, summaries, atoms = [], [], [], []
    for year in sorted(source_config['campaign_jd_offsets']):
        series = {name: pd.read_csv(BASE / 'bins_2' / year / 'native' / f'{name}.csv')
                  for name in ['total', 'continuum', 'blue', 'core', 'red']}
        edges = {b: {n: np.r_[-np.inf, np.unique(np.quantile(f.value, np.arange(1, b) / b)), np.inf]
                     for n, f in series.items()} for b in (2, 3)}
        for subset, names in SUBSETS.items():
            for target in names:
                for scenario, offsets in source_config['unequal_offsets_days'].items():
                    key = dict(epoch=int(year), subset=subset, target=target, scenario=scenario)
                    matrices, keep, varying_counts = fixed_lag_tuples(
                        series, names, target, lags, offsets, source_config['maximum_interpolation_gap_days'])
                    count = int(keep.sum())
                    selected_dates = series[target].jd_offset.to_numpy()[keep]
                    date_hash = hashlib.sha256(np.asarray(selected_dates, dtype='<f8').tobytes()).hexdigest()
                    status = 'ok' if count >= source_config['minimum_native_targets'] else 'insufficient_support'
                    audits.append(dict(**key, native_target_count=len(keep), fixed_tuple_count=count,
                                       varying_tuple_minimum=int(varying_counts.min()),
                                       varying_tuple_maximum=int(varying_counts.max()),
                                       retained_fraction=count / len(keep),
                                       target_dates_sha256=date_hash, status=status))
                    dates.extend(dict(**key, jd_offset=float(t)) for t in selected_dates)
                    for bins in (2, 3):
                        results = calculate_curve(matrices, keep, names, target, edges[bins],
                                                  source_config['minimum_native_targets'])
                        for index, lag in enumerate(lags):
                            row = dict(**key, bins=bins, lag_days=lag, tuple_count=count,
                                       target_dates_sha256=date_hash, status=status)
                            if results:
                                result = results[index]
                                row.update({n: v for n, v in result.items() if n not in ['atoms', 'normalization_defined']})
                                if result['status'] == 'ok':
                                    for atom in result['atoms']:
                                        atoms.append(dict(**row, component=atom['kind'] + '(' + ','.join(atom['predictors']) + ')',
                                                          kind=atom['kind'], bits=atom['bits'], fraction=atom['fraction']))
                            summaries.append(row)
    audit = pd.DataFrame(audits)
    summary = pd.DataFrame(summaries)
    components = pd.DataFrame(atoms)
    audit.to_csv(OUT / 'curve_support.csv', index=False)
    pd.DataFrame(dates).to_csv(OUT / 'fixed_target_dates.csv', index=False)
    summary.to_csv(OUT / 'support_and_leakage.csv', index=False)
    components.to_csv(OUT / 'individual_components.csv', index=False)
    # Any degenerate target is explicitly reported. Only jointly defined estimates
    # enter a normalized comparison; no missing value is replaced with zero.
    jointly_defined = summary.query('bins == 2 and status == "ok"')[KEYS].merge(
        summary.query('bins == 3 and status == "ok"')[KEYS], on=KEYS, validate='one_to_one')
    if len(jointly_defined):
        matched_components = components.merge(jointly_defined, on=KEYS, validate='many_to_one')
        frames = [matched_components[matched_components.bins == b].drop(columns='bins') for b in (2, 3)]
        paired, changes = compare_atoms(frames, jointly_defined.assign(status='ok'))
        peaks = peak_review(paired)
        paired.to_csv(OUT / 'individual_atom_changes.csv', index=False)
        changes.to_csv(OUT / 'bin_comparison.csv', index=False)
        peaks.to_csv(OUT / 'component_peak_review.csv', index=False)
        review = changes.groupby(CURVES).agg(configurations=('lag_days', 'size'),
                                            median_total_variation=('total_variation', 'median'),
                                            maximum_total_variation=('total_variation', 'max'),
                                            median_absolute_leakage_change=('absolute_leakage_change', 'median')).reset_index()
        peak_summary = peaks.groupby(CURVES).agg(component_curves=('component', 'size'),
                                                fraction_peak_sets_overlap=('peak_sets_overlap', 'mean')).reset_index()
        review = review.merge(peak_summary, on=CURVES, validate='one_to_one')
        review = review.merge(audit[CURVES + ['fixed_tuple_count']], on=CURVES, validate='one_to_one')
        review.to_csv(OUT / 'curve_comparison.csv', index=False)
        support_changes = []
        for bins in (2, 3):
            reference = pd.read_csv(BASE / f'bins_{bins}' / 'individual_components.csv')
            fixed = matched_components[matched_components.bins == bins]
            comparison = fixed.merge(reference, on=KEYS + ['component', 'kind'],
                                     suffixes=('_fixed', '_varying'), validate='one_to_one')
            if len(comparison) != len(fixed):
                raise ValueError('Missing varying support reference atoms')
            comparison['absolute_change'] = abs(comparison.fraction_fixed - comparison.fraction_varying)
            change = (comparison.groupby(KEYS).absolute_change.sum() / 2).rename('total_variation').reset_index()
            leak = comparison.drop_duplicates(KEYS)[KEYS + ['normalized_leakage_fixed', 'normalized_leakage_varying']]
            change = change.merge(leak, on=KEYS, validate='one_to_one')
            change['absolute_leakage_change'] = abs(change.normalized_leakage_fixed - change.normalized_leakage_varying)
            change['bins'] = bins
            support_changes.append(change)
        support_change = pd.concat(support_changes, ignore_index=True)
        support_change.to_csv(OUT / 'support_change_comparison.csv', index=False)
        support_change.groupby(CURVES + ['bins']).agg(
            median_total_variation=('total_variation', 'median'),
            maximum_total_variation=('total_variation', 'max'),
            median_absolute_leakage_change=('absolute_leakage_change', 'median')).reset_index().to_csv(
                OUT / 'support_change_summary.csv', index=False)
        for key, group in matched_components.groupby(CURVES):
            plot_curve(key, group, OUT)
        print(review.to_string(index=False))
    inputs = sorted((BASE / 'bins_2').glob('*/native/*.csv')) + [
        BASE / f'bins_{b}' / name for b in (2, 3) for name in ('configuration.json', 'manifest.json', 'individual_components.csv')]
    manifest = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    configuration_sha256=hashlib.sha256((OUT / 'configuration.json').read_bytes()).hexdigest(),
                    official_source_sha256=HASHES, verified_native_files=native_count,
                    input_sha256={str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                    curve_count=len(audit), curves_passing_tuple_guard=int((audit.status == 'ok').sum()),
                    jointly_defined_comparison_rows=len(jointly_defined),
                    scan_rows=len(summary), usable_scan_rows=int((summary.status == 'ok').sum()), atom_rows=len(components),
                    status_counts=summary.status.value_counts().to_dict())
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(audit.groupby('epoch').fixed_tuple_count.agg(['min', 'max']).to_string())
    print(json.dumps({key: value for key, value in manifest.items() if key != 'input_sha256'}, indent=2))


if __name__ == '__main__':
    main()
