"""Audit every predictor order on the four retained astronomical curves."""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from collaborator_rebuild.fixed_support import SUBSETS, fixed_lag_tuples
from collaborator_rebuild.official_decomposition import HASHES, decompose_counts
from collaborator_rebuild.review_campaigns import BASE, CURVES, verify_inputs

OUT = BASE / 'predictor_order_audit'
SNAPSHOT = Path(__file__).resolve().parent / 'order_review'
TOLERANCE = 1e-8


def canonical_atoms(result, field='fraction'):
    return {(a['kind'], tuple(sorted(a['predictors']))): a[field] for a in result['atoms']}


def compare_named_atoms(original, reordered):
    first, second = canonical_atoms(original), canonical_atoms(reordered)
    if first.keys() != second.keys():
        raise ValueError('Named component inventories differ')
    differences = np.array([abs(first[key] - second[key]) for key in first])
    return float(differences.sum() / 2), float(differences.max())


def main():
    config, _, _ = verify_inputs(BASE)
    fixed = BASE / 'fixed_support'
    settings = json.loads((fixed / 'configuration.json').read_text())
    manifest = json.loads((fixed / 'manifest.json').read_text())
    for relative, expected in manifest['input_sha256'].items():
        if hashlib.sha256((BASE / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Fixed scan input changed: {relative}')
    components = pd.read_csv(fixed / 'individual_components.csv')
    saved_support = pd.read_csv(fixed / 'support_and_leakage.csv')
    support = pd.read_csv(fixed / 'curve_support.csv').query('status == "ok"')
    if len(support) != 4:
        raise ValueError('Audit is scoped to the four retained curves')
    OUT.mkdir(exist_ok=True)
    audit_config = dict(curves=support[CURVES].to_dict('records'), bins=settings['bins'],
                        lag_days=settings['lag_days'], predictor_orders=24,
                        comparison_tolerance=TOLERANCE,
                        tolerance_scope='Numerical comparison only, not a significance threshold',
                        invariant='Reorder predictor histogram axes and their names together; target axis stays fixed',
                        analysis_boundary='No new data, bin settings, lag ranges, offsets, null tests, or algorithm changes')
    (OUT / 'configuration.json').write_text(json.dumps(audit_config, indent=2) + '\n')
    rows, configurations, witnesses = [], [], {}
    for curve in support.itertuples():
        key = {name: getattr(curve, name) for name in CURVES}
        year, names = str(curve.epoch), SUBSETS[curve.subset]
        series = {name: pd.read_csv(BASE / 'bins_2' / year / 'native' / f'{name}.csv') for name in names}
        matrices, keep, _ = fixed_lag_tuples(series, names, curve.target, settings['lag_days'],
                                           config['unequal_offsets_days'][curve.scenario],
                                           config['maximum_interpolation_gap_days'])
        dates = np.asarray(series[curve.target].jd_offset.to_numpy()[keep], dtype='<f8')
        date_hash = hashlib.sha256(dates.tobytes()).hexdigest()
        if int(keep.sum()) != curve.fixed_tuple_count or date_hash != curve.target_dates_sha256:
            raise ValueError('Reconstructed target dates differ from fixed scan')
        for bins in settings['bins']:
            edges = {n: np.r_[-np.inf, np.unique(np.quantile(f.value, np.arange(1, bins) / bins)), np.inf]
                     for n, f in series.items()}
            for lag, data in zip(settings['lag_days'], matrices):
                case = dict(**key, bins=bins, lag_days=lag)
                counts, _ = np.histogramdd(data[keep], bins=[edges[n] for n in [curve.target] + names])
                original = decompose_counts(counts, names)
                mask = (components.epoch == curve.epoch) & (components.subset == curve.subset) & (
                    components.target == curve.target) & (components.scenario == curve.scenario) & (
                    components.bins == bins) & (components.lag_days == lag)
                reference = components[mask]
                if len(reference) != 26 or counts.sum() != curve.fixed_tuple_count:
                    raise ValueError('Reconstructed configuration inventory differs')
                actual = {a['kind'] + '(' + ','.join(a['predictors']) + ')': a['fraction'] for a in original['atoms']}
                reference_error = max(abs(actual[r.component] - r.fraction) for r in reference.itertuples())
                if reference_error > 1e-10:
                    raise ValueError('Original order fails to reproduce saved components')
                support_mask = (saved_support.epoch == curve.epoch) & (saved_support.subset == curve.subset) & (
                    saved_support.target == curve.target) & (saved_support.scenario == curve.scenario) & (
                    saved_support.bins == bins) & (saved_support.lag_days == lag)
                saved = saved_support[support_mask].iloc[0]
                if not np.isclose(original['joint_mi_bits'], saved.joint_mi_bits, atol=1e-10, rtol=0):
                    raise ValueError('Original joint information does not reproduce')
                order_rows = []
                for order in itertools.permutations(range(4)):
                    reordered_names = [names[i] for i in order]
                    result = decompose_counts(counts.transpose((0,) + tuple(i + 1 for i in order)), reordered_names)
                    if not result['normalization_defined']:
                        raise ValueError('Reordered normalization became undefined')
                    fractions = np.array([a['fraction'] for a in result['atoms']])
                    if not np.isclose(fractions.sum(), 1, atol=1e-10, rtol=0):
                        raise ValueError('Reordered components do not normalize')
                    tv, maximum = compare_named_atoms(original, result)
                    row = dict(**case, predictor_order=','.join(reordered_names), total_variation=tv,
                               maximum_absolute_component_change=maximum,
                               absolute_joint_mi_change_bits=abs(result['joint_mi_bits'] - original['joint_mi_bits']),
                               absolute_leakage_change=abs(result['normalized_leakage'] - original['normalized_leakage']))
                    rows.append(row)
                    order_rows.append(row)
                    witness_key = '_'.join(map(str, [curve.epoch, curve.scenario, bins]))
                    if witness_key not in witnesses or tv > witnesses[witness_key]['total_variation']:
                        witnesses[witness_key] = dict(**case, total_variation=tv, predictor_order=reordered_names,
                                                      original_predictor_order=names, counts=counts.astype(int).tolist(),
                                                      target_dates_sha256=date_hash, original=original, reordered=result)
                configs = pd.DataFrame(order_rows)
                configurations.append(dict(**case, tuple_count=int(keep.sum()), target_dates_sha256=date_hash,
                                           reference_maximum_error=reference_error,
                                           maximum_total_variation=configs.total_variation.max(),
                                           maximum_absolute_component_change=configs.maximum_absolute_component_change.max(),
                                           changed_orders=int((configs.total_variation > TOLERANCE).sum()),
                                           maximum_joint_mi_change_bits=configs.absolute_joint_mi_change_bits.max(),
                                           maximum_leakage_change=configs.absolute_leakage_change.max()))
            print(f'Completed {curve.epoch}, {curve.scenario}, {bins} bins', flush=True)
    result = pd.DataFrame(configurations)
    pd.DataFrame(rows).to_csv(OUT / 'all_orders.csv', index=False)
    result.to_csv(OUT / 'configuration_summary.csv', index=False)
    result['affected'] = result.maximum_total_variation > TOLERANCE
    summary = result.groupby(CURVES + ['bins']).agg(
        lag_configurations=('lag_days', 'size'), affected_configurations=('affected', 'sum'),
        median_maximum_total_variation=('maximum_total_variation', 'median'),
        maximum_total_variation=('maximum_total_variation', 'max'),
        maximum_absolute_component_change=('maximum_absolute_component_change', 'max'),
        maximum_joint_mi_change_bits=('maximum_joint_mi_change_bits', 'max'),
        maximum_leakage_change=('maximum_leakage_change', 'max')).reset_index()
    summary.to_csv(OUT / 'curve_summary.csv', index=False)
    (OUT / 'witnesses.json').write_text(json.dumps(witnesses, indent=2) + '\n')
    inputs = [fixed / name for name in ['configuration.json', 'manifest.json', 'individual_components.csv',
                                      'support_and_leakage.csv', 'curve_support.csv']]
    outputs = [OUT / name for name in ['configuration.json', 'all_orders.csv', 'configuration_summary.csv',
                                      'curve_summary.csv', 'witnesses.json']]
    audit_manifest = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                          official_source_sha256=HASHES, configurations=len(result), order_evaluations=len(rows),
                          affected_configurations=int(result.affected.sum()),
                          comparison_tolerance=TOLERANCE,
                          maximum_total_variation=float(result.maximum_total_variation.max()),
                          input_sha256={str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                          output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs})
    (OUT / 'manifest.json').write_text(json.dumps(audit_manifest, indent=2) + '\n')
    SNAPSHOT.mkdir(exist_ok=True)
    for name in ['configuration.json', 'configuration_summary.csv', 'curve_summary.csv', 'witnesses.json', 'manifest.json']:
        (SNAPSHOT / name).write_bytes((OUT / name).read_bytes())
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
