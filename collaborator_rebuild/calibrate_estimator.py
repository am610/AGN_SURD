"""Calibrate the unchanged SURD estimator against exact categorical laws."""
import hashlib
import itertools
import json
import inspect
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from collaborator_rebuild.official_decomposition import HASHES, decompose_counts, official_module

BASE = Path(__file__).resolve().parent
OUT = BASE / 'published_results' / 'estimator_calibration'
NAMES = ['continuum', 'blue', 'core', 'red']
MODELS = ['independent', 'unique', 'redundant', 'synergistic']
CONFIG = dict(bins=[2, 3], sample_sizes=[36, 78, 90, 512], repetitions=200,
              seed=20261007, predictors=NAMES, target='future continuum',
              sample_size_rule='36, 78, and 90 match retained curves; 512 is an artificial reference',
              sampling='Independent draws from an exact categorical joint probability law',
              bin_rule='Known discrete states; no fitted quantile boundaries',
              model_definitions={'independent': 'All five variables mutually independent and uniform',
                                 'unique': 'Target equals blue; all four predictors independent and uniform',
                                 'redundant': 'Target equals blue equals core; continuum history and red independent',
                                 'synergistic': 'Target is blue plus core modulo bin count; predictors independent'},
              interpretation='Ideal finite sample calibration only; no AGN null test or delay recovery claim',
              limitations=['No time dependence, irregular observing dates, interpolation, or measurement errors',
                           'Two and three state laws are separate benchmarks, not rebinning of one continuous dataset',
                           'Perfect categorical signals do not represent realistic AGN effect sizes'])


def probability_law(model, bins):
    """Exact joint law with axes target, continuum history, blue, core, red."""
    if model not in MODELS or bins not in (2, 3):
        raise ValueError('Unknown model or bin count')
    law = np.zeros((bins,) * 5)
    if model == 'independent':
        law.fill(1 / law.size)
    elif model == 'redundant':
        for history, blue, red in itertools.product(range(bins), repeat=3):
            law[blue, history, blue, blue, red] = 1 / bins ** 3
    else:
        for history, blue, core, red in itertools.product(range(bins), repeat=4):
            target = blue if model == 'unique' else (blue + core) % bins
            law[target, history, blue, core, red] = 1 / bins ** 4
    if not np.isclose(law.sum(), 1):
        raise ValueError('Probability law does not sum to one')
    return law


def atom_label(atom):
    return atom['kind'] + '(' + ','.join(atom['predictors']) + ')'


def exact_reference(model, bins):
    result = decompose_counts(probability_law(model, bins), NAMES)
    expected_mi = 0. if model == 'independent' else np.log2(bins)
    expected_leak = 1. if model == 'independent' else 0.
    if not np.isclose(result['joint_mi_bits'], expected_mi, atol=1e-10):
        raise ValueError('Exact joint information disagrees with model truth')
    if not np.isclose(result['normalized_leakage'], expected_leak, atol=1e-10):
        raise ValueError('Exact leakage disagrees with model truth')
    expected = {'unique': 'U(blue)', 'redundant': 'R(blue,core)', 'synergistic': 'S(blue,core)'}.get(model)
    if expected is not None:
        atom = next(a for a in result['atoms'] if atom_label(a) == expected)
        if not np.isclose(atom['fraction'], 1., atol=1e-10):
            raise ValueError('Exact component disagrees with model truth')
    return result, expected


def permutation_audit(counts):
    """Reorder axes and labels together without changing the probability law."""
    def fractions(result):
        return {(a['kind'], tuple(sorted(a['predictors']))): a['fraction'] for a in result['atoms']}
    original = decompose_counts(counts, NAMES)
    baseline = fractions(original)
    rows = []
    for order in itertools.permutations(range(4)):
        names = [NAMES[i] for i in order]
        result = decompose_counts(counts.transpose((0,) + tuple(i + 1 for i in order)), names)
        changed = fractions(result)
        rows.append(dict(predictor_order=names,
                         total_variation_from_original=float(sum(abs(baseline[k] - changed[k]) for k in baseline) / 2),
                         joint_mi_change_bits=float(result['joint_mi_bits'] - original['joint_mi_bits'])))
    return rows


def trace_strict_resets(counts):
    """Observe the official routine's strict comparison without changing it."""
    routine = official_module().surd
    lines, start = inspect.getsourcelines(routine)
    reset_line = start + next(i for i, line in enumerate(lines) if 'I1[inds_] = 0' in line)
    observations = []
    def trace(frame, event, arg):
        if frame.f_code is routine.__code__ and event == 'line' and frame.f_lineno == reset_line:
            values = frame.f_locals
            for index in values['inds_']:
                value = float(values['I1'][index])
                if value > .1:
                    observations.append(dict(target_state=int(values['t']),
                                             lower_order=int(values['l']),
                                             reset_specific_information=value,
                                             comparison_maximum=float(values['Il1max']),
                                             difference=float(values['Il1max'] - value)))
        return trace
    prior = sys.gettrace()
    try:
        sys.settrace(trace)
        routine(np.asarray(counts, dtype=float).copy())
    finally:
        sys.settrace(prior)
    return observations


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'configuration.json').write_text(json.dumps(CONFIG, indent=2) + '\n')
    rng = np.random.default_rng(CONFIG['seed'])
    rows, references, witnesses = [], [], {}
    for model in MODELS:
        for bins in CONFIG['bins']:
            law = probability_law(model, bins)
            truth, expected = exact_reference(model, bins)
            references.append(dict(model=model, bins=bins, expected_component=expected,
                                   joint_mi_bits=truth['joint_mi_bits'], normalized_leakage=truth['normalized_leakage'],
                                   atoms=truth['atoms']))
            truth_bits = np.array([a['bits'] for a in truth['atoms']])
            for sample_size in CONFIG['sample_sizes']:
                for repetition in range(CONFIG['repetitions']):
                    counts = rng.multinomial(sample_size, law.ravel()).reshape(law.shape)
                    key = dict(model=model, bins=bins, sample_size=sample_size, repetition=repetition)
                    if np.count_nonzero(counts.sum(axis=(1, 2, 3, 4))) < 2:
                        rows.append(dict(**key, status='degenerate_target'))
                        continue
                    result = decompose_counts(counts, NAMES)
                    estimated_bits = np.array([a['bits'] for a in result['atoms']])
                    row = dict(**key, status='ok' if result['normalization_defined'] else 'undefined_normalization',
                               joint_mi_bits=result['joint_mi_bits'],
                               joint_mi_error_bits=result['joint_mi_bits'] - truth['joint_mi_bits'],
                               normalized_leakage=result['normalized_leakage'],
                               leakage_error=result['normalized_leakage'] - truth['normalized_leakage'],
                               absolute_atom_error_bits=float(abs(estimated_bits - truth_bits).sum()),
                               occupied_cells=int(np.count_nonzero(counts)),
                               singleton_cells=int(np.count_nonzero(counts == 1)))
                    if result['normalization_defined']:
                        fractions = np.array([a['fraction'] for a in result['atoms']])
                        if not np.isclose(fractions.sum(), 1, atol=1e-10):
                            raise ValueError('Sample atom fractions do not sum to one')
                        if expected is not None:
                            index = next(i for i, a in enumerate(result['atoms']) if atom_label(a) == expected)
                            row['expected_component_fraction'] = float(fractions[index])
                            row['fraction_total_variation'] = float(abs(fractions - np.array([a['fraction'] for a in truth['atoms']])).sum() / 2)
                            top = np.isclose(fractions, fractions.max(), atol=1e-12, rtol=0)
                            row['expected_component_is_largest'] = bool(top[index])
                            row['expected_component_is_sole_largest'] = bool(top[index] and top.sum() == 1)
                            witness_key = f'{model}_{bins}_{sample_size}'
                            if bins == 2 and sample_size in (90, 512) and not top[index] and witness_key not in witnesses:
                                witnesses[witness_key] = dict(**key, expected_component=expected,
                                                             counts=counts.astype(int).tolist(), result=result,
                                                             permutations=permutation_audit(counts),
                                                             strict_resets=trace_strict_resets(counts))
                    rows.append(row)
                print(f'Completed {model}, {bins} bins, {sample_size} tuples', flush=True)
    samples = pd.DataFrame(rows)
    samples.to_csv(OUT / 'repetitions.csv', index=False)
    (OUT / 'exact_references.json').write_text(json.dumps(references, indent=2) + '\n')
    summaries = []
    for key, group in samples.groupby(['model', 'bins', 'sample_size']):
        defined = group[group.status == 'ok']
        row = dict(zip(['model', 'bins', 'sample_size'], key), repetitions=len(group),
                   defined_repetitions=len(defined), median_joint_mi_bits=group.joint_mi_bits.median(),
                   q05_joint_mi_bits=group.joint_mi_bits.quantile(.05), q95_joint_mi_bits=group.joint_mi_bits.quantile(.95),
                   median_joint_mi_error_bits=group.joint_mi_error_bits.median(),
                   median_leakage=group.normalized_leakage.median(),
                   median_absolute_atom_error_bits=group.absolute_atom_error_bits.median(),
                   median_occupied_cells=group.occupied_cells.median(),
                   median_singleton_cells=group.singleton_cells.median())
        if key[0] != 'independent':
            row.update(median_fraction_total_variation=defined.fraction_total_variation.median(),
                       median_expected_component_fraction=defined.expected_component_fraction.median(),
                       q05_expected_component_fraction=defined.expected_component_fraction.quantile(.05),
                       q95_fraction_total_variation=defined.fraction_total_variation.quantile(.95),
                       fraction_expected_component_largest=defined.expected_component_is_largest.mean(),
                       fraction_expected_component_sole_largest=defined.expected_component_is_sole_largest.mean())
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / 'summary.csv', index=False)
    (OUT / 'allocation_witnesses.json').write_text(json.dumps(witnesses, indent=2) + '\n')
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for ax, model in zip(axes.ravel(), MODELS):
        for bins, color in [(2, 'tab:blue'), (3, 'tab:orange')]:
            g = summary[(summary.model == model) & (summary.bins == bins)].sort_values('sample_size')
            ax.plot(g.sample_size, g.median_joint_mi_bits, marker='o', color=color, label=f'{bins} states')
            ax.fill_between(g.sample_size, g.q05_joint_mi_bits, g.q95_joint_mi_bits, color=color, alpha=.15)
            ax.axhline(0 if model == 'independent' else np.log2(bins), color=color, linestyle=':')
        ax.set_title(model.capitalize())
        ax.set_xscale('log')
        ax.set_ylabel('Estimated joint information in bits')
        ax.set_xlabel('Independent tuples')
        ax.legend(fontsize=8)
    fig.suptitle('Ideal categorical calibration: median and central 90 percent of 200 draws\nDotted lines show exact truth; 512 tuples is an artificial reference')
    fig.tight_layout()
    fig.savefig(OUT / 'joint_information_calibration.png', dpi=180)
    plt.close(fig)
    files = [OUT / name for name in ['configuration.json', 'repetitions.csv', 'summary.csv', 'exact_references.json', 'allocation_witnesses.json', 'joint_information_calibration.png']]
    manifest = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    official_source_sha256=HASHES, sample_rows=len(samples), exact_laws=len(references),
                    status_counts=samples.status.value_counts().to_dict(),
                    output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
