"""SURD scans with explicit predictors, histogram support and information units."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .sampling import tuples

UTILS = Path(__file__).resolve().parents[1] / 'SURD/utils'
sys.path.insert(0, str(UTILS))
import surd


def marginal_edges(values, bins):
    interior = np.unique(np.quantile(values, np.arange(1, bins) / bins))
    return np.r_[-np.inf, interior, np.inf]


def predictors(target):
    if target not in ('blue', 'core', 'red'):
        raise ValueError('Unknown line target')
    # The target component is excluded at every lag, keeping the estimand
    # identical across the scan and preventing contemporaneous self leakage.
    return ['continuum'] + [name for name in ('blue', 'core', 'red') if name != target]


def decompose(histogram):
    counts = np.asarray(histogram, dtype=float)
    if counts.sum() <= 0 or not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError('Invalid histogram')
    prob = counts / counts.sum()
    py = prob.sum(axis=tuple(range(1, prob.ndim)))
    entropy = float(-np.sum(py[py > 0] * np.log2(py[py > 0])))
    if entropy <= 1e-12:
        return None
    red, syn, mi, leak = surd.surd(prob.copy())
    joint = float(mi[tuple(range(1, counts.ndim))])
    total = float(sum(red.values()) + sum(syn.values()))
    if not np.isclose(total, joint, atol=1e-8):
        raise ValueError('SURD information accounting failed')
    return red, syn, joint, float(leak), entropy


def scan(lines, cont, sampling, config):
    edges = {name: marginal_edges(lines[name], config['bins']) for name in ('blue', 'core', 'red')}
    edges['continuum'] = marginal_edges(cont.continuum, config['bins'])
    rows = []
    for method in sampling['methods']:
        for lag in range(sampling['lag_min_days'], sampling['lag_max_days'] + 1, sampling['lag_step_days']):
            frame = tuples(lines, cont, lag, method, sampling)
            for target in config['targets']:
                names = predictors(target)
                columns = ['target_' + target] + ['continuum' if n == 'continuum' else 'predictor_' + n for n in names]
                row = dict(method=method, lag_days=lag, target=target,
                           predictors=';'.join(names), tuple_count=len(frame),
                           line_observed_dates=len(lines), continuum_observed_dates=len(cont),
                           status='insufficient_tuples')
                if len(frame) >= config['minimum_tuples']:
                    hist, _ = np.histogramdd(frame[columns].to_numpy(), bins=[edges[n] for n in [target] + names])
                    result = decompose(hist)
                    row.update(histogram_cells=hist.size, occupied_cells=np.count_nonzero(hist),
                               empty_fraction=float(np.mean(hist == 0)),
                               singleton_tuple_fraction=float(np.sum(hist == 1) / len(frame)))
                    if result is None:
                        row['status'] = 'constant_target'
                    else:
                        red, syn, joint, leak, entropy = result
                        synergy = float(sum(syn.values()))
                        row.update(status='ok', joint_mi_bits=joint, synergy_bits=synergy,
                                   unique_bits=float(sum(v for k, v in red.items() if len(k) == 1)),
                                   redundancy_bits=float(sum(v for k, v in red.items() if len(k) > 1)),
                                   target_entropy_bits=entropy, information_leak_fraction=leak,
                                   synergy_fraction=synergy / joint if joint > config['minimum_joint_mi_bits_for_ratio'] else np.nan)
                        for prefix, atoms in [('R', red), ('S', syn)]:
                            for key, value in atoms.items():
                                row[prefix + '_' + '_'.join(names[i - 1] for i in key) + '_bits'] = float(value)
                rows.append(row)
    return pd.DataFrame(rows)


def peaks(curves):
    rows = []
    for (method, target), group in curves.groupby(['method', 'target']):
        for statistic in ('joint_mi_bits', 'synergy_bits', 'synergy_fraction'):
            valid = group[group.status == 'ok'].dropna(subset=[statistic])
            if valid.empty:
                continue
            row = valid.loc[valid[statistic].idxmax()]
            rows.append(dict(method=method, target=target, statistic=statistic,
                             peak_lag_days=int(row.lag_days), peak_value=float(row[statistic]),
                             peak_tuple_count=int(row.tuple_count), scanned_lags=len(valid)))
    return pd.DataFrame(rows)
