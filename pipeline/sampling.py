"""Construct identical sampling rules for observations and future simulations."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from .adopted_data import ROOT, load_adopted
except ImportError:
    from adopted_data import ROOT, load_adopted

CONFIG = Path(__file__).with_name('sampling_config.json')


def interpolate(times, values, query, maximum_gap=np.inf):
    """No extrapolation; exact observations remain valid beside long gaps."""
    times, values, query = map(np.asarray, (times, values, query))
    right = np.searchsorted(times, query)
    safe = np.clip(right, 0, len(times) - 1)
    exact = times[safe] == query
    left = np.clip(right - 1, 0, len(times) - 1)
    valid = exact | ((right > 0) & (right < len(times)) &
                     ((times[safe] - times[left]) <= maximum_gap))
    result = np.interp(query, times, values)
    result[~valid] = np.nan
    return result


def nearest(times, query, tolerance):
    times, query = np.asarray(times), np.asarray(query)
    right = np.clip(np.searchsorted(times, query), 0, len(times) - 1)
    left = np.maximum(right - 1, 0)
    # Ties select the earlier observation deterministically.
    idx = np.where(abs(times[left] - query) <= abs(times[right] - query), left, right)
    valid = (abs(times[idx] - query) <= tolerance) & (query >= times[0]) & (query <= times[-1])
    return idx, valid


def tuples(lines, cont, lag, method, config):
    """Target line at t; continuum and all line predictors at t less lag.

    Keep every component so each target can use the same support. Predictor
    history and incremental nulls belong to the later inference stage.
    """
    lt, ct = lines.jd_offset.to_numpy(), cont.jd_offset.to_numpy()
    if method == 'native':
        t = lt.copy()
        query = t - lag
        li, lv = nearest(lt, query, config['native_pair_tolerance_days'])
        ci, cv = nearest(ct, query, config['native_pair_tolerance_days'])
        valid = lv & cv
        frame = pd.DataFrame({'target_time': t, 'predictor_time': query,
                              'line_predictor_observed_time': lt[li],
                              'continuum_predictor_observed_time': ct[ci],
                              'continuum': cont.continuum.to_numpy()[ci]})
        for name in ['blue', 'core', 'red']:
            frame['target_' + name] = lines[name].to_numpy()
            frame['predictor_' + name] = lines[name].to_numpy()[li]
        return frame[valid].reset_index(drop=True)
    if method not in {'daily', 'gap_limited', 'within_season'}:
        raise ValueError('Unknown sampling method')
    t = np.arange(np.ceil(lt[0]), np.floor(lt[-1]) + 1, config['grid_step_days'])
    query = t - lag
    gap = np.inf if method == 'daily' else config['maximum_interpolation_gap_days']
    frame = pd.DataFrame({'target_time': t, 'predictor_time': query,
                          'continuum': interpolate(ct, cont.continuum, query, gap)})
    for name in ['blue', 'core', 'red']:
        frame['target_' + name] = interpolate(lt, lines[name], t, gap)
        frame['predictor_' + name] = interpolate(lt, lines[name], query, gap)
    if method == 'within_season':
        # A gap in either required source ends the common observing season.
        # Midpoints only label seasons; interpolation support still enforces
        # the stricter maximum gap for every variable.
        edges = np.unique(np.concatenate([
            (times[:-1] + times[1:])[np.diff(times) > config['season_boundary_gap_days']] / 2
            for times in [lt, ct]]))
        valid = np.searchsorted(edges, t) == np.searchsorted(edges, query)
        frame = frame[valid]
    return frame.dropna().reset_index(drop=True)


def main():
    config = json.loads(CONFIG.read_text())
    lines, cont = load_adopted(config['continuum'])
    rows = []
    for lag in range(config['lag_min_days'], config['lag_max_days'] + 1, config['lag_step_days']):
        counts = {}
        for method in config['methods']:
            frame = tuples(lines, cont, lag, method, config)
            counts[method] = len(frame)
            rows.append({'method': method, 'lag_days': lag, 'tuple_count': len(frame),
                         'line_observed_dates': len(lines), 'continuum_observed_dates': len(cont),
                         'unique_target_times': frame.target_time.nunique(),
                         'unique_line_predictor_dates': frame.line_predictor_observed_time.nunique() if method == 'native' else None,
                         'unique_continuum_predictor_dates': frame.continuum_predictor_observed_time.nunique() if method == 'native' else None})
        assert counts['within_season'] <= counts['gap_limited'] <= counts['daily']
        assert counts['native'] <= len(lines)
    out = ROOT / 'agn_surd_project/processed/reconciled/sampling'
    out.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(out / 'support_counts.csv', index=False)
    manifest = {'status': 'Support audit only; no information statistic or significance calculated',
                'configuration': config,
                'configuration_sha256': hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
                'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'adopted_manifest_sha256': hashlib.sha256((out.parent / 'adopted/manifest.json').read_bytes()).hexdigest(),
                'native_reuse_policy': 'One tuple per target observation; predictor observations may be reused and are not independent samples',
                'lag_convention': 'Positive lag means target follows predictors',
                'error_policy': 'Flux errors stay in native input tables; sampling does not treat interpolated values as independent measurements'}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(result[result.lag_days.isin([0, 30, 100, 200])].to_string(index=False))


if __name__ == '__main__':
    main()
