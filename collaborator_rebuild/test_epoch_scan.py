import numpy as np
import pandas as pd
import pytest

from collaborator_rebuild.epoch_scan import CONFIGS, epoch_grid, lag_tuples, quantile_edges
from collaborator_rebuild.official_decomposition import decompose_counts


def frame(times, values):
    return pd.DataFrame({'jd_offset': times, 'value': values})


def toy_series():
    t = np.arange(0., 60., 2.)
    return {n: frame(t, np.sin(t / (6. + k)) + k) for k, n in enumerate(['continuum', 'total', 'blue', 'core', 'red'])}


def test_grid_spans_common_support():
    series = toy_series()
    series['red'] = frame(np.arange(10., 50., 2.), np.zeros(20))
    grid = epoch_grid(series)
    assert grid[0] == 10 and grid[-1] == 48


def test_no_extrapolation_and_gap_limit():
    series = toy_series()
    series['blue'] = frame(np.r_[0., 2., 40., 42.], np.ones(4))
    names = ['continuum', 'blue', 'core', 'red']
    data, dates = lag_tuples(series, names, 'core', epoch_grid(series), 3, CONFIGS['common'], 15)
    assert np.isfinite(data).all()
    blue_times = dates - 3
    assert ((blue_times <= 2) | (blue_times >= 40)).all()
    assert (blue_times >= 0).all()


def test_offset_moves_only_that_predictor():
    series = {n: frame(np.arange(0., 80.), np.arange(0., 80.) + k) for k, n in
              enumerate(['continuum', 'total', 'blue', 'core', 'red'])}
    names = ['continuum', 'blue', 'core', 'red']
    data, dates = lag_tuples(series, names, 'core', epoch_grid(series), 4, CONFIGS['blue_leading'], 15)
    # columns: target core, continuum, blue, core, red; values equal time plus series index
    assert np.allclose(data[:, 0], dates + 3)
    assert np.allclose(data[:, 1], dates - 4)
    assert np.allclose(data[:, 2], dates - 9 + 2)
    assert np.allclose(data[:, 3], dates - 4 + 3)
    assert np.allclose(data[:, 4], dates - 4 + 4)


def test_scan_row_has_26_normalized_atoms():
    rng = np.random.default_rng(3)
    series = {n: frame(np.arange(0., 120.), rng.normal(size=120)) for n in
              ['continuum', 'total', 'blue', 'core', 'red']}
    names = ['total', 'blue', 'core', 'red']
    data, _ = lag_tuples(series, names, 'total', epoch_grid(series), 5, CONFIGS['common'], 15)
    edges = quantile_edges(series, 2)
    counts, _ = np.histogramdd(data, bins=[edges[n] for n in ['total'] + names])
    result = decompose_counts(counts, names)
    assert len(result['atoms']) == 26
    assert sum(a['fraction'] for a in result['atoms']) == pytest.approx(1)
