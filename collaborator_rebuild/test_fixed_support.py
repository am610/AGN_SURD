"""Verify date intersection, declared lag retention, and the fixed tuple guard."""
import numpy as np
import pandas as pd

from collaborator_rebuild.fixed_support import calculate_curve, fixed_lag_tuples


def test_same_target_dates_at_every_lag_without_extrapolation():
    times = np.array([0., 1., 2., 10., 11., 12.])
    frame = pd.DataFrame(dict(jd_offset=times, value=times + 100))
    series = {name: frame.copy() for name in ['total', 'blue', 'core', 'red']}
    names = list(series)
    matrices, keep, counts = fixed_lag_tuples(series, names, 'total', [1, 2], dict.fromkeys(names, 0), 2)
    np.testing.assert_array_equal(keep, [False, False, True, False, False, True])
    np.testing.assert_array_equal(counts, [4, 2])
    for lag, matrix in zip([1, 2], matrices):
        np.testing.assert_array_equal(matrix[keep, 0], [102, 112])
        np.testing.assert_array_equal(matrix[keep, 1], np.array([102, 112]) - lag)
        assert np.isfinite(matrix[keep]).all()


def test_unequal_offsets_participate_in_intersection():
    times = np.arange(6.)
    series = {name: pd.DataFrame(dict(jd_offset=times, value=times)) for name in ['total', 'blue', 'core', 'red']}
    names = list(series)
    offsets = dict.fromkeys(names, 0)
    offsets['blue'] = 1
    matrices, keep, counts = fixed_lag_tuples(series, names, 'total', [1, 2], offsets, 2)
    np.testing.assert_array_equal(keep, [False, False, False, True, True, True])
    np.testing.assert_array_equal(counts, [4, 3])
    np.testing.assert_array_equal(matrices[1][keep, 2], [0, 1, 2])


def test_guard_applies_to_intersection_not_individual_lags():
    # Each lag has three finite tuples, but only two target dates are shared.
    matrices = [np.ones((4, 5)), np.ones((4, 5))]
    matrices[0][0, 1] = np.nan
    matrices[1][1, 1] = np.nan
    keep = np.isfinite(np.stack(matrices)).all(axis=(0, 2))
    assert keep.sum() == 2
    assert calculate_curve(matrices, keep, ['a', 'b', 'c', 'd'], 'a', {}, minimum_tuples=3) == []


def test_degenerate_target_is_reported_without_fabricated_fractions():
    matrices = [np.ones((4, 5))]
    edges = {name: np.array([-np.inf, .5, np.inf]) for name in ['a', 'b', 'c', 'd']}
    results = calculate_curve(matrices, np.ones(4, dtype=bool), list(edges), 'a', edges, minimum_tuples=3)
    assert results[0]['status'] == 'degenerate_target'
    assert 'atoms' not in results[0]
