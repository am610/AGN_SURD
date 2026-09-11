"""Boundary checks for sampling behavior that can bias lag inference."""
import unittest

import numpy as np
import pandas as pd

from pipeline.sampling import interpolate, nearest, tuples


class SamplingTests(unittest.TestCase):
    def test_exact_measurements_survive_large_gap(self):
        result = interpolate([0, 2, 100], [0, 2, 100], [-1, 0, 1, 2, 50, 100, 101], 3)
        np.testing.assert_allclose(result, [np.nan, 0, 1, 2, np.nan, 100, np.nan], equal_nan=True)

    def test_nearest_ties_and_no_extrapolation(self):
        idx, valid = nearest([0, 2], [-0.5, 1, 2, 2.5], 1)
        np.testing.assert_array_equal(idx, [0, 0, 1, 1])
        np.testing.assert_array_equal(valid, [False, True, True, False])

    def test_season_crossing_removed_even_with_exact_endpoints(self):
        t = np.array([0., 1., 100., 101.])
        lines = pd.DataFrame({'jd_offset': t, 'blue': t, 'core': t, 'red': t})
        cont = pd.DataFrame({'jd_offset': t, 'continuum': t})
        config = {'grid_step_days': 1, 'maximum_interpolation_gap_days': 30,
                  'season_boundary_gap_days': 60, 'native_pair_tolerance_days': 1}
        self.assertEqual(len(tuples(lines, cont, 100, 'gap_limited', config)), 2)
        self.assertEqual(len(tuples(lines, cont, 100, 'within_season', config)), 0)
        native = tuples(lines, cont, 100, 'native', config)
        np.testing.assert_array_equal(native.target_time - native.line_predictor_observed_time, [100, 100])


if __name__ == '__main__':
    unittest.main()
