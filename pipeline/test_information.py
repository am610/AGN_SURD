import unittest

import numpy as np

from pipeline.information import decompose, predictors
from pipeline.benchmark_information import simulate
from pipeline.adopted_data import load_adopted
from pathlib import Path
import json


class InformationTests(unittest.TestCase):
    def test_independence(self):
        red, syn, joint, leak, entropy = decompose(np.ones((2, 2, 2)))
        self.assertAlmostEqual(joint, 0, places=9)
        self.assertAlmostEqual(leak, 1, places=9)

    def test_xor_is_synergy(self):
        hist = np.zeros((2, 2, 2))
        for a in range(2):
            for b in range(2):
                hist[a ^ b, a, b] = 1
        red, syn, joint, leak, entropy = decompose(hist)
        self.assertAlmostEqual(joint, 1, places=8)
        self.assertAlmostEqual(syn[(1, 2)], 1, places=8)
        self.assertAlmostEqual(sum(red.values()), 0, places=8)

    def test_duplicate_predictors_are_redundant(self):
        hist = np.zeros((2, 2, 2))
        hist[0, 0, 0] = hist[1, 1, 1] = 1
        red, syn, joint, leak, entropy = decompose(hist)
        self.assertAlmostEqual(red[(1, 2)], 1, places=8)
        self.assertAlmostEqual(sum(syn.values()), 0, places=8)

    def test_joint_mi_matches_direct_calculation(self):
        hist = np.random.default_rng(3).integers(0, 10, (3, 3, 3, 3))
        prob = hist / hist.sum()
        py = prob.sum(axis=(1, 2, 3), keepdims=True)
        px = prob.sum(axis=0, keepdims=True)
        mask = prob > 0
        expected = np.sum(prob[mask] * np.log2((prob / (py * px))[mask]))
        self.assertAlmostEqual(decompose(hist)[2], expected, places=8)

    def test_constant_target_is_not_a_detection(self):
        self.assertIsNone(decompose(np.array([[1, 1], [0, 0]])))

    def test_target_never_predicts_itself(self):
        for target in ('blue', 'core', 'red'):
            self.assertNotIn(target, predictors(target))
            self.assertEqual(len(predictors(target)), 3)

    def test_simulation_retains_times_errors_and_seed(self):
        lines, cont = load_adopted()
        config = json.loads(Path('pipeline/information_config.json').read_text())
        a = simulate(lines, cont, config, 'delayed_shared_driver', 12)
        b = simulate(lines, cont, config, 'delayed_shared_driver', 12)
        for source, first, second in zip([lines, cont], a, b):
            np.testing.assert_array_equal(source.jd_offset, first.jd_offset)
            np.testing.assert_array_equal(first, second)
            for name in source.columns:
                if name.endswith('_error'):
                    np.testing.assert_array_equal(source[name], first[name])

    def test_continuum_baseline_decomposition(self):
        from pipeline.continuum_baseline import decompose_1d
        # Perfect correlation
        hist = np.diag([10, 10])
        mi, leak, entropy = decompose_1d(hist)
        self.assertAlmostEqual(mi, 1.0, places=8)
        self.assertAlmostEqual(leak, 0.0, places=8)
        self.assertAlmostEqual(entropy, 1.0, places=8)
        # Complete independence
        hist = np.ones((2, 2)) * 10
        mi, leak, entropy = decompose_1d(hist)
        self.assertAlmostEqual(mi, 0.0, places=8)
        self.assertAlmostEqual(leak, 1.0, places=8)


if __name__ == '__main__':
    unittest.main()

