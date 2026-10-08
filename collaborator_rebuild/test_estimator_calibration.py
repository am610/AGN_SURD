"""Verify exact laws before finite sample calibration."""
import numpy as np
import pytest

from collaborator_rebuild.calibrate_estimator import MODELS, exact_reference, probability_law


@pytest.mark.parametrize('bins', [2, 3])
@pytest.mark.parametrize('model', MODELS)
def test_exact_laws_match_known_information_and_atom(model, bins):
    law = probability_law(model, bins)
    assert law.shape == (bins,) * 5
    assert np.isclose(law.sum(), 1)
    np.testing.assert_allclose(law.sum(axis=(1, 2, 3, 4)), np.full(bins, 1 / bins))
    truth, expected = exact_reference(model, bins)
    assert len(truth['atoms']) == 26
    assert truth['normalization_defined'] == (model != 'independent')
    if model == 'independent':
        assert expected is None
        assert all(atom['fraction'] is None for atom in truth['atoms'])


def test_redundant_law_keeps_duplicate_predictors():
    law = probability_law('redundant', 3)
    for y, history, blue, core, red in np.argwhere(law > 0):
        assert y == blue == core


def test_modular_synergy_is_uninformative_from_either_driver_alone():
    law = probability_law('synergistic', 3)
    np.testing.assert_allclose(law.sum(axis=(1, 3, 4)), np.full((3, 3), 1 / 9))
    np.testing.assert_allclose(law.sum(axis=(1, 2, 4)), np.full((3, 3), 1 / 9))
