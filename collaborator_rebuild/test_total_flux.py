import numpy as np
from collaborator_rebuild.audit_total_flux import integrate_profile
from collaborator_rebuild.prepare_historical import published_lines


def test_interval_centres_and_signed_flux():
    a = np.array([[4870., 2., .1], [4872., -1., .2], [4874., 50., .3]])
    flux, error = integrate_profile(a, 4870, 4874)
    assert flux == 2.
    np.testing.assert_allclose(error, 2*np.sqrt(.05))


def test_corrected_published_flux_columns():
    frame = published_lines()
    first = frame.iloc[0]
    assert first.total == 67.54
    assert first.blue == 19.04
    assert first.core == 30.30
    assert first.red == 18.22
    assert frame.jd_offset.nunique() == 224
    assert len(frame) == 242
