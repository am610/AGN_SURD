import numpy as np
import pandas as pd
import pytest
from collaborator_rebuild.lag_design import campaign_labels, expanded_predictors


def test_published_campaigns():
    assert campaign_labels([47509, 47809, 47861, 49255]).tolist() == [1989, 1989, 1990, 1993]
    with pytest.raises(ValueError):
        campaign_labels([48000, 48599])


def test_equation_318_keeps_both_copies():
    series = {'a': pd.DataFrame({'jd_offset': np.arange(20.), 'value': np.arange(20.)}),
              'b': pd.DataFrame({'jd_offset': np.arange(20.), 'value': 2*np.arange(20.)})}
    x, labels = expanded_predictors(series, ['a', 'b'], [10., 12.], 2, [5])
    assert labels == ['a@0', 'a@5', 'b@0', 'b@5']
    np.testing.assert_allclose(x, [[8, 3, 16, 6], [10, 5, 20, 10]])


def test_no_extrapolation():
    series = {'a': pd.DataFrame({'jd_offset': [0., 1., 20.], 'value': [0., 1., 20.]})}
    x, _ = expanded_predictors(series, ['a'], [-1., 10.], 0, [])
    assert np.isnan(x).all()
