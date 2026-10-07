"""Checks for missing atoms, invalid fractions, and ambiguous peak lags."""
import numpy as np
import pandas as pd
import pytest

from collaborator_rebuild.review_campaigns import compare_atoms, peak_review


def example():
    rows = [dict(epoch=1993, subset='total_wings', target='total', scenario='common',
                 lag_days=lag, component=f'atom_{atom}', kind='U', fraction=1 / 26,
                 normalized_leakage=.2) for lag in (1, 2) for atom in range(26)]
    frame = pd.DataFrame(rows)
    support = frame.drop_duplicates('lag_days').assign(status='ok')
    return frame, support


def test_missing_atom_is_rejected():
    frame, support = example()
    with pytest.raises(ValueError, match='inventories differ'):
        compare_atoms([frame, frame.iloc[1:]], support)


def test_invalid_normalization_is_rejected():
    frame, support = example()
    bad = frame.copy()
    bad.loc[0, 'fraction'] = np.nan
    with pytest.raises(ValueError, match='Undefined'):
        compare_atoms([frame, bad], support)
    bad.loc[0, 'fraction'] = .5
    with pytest.raises(ValueError, match='sum to one'):
        compare_atoms([frame, bad], support)


def test_tied_peaks_preserve_all_lags():
    frame, support = example()
    paired, changes = compare_atoms([frame, frame.copy()], support)
    assert (changes.total_variation == 0).all()
    peaks = peak_review(paired)
    assert (peaks.peak_lags_2 == '[1, 2]').all()
    assert peaks.peak_sets_overlap.all()
    assert (peaks.minimum_peak_separation_days == 0).all()
