"""Explicit historical campaigns and equation 3.18 predictor construction."""
import numpy as np
from pipeline.sampling import interpolate

CAMPAIGNS = {1989: (47509, 47809), 1990: (47861, 48179),
             1991: (48225, 48534), 1992: (48623, 48898),
             1993: (48954, 49255)}


def campaign_labels(times):
    times = np.asarray(times, dtype=float)
    labels = np.zeros(times.shape, dtype=int)
    for year, (lo, hi) in CAMPAIGNS.items():
        labels[(times >= lo) & (times <= hi)] = year
    if (labels == 0).any():
        raise ValueError('Observation outside published campaigns')
    return labels


def expanded_predictors(series, names, target_times, base_lag, extra_delays,
                        maximum_gap=10):
    """Q_i(t), Q_i(t minus delay), with t = target time minus base lag.

    This is equation 3.18, not one unequal lag assigned to each variable.
    Positive extra delays are analysis choices, not measured physical delays.
    Input series must already be restricted to one campaign.
    """
    delays = tuple(float(d) for d in extra_delays)
    if base_lag < 0 or any(d <= 0 or not np.isfinite(d) for d in delays):
        raise ValueError('Base lag must be nonnegative and extra delays positive')
    if len(set(delays)) != len(delays):
        raise ValueError('Repeated extra delays')
    columns, labels = [], []
    for name in names:
        frame = series[name]
        for delay in (0.,) + delays:
            columns.append(interpolate(frame.jd_offset.to_numpy(), frame.value.to_numpy(),
                                       np.asarray(target_times) - base_lag - delay,
                                       maximum_gap))
            labels.append(f'{name}@{delay:g}')
    return np.column_stack(columns), labels
