"""Compare atom distributions on exactly matched descriptive configurations."""
import json
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent / 'published_results'
KEYS = ['epoch', 'subset', 'target', 'scenario', 'lag_days']


def main():
    frames = [pd.read_csv(BASE / f'bins_{bins}' / 'individual_components.csv') for bins in [2, 3]]
    paired = frames[0].merge(frames[1], on=KEYS + ['component', 'kind'],
                             suffixes=('_2', '_3'), validate='one_to_one')
    paired = paired[paired.lag_days > 0].copy()
    if paired[['fraction_2', 'fraction_3']].isna().any().any():
        raise ValueError('Undefined normalization in matched comparison')
    paired['absolute_change'] = abs(paired.fraction_2 - paired.fraction_3)
    change = (paired.groupby(KEYS).absolute_change.sum()/2).rename('total_variation').reset_index()
    leakage = paired.drop_duplicates(KEYS)[KEYS + ['normalized_leakage_2', 'normalized_leakage_3']]
    change = change.merge(leakage, on=KEYS, validate='one_to_one')
    change['absolute_leakage_change'] = abs(change.normalized_leakage_2 - change.normalized_leakage_3)
    change.to_csv(BASE / 'bin_comparison.csv', index=False)
    summary = dict(positive_lag_configurations=len(change),
                   median_total_variation=float(change.total_variation.median()),
                   maximum_total_variation=float(change.total_variation.max()),
                   fraction_total_variation_above_point_two=float((change.total_variation > .2).mean()),
                   median_absolute_leakage_change=float(change.absolute_leakage_change.median()),
                   interpretation='Descriptive bin sensitivity, not a significance test')
    (BASE / 'bin_comparison.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
