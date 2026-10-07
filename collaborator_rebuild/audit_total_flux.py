"""Compare documented integration definitions without fitting a calibration."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from pipeline.adopted_data import corrected_table, load_adopted
from pipeline.reconcile_data import DATA, combine_dates
from collaborator_rebuild.lag_design import campaign_labels

BASE = Path(__file__).resolve().parent
INTERVALS = {'paper_section_3_2': (4870, 5010),
             'caption_alternative': (4840, 5010),
             'complete_file': (4840, 5162)}


def integrate_profile(array, lower, upper):
    """Sum pixel centres in the stated interval, retaining signed flux."""
    if not np.allclose(np.diff(array[:, 0]), 2):
        raise ValueError('Unexpected pixel spacing')
    selected = (array[:, 0] >= lower) & (array[:, 0] < upper)
    if not selected.any():
        raise ValueError('Empty integration interval')
    return float(array[selected, 1].sum()*2), float(np.sqrt((array[selected, 2]**2).sum())*2)


def main():
    out = BASE / 'support_review'
    out.mkdir(exist_ok=True)
    source = BASE / 'sources' / 'historical_revised_hbeta.dat'
    archive = pd.read_csv(source, sep=r'\s+', names=['date', 'flux', 'error'])
    archive['jd_offset'] = archive.date + 40000
    archive = archive.groupby('jd_offset').agg(archive_total=('flux', 'mean'),
                                               archive_measurements=('flux', 'size'))
    rows = []
    for observation in corrected_table().itertuples():
        a = np.loadtxt(DATA / 'hb_profiles_extracted' / observation.filename)
        row = dict(jd_offset=observation.table_date + 40000)
        for name, (lo, hi) in INTERVALS.items():
            row[name], row[name + '_error'] = integrate_profile(a, lo, hi)
        rows.append(row)
    profiles = combine_dates(pd.DataFrame(rows), list(INTERVALS))
    adopted, _ = load_adopted()
    matched = profiles.merge(archive, on='jd_offset', validate='one_to_one').merge(
        adopted[['jd_offset', 'profile_total']], on='jd_offset', validate='one_to_one')
    if len(matched) != len(adopted):
        raise ValueError('Incomplete matched sample')
    matched['epoch'] = campaign_labels(matched.jd_offset)
    summaries = []
    for epoch, frame in [('all', matched)] + list(matched.groupby('epoch')):
        for name in list(INTERVALS) + ['profile_total']:
            relative = (frame[name] - frame.archive_total) / frame.archive_total
            summaries.append(dict(epoch=epoch, definition=name, dates=len(frame),
                                  correlation=float(frame[name].corr(frame.archive_total)),
                                  median_relative_difference=float(relative.median()),
                                  maximum_absolute_relative_difference=float(relative.abs().max())))
    matched.to_csv(out / 'integration_comparison.csv', index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(out / 'integration_summary.csv', index=False)
    manifest = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    intervals_angstrom=INTERVALS, pixel_policy='Lower inclusive, upper exclusive',
                    status='Diagnostic only; no fitted flux correction or scientific input replacement',
                    caveats=['Archive duplicates cannot be linked to individual instruments by date alone',
                             'Profile units are arbitrary; similarity is not proof of absolute calibration',
                             'Source section and caption disagree on the lower boundary'])
    (out / 'integration_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(summary[summary.epoch.eq('all')].to_string(index=False))


if __name__ == '__main__':
    main()
