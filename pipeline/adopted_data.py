"""Select observations with independently checked date and instrument provenance."""
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from .reconcile_data import ROOT, DATA, WINDOWS, combine_dates, load_native
except ImportError:
    from reconcile_data import ROOT, DATA, WINDOWS, combine_dates, load_native

OUT = ROOT / 'agn_surd_project/processed/reconciled/adopted'
TABLE = ROOT / 'audit/iw_bp_96_tab.txt'
REVISED = DATA / 'c5150_revised.dat'


def corrected_table():
    text = TABLE.read_text().split('Table 2.')[0]
    rows = [(int(m[1]), m[2], float(m[3])) for m in re.finditer(
        r'(?m)^(\d{4})\s+([a-z]+)\s+(\d+\.\d+)', text)]
    frame = pd.DataFrame(rows, columns=['table_date', 'instrument', 'table_continuum'])
    if len(frame) != 242 or frame.duplicated(['table_date', 'instrument']).any():
        raise ValueError('Unexpected corrected table structure')
    frame['filename'] = [f'n5{t:04d}{i}.spc' for t, i, _ in rows]
    return frame


def load_adopted(continuum='original'):
    spectra, _, original, paths = load_native()
    table = corrected_table()
    selected = spectra.merge(table, on='filename', validate='one_to_one')
    if len(selected) != len(table):
        raise ValueError('A corrected table spectrum is missing')
    np.testing.assert_array_equal(selected.jd_offset, selected.table_date + 40000)
    lines = combine_dates(selected, [*WINDOWS, 'profile_total'])
    if continuum == 'original':
        cont = original
    elif continuum == 'revised':
        raw, _ = reconcile_revised(table)
        cont = combine_dates(raw[raw.table_matched], ['continuum'])
    else:
        raise ValueError('Unknown continuum choice')
    # Keep actual endpoints. Pairing must enforce each series support at the
    # requested lag instead of repeatedly shrinking the observation window.
    return lines.copy(), cont.copy()


def reconcile_revised(table):
    raw = pd.read_csv(REVISED, sep=r'\s+', names=['table_date', 'continuum', 'continuum_error'])
    if not np.isfinite(raw.to_numpy()).all() or (raw.continuum_error < 0).any():
        raise ValueError('Invalid revised continuum')
    raw['jd_offset'] = raw.table_date + 40000
    raw['table_matched'] = False
    raw['filename'] = ''
    matches = []
    # Table values have three decimals. Allow the observed 0.000504 rounding
    # discrepancy, with a fixed tolerance chosen from printed precision.
    for row in table.itertuples():
        delta = abs(raw.continuum - row.table_continuum)
        hits = raw.index[(raw.table_date == row.table_date) &
                         (delta <= 0.00051) & ~raw.table_matched]
        if len(hits) != 1:
            raise ValueError(f'Nonunique continuum match: {row.filename}')
        idx = hits[0]
        raw.loc[idx, ['table_matched', 'filename']] = [True, row.filename]
        matches.append({'filename': row.filename, 'jd_offset': raw.loc[idx, 'jd_offset'],
                        'absolute_flux_difference': delta.loc[idx]})
    return raw, pd.DataFrame(matches)


def main():
    table = corrected_table()
    spectra, _, _, paths = load_native()
    revised, matches = reconcile_revised(table)
    excluded = spectra[~spectra.filename.isin(table.filename)].copy()
    np.testing.assert_array_equal(np.sort(excluded.jd_offset),
                                  np.sort(revised.loc[~revised.table_matched, 'jd_offset']))
    if len(excluded) != 5:
        raise ValueError('Unexpected exclusion count')
    lines, cont = load_adopted()
    revised_lines, revised_cont = load_adopted('revised')
    np.testing.assert_array_equal(lines.jd_offset, revised_lines.jd_offset)
    np.testing.assert_array_equal(lines.jd_offset, revised_cont.jd_offset)
    np.testing.assert_allclose(lines.profile_total, lines[list(WINDOWS)].sum(axis=1))
    OUT.mkdir(parents=True, exist_ok=True)
    for name, frame in [('line_native', lines), ('continuum_native', cont),
                        ('continuum_revised_native', revised_cont),
                        ('excluded_spectra', excluded), ('revised_continuum_audit', revised),
                        ('continuum_table_matches', matches)]:
        frame.to_csv(OUT / f'{name}.csv', index=False)
    summary = []
    for name, frame in [('line', lines), ('continuum_original', cont), ('continuum_revised', revised_cont)]:
        summary.append(dict(series=name, native_dates=len(frame),
                            source_observations=int(frame.n_observations.sum()),
                            first_jd_offset=float(frame.jd_offset.min()),
                            last_jd_offset=float(frame.jd_offset.max()),
                            largest_gap_days=float(np.diff(frame.jd_offset).max())))
    pd.DataFrame(summary).to_csv(OUT / 'data_table.csv', index=False)
    inputs = paths + [DATA / 'c5100.dat', REVISED, TABLE,
                      ROOT / 'audit/wanders_peterson_1996.pdf',
                      Path(__file__).resolve(), ROOT / 'pipeline/reconcile_data.py']
    manifest = {
        'status': 'Adopted sample under explicit exclusion and continuity policies',
        'primary_continuum': 'Original 5100 archive series; preserves prior analysis and denser observed sampling',
        'sensitivity_continuum': 'Revised archive series labelled 5150; paper section 3.1 describes 5188 Angstrom observed',
        'continuum_units': {'original': '1e-15 erg / s / cm^2 / Angstrom',
                            'revised': '1e-14 erg / s / cm^2 / Angstrom'},
        'time_coordinate': 'JD less 2400000',
        'revised_time_conversion': 'Add 40000 to file dates; validated against all 242 corrected table rows',
        'revised_flux_match_tolerance': 0.00051,
        'revised_flux_match_max_difference': float(matches.absolute_flux_difference.max()),
        'spectral_date_precision': 'Integer JD only; fractional exposure times unavailable',
        'time_scale_limit': 'UTC, TT and heliocentric correction not established',
        'exclusion_policy': 'Exclude every spectrum without an exact date and instrument match in corrected Table 1',
        'excluded_spectra': excluded.filename.tolist(),
        'wavelength_frame': 'Observed wavelength, supported by source paper section 3.2',
        'velocity_km_s': WINDOWS,
        'reference_wavelength_angstrom': 4861.33,
        'redshift_adopted': 0.017175,
        'redshift_policy': 'Numerical reference retained from prior reduction; not the exact value used in the source paper',
        'source_paper_redshifts': {'broad_peak': 0.0165, 'narrow_lines': 0.0174, 'radio_21cm': 0.0172},
        'source_paper_section_3_2_intervals_angstrom': [[4870, 4920], [4920, 4960], [4960, 5010]],
        'source_paper_caption_discrepancy': 'Figure 1 caption prints 4840 for lower boundary; section 3.2 prints 4870',
        'line_units': 'Arbitrary flux density units times Angstrom',
        'integration': 'Sum pixels times 2 Angstrom; lower inclusive and upper exclusive velocity boundaries',
        'duplicate_policy': 'Equal arithmetic means; errors propagated in quadrature',
        'covariance_limit': 'Shared calibration and pixel covariance unavailable',
        'subtraction_conventions': 'Profiles remove continuum and narrow lines; primary continuum retains host starlight; revised continuum subtracts host contribution with aperture residuals discussed in section 3.1',
        'sources': ['https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/',
                    'https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/iw-ct.lcv',
                    'https://www.asc.ohio-state.edu/astronomy/agnwatch/papers/iw_bp_96_tab.ps.gz',
                    'https://articles.adsabs.harvard.edu/pdf/1996ApJ...466..174W'],
        'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'versions': {'numpy': np.__version__, 'pandas': pd.__version__},
    }
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(summary).to_string(index=False))
    print('All corrected table dates and revised continuum matches passed')


if __name__ == '__main__':
    main()
