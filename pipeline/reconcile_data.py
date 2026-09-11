"""Rebuild native observations with explicit units and provenance."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'agn_surd_project/agn_data/ngc5548_agnwatch'
OUT = ROOT / 'agn_surd_project/processed/reconciled'
WINDOWS = {'blue': (-6000, -2000), 'core': (-2000, 2000), 'red': (2000, 6000)}

def combine_dates(frame, names):
    """Use equal weights so the total remains the sum of its components."""
    rows = []
    for time, group in frame.groupby('jd_offset', sort=True):
        row = {'jd_offset': time, 'n_observations': len(group)}
        for name in names:
            row[name] = float(group[name].mean())
            row[name + '_error'] = float(np.sqrt(np.sum(group[name + '_error'] ** 2)) / len(group))
        rows.append(row)
    return pd.DataFrame(rows)

def load_native():
    paths = sorted((DATA / 'hb_profiles_extracted').glob('*.spc'))
    if not paths:
        raise ValueError('No spectra found')
    rows = []
    for path in paths:
        a = np.loadtxt(path)
        if a.ndim != 2 or a.shape[1] != 3 or not np.isfinite(a).all() or np.any(a[:, 2] < 0):
            raise ValueError(f'Invalid spectrum: {path}')
        spacing = np.diff(a[:, 0])
        if not np.allclose(spacing, 2):
            raise ValueError(f'Unexpected wavelength grid: {path}')
        v = (a[:, 0] / (4861.33 * 1.017175) - 1) * 299792.458
        row = {'filename': path.name, 'jd_offset': float(re.fullmatch(r'n(\d{5})[a-z]+\.spc', path.name)[1]) - 10000}
        for name, (low, high) in WINDOWS.items():
            selected = (v >= low) & (v < high)
            if not selected.any():
                raise ValueError(f'Empty velocity interval: {path}')
            row[name] = float(np.sum(a[selected, 1] * 2))
            row[name + '_error'] = float(np.sqrt(np.sum((a[selected, 2] * 2) ** 2)))
        row['profile_total'] = sum(row[n] for n in WINDOWS)
        row['profile_total_error'] = float(np.sqrt(sum(row[n + '_error'] ** 2 for n in WINDOWS)))
        rows.append(row)
    spectra = pd.DataFrame(rows).sort_values(['jd_offset', 'filename'])
    cont = pd.read_csv(DATA / 'c5100.dat', sep=r'\s+', names=['jd_offset', 'continuum', 'continuum_error'])
    if not np.isfinite(cont.to_numpy()).all() or (cont.continuum_error < 0).any():
        raise ValueError('Invalid continuum')
    lower = max(spectra.jd_offset.min(), cont.jd_offset.min())
    upper = min(spectra.jd_offset.max(), cont.jd_offset.max())
    cont = cont[cont.jd_offset.between(lower, upper)]
    spectra = spectra[spectra.jd_offset.between(lower, upper)]
    return spectra, combine_dates(spectra, [*WINDOWS, 'profile_total']), combine_dates(cont, ['continuum']), paths

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    spectra, lines, cont, paths = load_native()
    spectra.to_csv(OUT / 'spectra.csv', index=False)
    lines.to_csv(OUT / 'line_native.csv', index=False)
    cont.to_csv(OUT / 'continuum_native.csv', index=False)
    legacy = pd.read_csv(ROOT / 'agn_surd_project/processed/ngc5548_hb_velocity_bins.csv')
    oldcols = ['blue_wing_flux', 'core_flux', 'red_wing_flux']
    old = legacy.sort_values(['mjd'] + oldcols)[oldcols].to_numpy()
    new = spectra.sort_values(['jd_offset', 'blue', 'core', 'red'])[['blue', 'core', 'red']].to_numpy()
    np.testing.assert_allclose(new, 2 * old, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(lines.profile_total, lines[['blue', 'core', 'red']].sum(axis=1), atol=1e-12)
    summary = []
    for name, frame in [('continuum', cont), ('line', lines)]:
        summary.append({'series': name, 'native_dates': len(frame), 'source_observations': int(frame.n_observations.sum()), 'first_jd_offset': frame.jd_offset.min(), 'last_jd_offset': frame.jd_offset.max(), 'largest_gap_days': np.diff(frame.jd_offset).max()})
    pd.DataFrame(summary).to_csv(OUT / 'data_table.csv', index=False)
    table_path = ROOT / 'audit/iw_bp_96_tab.txt'
    table = table_path.read_text().split('Table 2.')[0]
    keys = {(int(m[1]), m[2]) for m in re.finditer(r'(?m)^(\d{4})\s+([a-z]+)\s+\d', table)}
    unmatched = []
    for path in paths:
        m = re.fullmatch(r'n5(\d{4})([a-z]+)\.spc', path.name)
        if (int(m[1]), m[2]) not in keys:
            unmatched.append(path.name)
    manifest = {
        'status': 'Candidate pipeline; spectral date convention requires archive observation log confirmation',
        'time_coordinate': 'JD less 2400000; not MJD',
        'spectral_date_rule': 'Filename numeric field less 10000, inherited from existing reduction; no fractional times available',
        'time_scale': 'Archive labels JD; UTC, TT and barycentric correction unspecified',
        'velocity_km_s': WINDOWS, 'redshift_adopted': 0.017175, 'reference_wavelength_angstrom': 4861.33,
        'integration': 'Sum flux density times 2 Angstrom for pixels whose centres fall in each interval; lower inclusive and upper exclusive',
        'line_units': 'Arbitrary profile flux density units times Angstrom',
        'continuum_units': '1e-15 erg / s / cm^2 / Angstrom',
        'duplicate_policy': 'Arithmetic mean at identical dates; propagated errors assume independent exposures',
        'uncertainty_limit': 'Pixel covariance and shared calibration covariance unavailable; propagated errors omit them',
        'total_definition': 'Sum of adopted components; distinct from archive integrated Hbeta with narrow line contribution',
        'legacy_comparison_max_error': float(np.max(np.abs(new - 2 * old))),
        'corrected_table_date_matches': len(paths) - len(unmatched),
        'spectra_absent_from_corrected_table': unmatched,
        'sources': ['https://www.asc.ohio-state.edu/astronomy/agnwatch/papers/iw_bp_96_tab.ps.gz', 'https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/spectra/', 'https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/'],
        'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths + [DATA / 'c5100.dat', table_path, ROOT / 'audit/iw_bp_96_tab.ps', Path(__file__).resolve()]},
        'versions': {'numpy': np.__version__, 'pandas': pd.__version__},
    }
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(summary).to_string(index=False))
    print('Legacy flux scaling and component sum checks passed')

if __name__ == '__main__':
    main()
