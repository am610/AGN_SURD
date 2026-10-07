"""Prepare published corrected line fluxes separately from exploratory SURD."""
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from pipeline.adopted_data import TABLE, corrected_table, load_adopted
from pipeline.reconcile_data import combine_dates
from collaborator_rebuild.lag_design import campaign_labels

BASE = Path(__file__).resolve().parent
OUT = BASE / 'historical_prepared'
SOURCE = 'https://www.asc.ohio-state.edu/astronomy/agnwatch/papers/iw_bp_96_tab.ps.gz'


def published_lines(text=None):
    text = TABLE.read_text() if text is None else text
    rows = []
    for line in text.split('Table 2.')[0].splitlines():
        match = re.match(r'^(\d{4})\s+([a-z]+)\s+(.*)', line)
        if not match:
            continue
        values = [float(v) for v in re.findall(r'[-+]?\d+\.\d+', match[3])]
        if len(values) != 13:
            raise ValueError('Unexpected corrected table columns')
        row = dict(jd_offset=int(match[1]) + 40000, instrument=match[2],
                   table_date=int(match[1]), filename=f'n5{match[1]}{match[2]}.spc')
        for name, index in [('table_continuum', 0), ('total', 2), ('blue', 4), ('core', 7), ('red', 10)]:
            row[name], row[name + '_error'] = values[index:index+2]
        rows.append(row)
    frame = pd.DataFrame(rows)
    if len(frame) != 242 or frame.filename.duplicated().any():
        raise ValueError('Expected 242 distinct corrected observations')
    if not np.isfinite(frame.select_dtypes('number').to_numpy()).all():
        raise ValueError('Nonfinite published values')
    if (frame.filter(like='_error') <= 0).any().any():
        raise ValueError('Nonpositive published uncertainties')
    residual = frame.total - frame[['blue', 'core', 'red']].sum(axis=1)
    if abs(residual).max() > .031:
        raise ValueError('Published total differs from component sum beyond printed rounding')
    expected = corrected_table()
    np.testing.assert_array_equal(frame.filename, expected.filename)
    np.testing.assert_allclose(frame.table_continuum, expected.table_continuum, atol=1e-12)
    return frame


def prepare():
    observations = published_lines()
    _, continuum = load_adopted()
    lines = combine_dates(observations, ['total', 'blue', 'core', 'red'])
    if len(lines) != 224:
        raise ValueError('Expected 224 distinct line dates')
    sources = {'continuum': continuum.rename(columns={'continuum': 'value', 'continuum_error': 'error'})}
    for name in ['total', 'blue', 'core', 'red']:
        sources[name] = lines.rename(columns={name: 'value', name + '_error': 'error'})
    sources = {name: frame[['jd_offset', 'value', 'error', 'n_observations']].assign(
        epoch=campaign_labels(frame.jd_offset)) for name, frame in sources.items()}
    OUT.mkdir(exist_ok=True)
    observations.to_csv(OUT / 'published_observations.csv', index=False)
    inventory = []
    for name, frame in sources.items():
        for year, selected in frame.groupby('epoch'):
            directory = OUT / str(year)
            directory.mkdir(exist_ok=True)
            selected.to_csv(directory / f'{name}.csv', index=False)
            inventory.append(dict(epoch=int(year), series=name, native_dates=len(selected),
                                  source_observations=int(selected.n_observations.sum())))
    pd.DataFrame(inventory).to_csv(OUT / 'inventory.csv', index=False)
    manifest = dict(line_source=SOURCE, line_table_sha256=hashlib.sha256(TABLE.read_bytes()).hexdigest(),
                    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    definitions='Published corrected Table 1 values, not fitted profile reconstruction',
                    line_units='1e-14 erg / s / cm^2',
                    continuum_units='1e-15 erg / s / cm^2 / Angstrom',
                    continuum_choice='Original archive F5100, host retained; revised continuum is a separate sensitivity',
                    intervals_observed_angstrom={'blue': [4870, 4920], 'core': [4920, 4960], 'red': [4960, 5010]},
                    interval_authority='Wanders and Peterson 1996 section 3.2, page 180; figure caption differs',
                    time_coordinate='JD less 2400000; no fractional line exposure times established',
                    duplicate_policy='Equal date means, uncertainties in quadrature; shared covariance unavailable',
                    line_observations=len(observations), line_dates=len(lines),
                    total_policy='Use published total column, retain rounding rather than force exact sum',
                    interpretation='Total and components share spectra and nearly sum; not independent measurements')
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(inventory).to_string(index=False))
    return sources


if __name__ == '__main__':
    prepare()
