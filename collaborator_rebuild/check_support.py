"""Save provenance and sampling diagnostics without changing scientific inputs."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
import numpy as np
import pandas as pd
from pipeline.adopted_data import load_adopted
from collaborator_rebuild.lag_design import campaign_labels, expanded_predictors

BASE = Path(__file__).resolve().parent
URL = 'https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/iw-hb.lcv'


def main():
    out = BASE / 'support_review'
    out.mkdir(exist_ok=True)
    raw = urlopen(URL, timeout=30).read()
    source = BASE / 'sources' / 'historical_revised_hbeta.dat'
    source.write_bytes(raw)
    revised = pd.read_csv(source, sep=r'\s+', names=['archive_date', 'value', 'error'])
    if not np.isfinite(revised.to_numpy()).all() or (revised.error < 0).any():
        raise ValueError('Invalid archive series')
    revised['jd_offset'] = revised.archive_date + 40000
    averaged = revised.groupby('jd_offset').agg(value=('value', 'mean'),
                                               observations=('value', 'size'))
    line, continuum = load_adopted()
    comparison = line[['jd_offset', 'profile_total']].merge(averaged, on='jd_offset',
                                                         validate='one_to_one')
    if len(comparison) != len(line):
        raise ValueError('Missing adopted dates in revised line series')
    comparison['flux_ratio'] = comparison.value / comparison.profile_total
    comparison.to_csv(out / 'total_comparison.csv', index=False)
    sources = {'continuum': continuum.rename(columns={'continuum': 'value'})}
    for name, column in [('total', 'profile_total'), ('blue', 'blue'), ('core', 'core'), ('red', 'red')]:
        sources[name] = line.rename(columns={column: 'value'})
    sources = {n: f[['jd_offset', 'value']].assign(epoch=campaign_labels(f.jd_offset))
               for n, f in sources.items()}
    rows = []
    # Five days is a declared exploratory spacing, not an inferred physical delay.
    for year in range(1989, 1994):
        series = {n: f[f.epoch == year] for n, f in sources.items()}
        for subset, names in [('total_wings', ['total', 'blue', 'core', 'red']),
                              ('continuum_wings', ['continuum', 'blue', 'core', 'red'])]:
            for target in names:
                for lag in range(31):
                    for expanded in [False, True]:
                        x, labels = expanded_predictors(series, names, series[target].jd_offset,
                                                        lag, [5] if expanded else [])
                        valid = np.isfinite(x).all(axis=1)
                        for bins in [2, 3]:
                            edges = {n: np.r_[-np.inf, np.unique(np.quantile(series[n].value,
                                                    np.arange(1, bins) / bins)), np.inf] for n in names}
                            data = np.column_stack([series[target].value.to_numpy(), x])[valid]
                            hist, _ = np.histogramdd(data, bins=[edges[target]] +
                                                    [edges[label.split('@')[0]] for label in labels])
                            occupied = np.count_nonzero(hist)
                            rows.append(dict(epoch=year, subset=subset, target=target, lag_days=lag,
                                             expanded=expanded, extra_delay_days=5 if expanded else 0,
                                             bins=bins, predictors=len(labels), tuples=int(valid.sum()),
                                             cells=hist.size, occupied_cells=int(occupied),
                                             singleton_cells=int(np.count_nonzero(hist == 1))))
    pd.DataFrame(rows).to_csv(out / 'sampling.csv', index=False)
    metadata = dict(source_url=URL, source_sha256=hashlib.sha256(raw).hexdigest(),
                    archive_rows=len(revised), matched_adopted_dates=len(comparison),
                    correlation=float(comparison.value.corr(comparison.profile_total)),
                    ratio_min=float(comparison.flux_ratio.min()),
                    ratio_max=float(comparison.flux_ratio.max()),
                    selection='Candidate only; no scientific total replacement',
                    date_policy='Add 40000, validated against all adopted dates; archive header differs',
                    expanded_policy='Equation 3.18 with one extra five day copy per variable; support only')
    (out / 'provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
