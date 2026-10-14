"""Tables and descriptive extrema for the per epoch SURD scans."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from collaborator_rebuild.epoch_scan import OUT, SUBSETS
from collaborator_rebuild.lag_design import CAMPAIGNS

SECTION = Path(__file__).resolve().parent / 'results_section'
NAMES = {'total': 'H', 'continuum': 'F', 'blue': 'B', 'core': 'C', 'red': 'R'}
TARGETS = {'total': 'Total', 'continuum': 'F5100', 'blue': 'Blue', 'core': 'Core', 'red': 'Red'}
KEYS = ['epoch', 'subset', 'target', 'config']


def write_tabular(path, columns, header, rows):
    """Write a complete booktabs tabular; a row may be the literal string \\midrule."""
    lines = ['\\begin{tabular}{' + columns + '}', '\\toprule', header, '\\midrule'] + [r for r in rows if r]
    lines += ['\\bottomrule', '\\end{tabular}']
    Path(path).write_text('\n'.join(lines) + '\n')


def tex_component(component):
    kind, rest = component.split('(')
    return '$' + kind + '(' + ','.join(NAMES[n] for n in rest[:-1].split(',')) + ')$'


def extrema(components, summary):
    rows = []
    for key, group in components.groupby(KEYS + ['bins']):
        leak = summary[(summary.epoch == key[0]) & (summary.subset == key[1]) & (summary.target == key[2]) &
                       (summary.config == key[3]) & (summary.bins == key[4])].sort_values('lag_days')
        low, high = leak.loc[leak.normalized_leakage.idxmin()], leak.loc[leak.normalized_leakage.idxmax()]
        row = dict(zip(KEYS + ['bins'], key), leakage_min=low.normalized_leakage, leakage_min_lag=int(low.lag_days),
                   leakage_max=high.normalized_leakage, leakage_max_lag=int(high.lag_days),
                   leakage_mean=leak.normalized_leakage.mean(), tuples_min=int(leak.tuple_count.min()),
                   tuples_max=int(leak.tuple_count.max()))
        for kind in 'URS':
            part = group[group.kind == kind]
            best = part.loc[part.fraction.idxmax()]
            row.update({f'top_{kind}': best.component, f'top_{kind}_lag': int(best.lag_days),
                        f'top_{kind}_fraction': best.fraction})
            at_low = part[part.lag_days == low.lag_days]
            row[f'{kind}_at_leakage_min'] = at_low.fraction.sum()
        at_low = group[group.lag_days == low.lag_days].sort_values('fraction')
        row.update(largest_at_leakage_min=at_low.iloc[-1].component, largest_at_leakage_min_fraction=at_low.iloc[-1].fraction)
        at_high = group[group.lag_days == high.lag_days].sort_values('fraction')
        row.update(largest_at_leakage_max=at_high.iloc[-1].component, largest_at_leakage_max_fraction=at_high.iloc[-1].fraction)
        rows.append(row)
    return pd.DataFrame(rows)


def bin_sensitivity(components):
    wide = components.pivot_table(index=KEYS + ['lag_days', 'component'], columns='bins', values='fraction')
    wide['absolute'] = (wide[2] - wide[3]).abs()
    tv = wide.groupby(KEYS + ['lag_days']).absolute.sum() / 2
    return tv.groupby(KEYS).median().rename('median_total_variation').reset_index()


LU_TABLE6 = {1989: ('Year 1', 19.7, 1.5, 1.5), 1990: ('Year 2', 18.6, 2.1, 2.3), 1991: ('Year 3', 15.9, 2.9, 2.5),
             1992: ('Year 4', 11.0, 1.9, 2.0), 1993: ('Year 5', 13.0, 1.6, 1.4)}


def continuum_peaks(components):
    """Lag of the largest unique continuum fraction for each line target."""
    part = components[(components.subset == 'continuum_wings') & (components.component == 'U(continuum)') &
                      (components.config == 'common') & components.target.isin(['blue', 'core', 'red'])]
    rows = []
    for key, group in part.groupby(['epoch', 'target', 'bins']):
        best = group.loc[group.fraction.idxmax()]
        rows.append(dict(epoch=key[0], target=key[1], bins=key[2], peak_lag=int(best.lag_days),
                         peak_fraction=best.fraction, detected=bool(best.fraction >= 0.05)))
    return pd.DataFrame(rows)


CONFIG_TEXT = {'common': 'common predictor lags',
               'blue_leading': 'uncommon predictor lags with the blue wing five days further back',
               'red_leading': 'uncommon predictor lags with the red wing five days further back'}
SET_TEXT = {'total_wings': 'total, blue, core and red H$\\beta$', 'continuum_wings': 'F5100, blue, core and red H$\\beta$'}


SELECTED = [('common', 1989, 'continuum_wings', 'core'), ('common', 1992, 'continuum_wings', 'red'),
            ('common', 1993, 'total_wings', 'core'), ('blue_leading', 1992, 'continuum_wings', 'core')]


def figure_list(summary):
    """Example figures from the atlas pages; page order matches plot_epoch.plot_all."""
    lines = []
    for config, want_epoch, want_subset, want_target in SELECTED:
        keys = sorted(summary[summary.config == config][['epoch', 'subset', 'target']].drop_duplicates().itertuples(index=False))
        for page, (epoch, subset, target) in enumerate(keys, start=1):
            if (epoch, subset, target) != (want_epoch, want_subset, want_target):
                continue
            rows = summary[(summary.epoch == epoch) & (summary.subset == subset) & (summary.target == target) & (summary.config == config)]
            tuples = rows.drop_duplicates('lag_days').tuple_count
            lines += ['\\begin{SURDfigure}', '\\centering',
                      f'\\includegraphics[page={page},width=\\textwidth]{{\\SURDEpochPath/figures/atlas_{config}.pdf}}',
                      f'\\caption{{Campaign {epoch}, variable set {SET_TEXT[subset]}, future target {TARGETS[target]}, {CONFIG_TEXT[config]}. '
                      f'The 26 normalized components (U, R, S panels) and normalized leakage against base lag. '
                      f'Solid curves use two bins and dotted curves three bins; the daily grid supplies {tuples.min()} to {tuples.max()} tuples per lag. '
                      f'B, C, R and H or F denote blue, core, red and total or continuum.}}',
                      f'\\label{{fig:ep_{config}_{epoch}_{subset}_{target}}}', '\\end{SURDfigure}']
    return '\n'.join(lines) + '\n'


def main():
    summary = pd.read_csv(OUT / 'support_and_leakage.csv')
    components = pd.read_csv(OUT / 'individual_components.csv.gz')
    inventory = pd.read_csv(OUT / 'native_inventory.csv')
    ext = extrema(components, summary)
    ext.to_csv(OUT / 'descriptive_extrema.csv', index=False)
    tv = bin_sensitivity(components)
    tv.to_csv(OUT / 'bin_sensitivity.csv', index=False)
    leak = summary.groupby(KEYS + ['bins']).normalized_leakage.mean().rename('mean_leakage').reset_index()
    leak.to_csv(OUT / 'mean_leakage.csv', index=False)
    SECTION.mkdir(exist_ok=True)
    # Support table
    rows = []
    for year, (lo, hi) in CAMPAIGNS.items():
        inv = inventory[inventory.epoch == year].set_index('series').native_dates
        sub = summary[(summary.epoch == year) & (summary.config == 'common') & (summary.bins == 2)]
        rows.append(f'{year} & {lo}--{hi} & {inv.continuum} & {inv.total} & '
                    f'{sub.tuple_count.min()}--{sub.tuple_count.max()} \\\\')
    write_tabular(SECTION / 'epoch_support_rows.tex', 'rlrrr',
                  'Epoch & JD$-2400000$ & $N_F$ & $N_{\\rm line}$ & Grid tuples \\\\', rows)
    # Extrema tables, one per subset
    common = ext[(ext.config == 'common') & (ext.bins == 2)]
    for subset in SUBSETS:
        rows = []
        for epoch in CAMPAIGNS:
            for target in SUBSETS[subset]:
                r = common[(common.epoch == epoch) & (common.subset == subset) & (common.target == target)].iloc[0]
                cells = [f'${r.leakage_min:.2f}$ ({r.leakage_min_lag})', f'${r.leakage_max:.2f}$ ({r.leakage_max_lag})']
                for kind in 'URS':
                    cells.append(f'{tex_component(r[f"top_{kind}"])} ${r[f"top_{kind}_fraction"]:.2f}$ ({r[f"top_{kind}_lag"]})')
                rows.append(f'{epoch} & {TARGETS[target]} & ' + ' & '.join(cells) + ' \\\\')
            rows.append('\\midrule' if epoch != 1993 else '')
        write_tabular(SECTION / f'epoch_extrema_{subset}.tex', 'rlcccccc',
                      'Epoch & Target & $\\ell_{\\min}$ (d) & $\\ell_{\\max}$ (d) & Largest $U$ & Largest $R$ & Largest $S$ \\\\', rows)
    # Config comparison
    two = leak[leak.bins == 2].pivot_table(index=['epoch', 'subset', 'target'], columns='config', values='mean_leakage')
    two['blue_minus_red'] = two.blue_leading - two.red_leading
    two.to_csv(OUT / 'config_leakage_comparison.csv')
    rows = []
    for epoch in CAMPAIGNS:
        for subset in SUBSETS:
            part = two.loc[(epoch, subset)]
            lines = part.drop(index=[t for t in ('continuum', 'total') if t in part.index])
            rows.append(f'{epoch} & {subset.replace("_wings", "").replace("continuum", "F5100").replace("total", "Total")} & '
                        f'${part.common.mean():.3f}$ & ${part.blue_leading.mean():.3f}$ & ${part.red_leading.mean():.3f}$ & '
                        f'${lines.blue_minus_red.mean():+.3f}$ & {int((lines.blue_minus_red < 0).sum())} of 3 \\\\')
    write_tabular(SECTION / 'epoch_config_rows.tex', 'rlccccc',
                  'Epoch & Set & Common & Blue lead & Red lead & $\\Delta\\ell$ & Targets \\\\', rows)
    # Bin sensitivity
    rows = []
    for epoch in CAMPAIGNS:
        cells = []
        for subset in SUBSETS:
            part = tv[(tv.epoch == epoch) & (tv.subset == subset)]
            cells += [f'${part[part.config == c].median_total_variation.median():.2f}$' for c in ('common', 'blue_leading', 'red_leading')]
        rows.append(f'{epoch} & ' + ' & '.join(cells) + ' \\\\')
    write_tabular(SECTION / 'epoch_bin_rows.tex', 'rcccccc',
                  ' & \\multicolumn{3}{c}{F5100 set} & \\multicolumn{3}{c}{Total set} \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\nEpoch & Common & Blue lead & Red lead & Common & Blue lead & Red lead \\\\', rows)
    (SECTION / 'epoch_figures.tex').write_text(figure_list(summary))
    peaks = continuum_peaks(components)
    peaks.to_csv(OUT / 'continuum_peak_lags.csv', index=False)
    rows = []
    for epoch in CAMPAIGNS:
        label, lag, up, down = LU_TABLE6[epoch]
        cells = []
        for bins in (2, 3):
            for target in ('blue', 'core', 'red'):
                r = peaks[(peaks.epoch == epoch) & (peaks.bins == bins) & (peaks.target == target)].iloc[0]
                cells.append(f'{r.peak_lag} ({r.peak_fraction:.2f})' if r.detected else '--')
        rows.append(f'{epoch} & ${lag:.1f}^{{+{up}}}_{{-{down}}}$ & ' + ' & '.join(cells) + ' \\\\')
    write_tabular(SECTION / 'epoch_lu_rows.tex', 'rcccccccc',
                  ' & & \\multicolumn{3}{c}{Two bins} & \\multicolumn{3}{c}{Three bins} \\\\\n\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\nEpoch & $\\tau_{\\rm H\\beta}$ (Lu) & Blue & Core & Red & Blue & Core & Red \\\\', rows)
    print(common.groupby(['subset', 'target'])[['leakage_min_lag', 'leakage_max_lag']].median())
    print(two.round(3).to_string())
    print(tv.groupby(['subset', 'config']).median_total_variation.median())


if __name__ == '__main__':
    main()
