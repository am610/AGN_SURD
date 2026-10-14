"""Figures of the 26 normalized SURD components and leakage against lag."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

from collaborator_rebuild.epoch_scan import OUT

plt.rcParams.update({'font.size': 12, 'axes.labelsize': 12, 'xtick.labelsize': 11, 'ytick.labelsize': 11})
LABELS = {'total': 'Htot', 'continuum': 'F5100', 'blue': 'Hb', 'core': 'Hc', 'red': 'Hr'}
TITLES = {'common': 'common lag', 'blue_leading': 'blue wing 5 d further back',
          'red_leading': 'red wing 5 d further back'}


def short(component):
    kind, rest = component.split('(')
    return kind + '(' + ','.join(LABELS[n] for n in rest[:-1].split(',')) + ')'


def draw(group, summary, key):
    epoch, subset, target, config = key
    fig, axes = plt.subplots(4, 1, figsize=(9, 9.6), sharex=True)
    for ax, kind in zip(axes[:3], ['U', 'R', 'S']):
        for index, (component, curve) in enumerate(group[group.kind == kind].groupby('component')):
            color = plt.get_cmap('tab20')(index)
            for bins, style, width in [(2, '-', 1.6), (3, ':', 1.2)]:
                selected = curve[curve.bins == bins]
                ax.plot(selected.lag_days, selected.fraction, linestyle=style, linewidth=width, color=color,
                        label=short(component) if bins == 2 else None)
        ax.set_ylabel(kind + ' / joint MI')
        ax.set_ylim(-.02, 1.02)
        ax.legend(fontsize=9, ncol=4, loc='upper right', handlelength=1.4, borderpad=.3, labelspacing=.25, framealpha=.85)
    for bins, style in [(2, '-'), (3, ':')]:
        selected = summary[summary.bins == bins]
        axes[3].plot(selected.lag_days, selected.normalized_leakage, linestyle=style, color='black',
                     label=f'{bins} bins')
    axes[3].set_ylim(-.02, 1.02)
    axes[3].set_ylabel('Normalized leakage')
    axes[3].set_xlabel('Base lag (days)')
    axes[3].legend()
    tuples = summary.drop_duplicates('lag_days').tuple_count
    fig.suptitle(f'{epoch}  {subset}  target {LABELS[target]}  {TITLES[config]}\n'
                 f'solid: 2 bins, dotted: 3 bins; daily grid tuples {tuples.min()} to {tuples.max()}')
    fig.tight_layout()
    return fig


def plot_all():
    summary = pd.read_csv(OUT / 'support_and_leakage.csv')
    components = pd.read_csv(OUT / 'individual_components.csv.gz')
    keys = ['epoch', 'subset', 'target', 'config']
    for config in summary.config.unique():
        (OUT / 'figures' / config).mkdir(parents=True, exist_ok=True)
        with PdfPages(OUT / 'figures' / f'atlas_{config}.pdf') as pdf:
            for key, group in components[components.config == config].groupby(keys, sort=True):
                rows = summary[(summary.epoch == key[0]) & (summary.subset == key[1]) &
                               (summary.target == key[2]) & (summary.config == key[3])]
                fig = draw(group, rows, key)
                fig.savefig(OUT / 'figures' / config / ('_'.join(map(str, key[:3])) + '.png'), dpi=100)
                pdf.savefig(fig)
                plt.close(fig)


if __name__ == '__main__':
    plot_all()
