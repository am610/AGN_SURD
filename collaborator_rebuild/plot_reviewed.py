"""Render existing descriptive atoms with explicit support and missing cases."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parent


def main():
    source = BASE / 'historical_results'
    out = BASE / 'reviewed_figures'
    out.mkdir(exist_ok=True)
    support = pd.read_csv(source / 'support_and_leakage.csv')
    atoms = pd.read_csv(source / 'individual_components.csv')
    keys = ['epoch', 'subset', 'target', 'scenario']
    inventory = []
    for key, scans in support.groupby(keys):
        mask = np.ones(len(atoms), dtype=bool)
        for column, value in zip(keys, key):
            mask &= atoms[column].eq(value).to_numpy()
        group = atoms[mask]
        fig, axs = plt.subplots(5, 1, figsize=(12, 15), sharex=True)
        for ax, kind in zip(axs[:3], ['U', 'R', 'S']):
            for index, (label, curve) in enumerate(group[group.kind == kind].groupby('component')):
                full = curve.set_index('lag_days').fraction.reindex(range(31))
                ax.plot(full.index, full, marker='.', markersize=4,
                        color=plt.cm.tab20(index / 20),
                        linestyle=['solid', 'dashed', 'dotted'][index % 3], label=label)
            ax.set_ylabel(kind + ' / joint MI')
            ax.set_ylim(-.02, 1.02)
            if ax.lines:
                ax.legend(fontsize=7, ncol=3, loc='upper right')
        scans = scans.set_index('lag_days').reindex(range(31))
        good = scans.status.eq('ok')
        axs[3].plot(scans.index, scans.normalized_leakage.where(good), '.-', color='black')
        axs[3].set_ylabel('Leakage / target entropy')
        axs[3].set_ylim(-.02, 1.02)
        axs[4].plot(scans.index, scans.tuple_count, '.-', color='black')
        axs[4].axhline(32, color='red', linestyle='dotted', label='Numerical guard only')
        axs[4].legend(fontsize=8)
        axs[4].set_ylabel('Native target tuples')
        axs[4].set_xlabel('Base lag in observed days; zero lag is an algebraic control')
        for ax in axs:
            ax.set_xlim(-.5, 30.5)
            for lag in scans.index[~good]:
                ax.axvspan(lag - .5, lag + .5, color='grey', alpha=.12)
        if not good.any():
            axs[1].text(.5, .5, 'No lag passes the numerical support guard',
                        transform=axs[1].transAxes, ha='center')
        fig.suptitle(' '.join(map(str, key)) + '\nDescriptive estimates; grey indicates unsupported lags')
        fig.tight_layout()
        filename = '_'.join(map(str, key)) + '.png'
        fig.savefig(out / filename, dpi=110)
        plt.close(fig)
        inventory.append(dict(zip(keys, key), file=filename, usable_lags=int(good.sum())))
    pd.DataFrame(inventory).to_csv(out / 'inventory.csv', index=False)
    print(f'Rendered {len(inventory)} explicit cases')


if __name__ == '__main__':
    main()
