"""Combine one factor detection power sensitivity runs into one audit table."""
import json
from pathlib import Path

import pandas as pd

from .adopted_data import ROOT


def summarize(root=None):
    root = Path(root) if root else (
        ROOT / 'agn_surd_project/processed/reconciled/detection_power_sensitivity'
    )
    frames = []
    for folder in sorted(path for path in root.iterdir() if path.is_dir()):
        manifest = json.loads((folder / 'manifest.json').read_text())
        frame = pd.read_csv(folder / 'detection_power_grid_summary.csv')
        frame.insert(0, 'sensitivity_case', folder.name)
        frame.insert(1, 'response_width_days', manifest['response_width_days'])
        frame.insert(2, 'process_timescale_days', manifest['process_timescale_days'])
        frame.insert(3, 'noise_scale', manifest['noise_scale'])
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    result.to_csv(root / 'sensitivity_summary.csv', index=False)
    summary = {
        'status': 'One factor sensitivity grid completed',
        'cases': sorted(result['sensitivity_case'].unique().tolist()),
        'interpretation': 'Each case changes one factor from the central width 5, timescale 50, noise scale 1 configuration',
        'realizations_per_cell': 25,
    }
    (root / 'manifest.json').write_text(json.dumps(summary, indent=2) + '\n')
    return result


if __name__ == '__main__':
    summarize()
