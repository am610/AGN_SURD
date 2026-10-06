import nbformat
import re

for filename in ['/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb']:
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    for cell in nb.cells:
        if cell.cell_type == 'code':
            if 'def run_surd_shuffle_test_extract_metrics(' in cell.source:
                cell.source = cell.source.replace("'p5_synergy': np.percentile(all_shuffled_s12, 5, axis=0),\n        'p95_synergy': np.percentile(all_shuffled_s12, 95, axis=0)",
                                                  "'p5_synergy': np.percentile(all_shuffled_s12, 5, axis=0),\n        'p95_synergy': np.percentile(all_shuffled_s12, 95, axis=0),\n        'p2_5_synergy': np.percentile(all_shuffled_s12, 2.5, axis=0),\n        'p16_synergy': np.percentile(all_shuffled_s12, 16, axis=0),\n        'p84_synergy': np.percentile(all_shuffled_s12, 84, axis=0),\n        'p97_5_synergy': np.percentile(all_shuffled_s12, 97.5, axis=0)")
                cell.source = cell.source.replace("'p5_leak': np.percentile(all_shuffled_leak, 5, axis=0),\n        'p95_leak': np.percentile(all_shuffled_leak, 95, axis=0)",
                                                  "'p5_leak': np.percentile(all_shuffled_leak, 5, axis=0),\n        'p95_leak': np.percentile(all_shuffled_leak, 95, axis=0),\n        'p2_5_leak': np.percentile(all_shuffled_leak, 2.5, axis=0),\n        'p16_leak': np.percentile(all_shuffled_leak, 16, axis=0),\n        'p84_leak': np.percentile(all_shuffled_leak, 84, axis=0),\n        'p97_5_leak': np.percentile(all_shuffled_leak, 97.5, axis=0)")

    with open(filename, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
print("Patched percentiles in shuffle function!")
