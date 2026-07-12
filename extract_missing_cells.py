import nbformat
import re

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V8.ipynb', 'r') as f:
    nb = nbformat.read(f, as_version=4)

missing_code = []
for cell in nb.cells:
    if cell.cell_type == 'code':
        source = cell.source
        # We need the cells that define synergy_p1_stats_core, etc.
        # Or even better, we can just grab all cells that call run_surd_shuffle_test_extract_metrics
        # and run_surd_block_shuffle_test_extract_metrics
        if 'synergy_' in source and 'run_surd_' in source:
            missing_code.append(source)
            
print(f"Found {len(missing_code)} cells with synergy computations.")
for i, code in enumerate(missing_code):
    print(f"--- Cell {i} ---")
    print(code[:200] + "...")

