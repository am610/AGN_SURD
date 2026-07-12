import nbformat
import sys

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

defined = set()
for idx, cell in enumerate(nb.cells):
    if cell.cell_type == 'code':
        source = cell.source
        
        # very simple check
        if 'metrics_extended_red = lag_scan_target3' in source:
            if 'X_real_red_target' not in defined:
                print(f"Cell {idx} uses X_real_red_target but it's not defined yet!")
        
        if 'X_real_red_target = np.vstack' in source:
            defined.add('X_real_red_target')
            
        if 'X_real_blue_target = np.vstack' in source:
            defined.add('X_real_blue_target')
            
        if 'metrics_extended_blue = lag_scan_target3' in source:
            if 'X_real_blue_target' not in defined:
                print(f"Cell {idx} uses X_real_blue_target but it's not defined yet!")
