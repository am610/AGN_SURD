import nbformat
import re

for filename in ['/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb']:
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    for cell in nb.cells:
        if cell.cell_type == 'code':
            source = cell.source
            source = source.replace('metrics_real_core', 'metrics_extended_core')
            source = source.replace('metrics_real_red', 'metrics_extended_red')
            source = source.replace('metrics_real_blue', 'metrics_extended_blue')
            cell.source = source

    with open(filename, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
print("Patched variable names!")
