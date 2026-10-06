import nbformat
import re

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

# We will just insert lags_fiducial_extended = np.arange(1, 121) at the beginning of the cell that first uses it, or just right after lags_fiducial is defined.
for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'lags_fiducial =' in cell.source or 'lags_fiducial_extended' in cell.source:
            if 'lags_fiducial_extended = np.arange(1, 121)' not in cell.source:
                # Add it at the top of the cell that first uses it if it's not defined
                # Or just put it in the cell that defines lags_fiducial
                pass

        # Let's just blindly inject lags_fiducial_extended = np.arange(1, 121) before it's used
        if 'X_real_core_target = np.vstack([' in cell.source:
            if 'lags_fiducial_extended =' not in cell.source:
                cell.source = "lags_fiducial_extended = np.arange(1, 121)\n" + cell.source
                
        # Also, check if there are other out-of-order variables!
        # metrics_extended_red is computed in cell 8004
        # metrics_extended_blue is computed in cell 8013

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb', 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'X_real_core_target = np.vstack([' in cell.source:
            if 'lags_fiducial_extended =' not in cell.source:
                cell.source = "lags_fiducial_extended = np.arange(1, 121)\n" + cell.source

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb', 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)
print("Patched out of order variables")
