import nbformat

for filename in ['/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb']:
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    for cell in nb.cells:
        if cell.cell_type == 'code':
            if '# --- Core Hβ Target Summary ---' in cell.source:
                cell.source = cell.source.replace('lags_fiducial[', 'lags_fiducial_extended[')
                cell.source = cell.source.replace('lags_fiducial ==', 'lags_fiducial_extended ==')
                cell.source = cell.source.replace('(lags_fiducial >=', '(lags_fiducial_extended >=')
                cell.source = cell.source.replace('(lags_fiducial >', '(lags_fiducial_extended >')
                cell.source = cell.source.replace('(lags_fiducial <=', '(lags_fiducial_extended <=')

    with open(filename, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
print("Patched summary cell!")
