import nbformat

for filename in ['/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb']:
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    for cell in nb.cells:
        if cell.cell_type == 'code':
            cell.source = cell.source.replace('/content/slide_figures', '/Users/ayan/Programs/SURD/agn_surd_project/slide_figures')
            cell.source = cell.source.replace('/content/', '/Users/ayan/Programs/SURD/agn_surd_project/')

    with open(filename, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
print("Patched /content/ references!")
