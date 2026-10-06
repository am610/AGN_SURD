import nbformat

for filename in ['/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb']:
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Filter out the cell that creates the google_drive_figures_dir
    filtered_cells = []
    for cell in nb.cells:
        if cell.cell_type == 'code' and 'google_drive_figures_dir =' in cell.source and 'shutil.copy' in cell.source:
            continue
        filtered_cells.append(cell)
        
    nb.cells = filtered_cells

    with open(filename, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
print("Removed Google Drive copy cells from both notebooks.")
