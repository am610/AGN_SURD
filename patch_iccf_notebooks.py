import nbformat

def patch_iccf_in_notebook(path):
    print(f"Patching ICCF calls in: {path}")
    with open(path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)
        
    patched = False
    for cell in nb.cells:
        if cell.cell_type == 'code' and 'compute_and_plot_iccf' in cell.source:
            code = cell.source
            # Simple string replacements
            code = code.replace(
                "cont_zscore_fiducial, blue_zscore_fiducial,",
                "blue_zscore_fiducial, cont_zscore_fiducial,"
            )
            code = code.replace(
                "cont_zscore_fiducial, core_zscore_fiducial,",
                "core_zscore_fiducial, cont_zscore_fiducial,"
            )
            code = code.replace(
                "cont_zscore_fiducial, red_zscore_fiducial,",
                "red_zscore_fiducial, cont_zscore_fiducial,"
            )
            cell.source = code
            patched = True
            print("Successfully patched code cell.")
            
    if patched:
        with open(path, 'w', encoding='utf-8') as f:
            nbformat.write(nb, f)
    else:
        print("Could not find cell to patch.")

patch_iccf_in_notebook('/Users/ayan/Programs/SURD/3_SURDS_AGN_V11_Strict_Overlap.ipynb')
patch_iccf_in_notebook('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb')
