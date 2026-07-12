import nbformat

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', 'r') as f:
    nb = nbformat.read(f, as_version=4)

for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'def run_surd_shuffle_test_extract_metrics(' in cell.source:
            print("--- Function definition ---")
            print(cell.source[-500:])
