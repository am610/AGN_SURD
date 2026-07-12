import nbformat

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V8.ipynb', 'r') as f:
    nb = nbformat.read(f, as_version=4)

for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'def run_surd_shuffle_test_extract_metrics' in cell.source:
            print("--- Found function def ---")
            lines = cell.source.split('\n')
            for line in lines:
                if 'synergy' in line and ':' in line and '{' in line or 'p' in line:
                    pass
            print('\n'.join([line for line in lines if 'synergy' in line and ':' in line]))
