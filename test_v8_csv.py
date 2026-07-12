import nbformat

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V8.ipynb', 'r') as f:
    nb = nbformat.read(f, as_version=4)

for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'df_robustness_core_combined = pd.DataFrame(' in cell.source:
            print("--- Found CSV generation cell ---")
            lines = cell.source.split('\n')
            for line in lines:
                if 'synergy' in line and ':' in line:
                    print(line)
