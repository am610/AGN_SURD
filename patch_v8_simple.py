import nbformat
import re

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V8.ipynb', 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

# Full Range
nb_full = nbformat.from_dict(nb.copy())
for cell in nb_full.cells:
    if cell.cell_type == 'code':
        source = cell.source
        if 'from google.colab import drive' in source or '!unzip -q' in source:
            cell.source = ""
            continue
        if '!git clone https://github.com/Computational-Turbulence-Group/SURD.git' in source:
            source = re.sub(r'!git clone .*', '# SURD repo should be available locally', source)
        source = source.replace('/content/agn_surd_project', '/Users/ayan/Programs/SURD/agn_surd_project')
        source = source.replace('/content/drive/MyDrive/SURDS/agn_surd_project.zip', '/Users/ayan/Programs/SURD/agn_surd_project.zip')
        source = source.replace('/content/SURD/utils', '/Users/ayan/Programs/SURD/SURD/utils')
        if 'start_mjd_strict = ' in source:
            source = re.sub(r'start_mjd_strict = .*', 'start_mjd_strict = -np.inf', source)
            source = re.sub(r'end_mjd_strict = .*', 'end_mjd_strict = np.inf', source)
        cell.source = source

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', 'w', encoding='utf-8') as f:
    nbformat.write(nb_full, f)

# Restricted Range
nb_restr = nbformat.from_dict(nb.copy())
for cell in nb_restr.cells:
    if cell.cell_type == 'code':
        source = cell.source
        if 'from google.colab import drive' in source or '!unzip -q' in source:
            cell.source = ""
            continue
        if '!git clone https://github.com/Computational-Turbulence-Group/SURD.git' in source:
            source = re.sub(r'!git clone .*', '# SURD repo should be available locally', source)
        source = source.replace('/content/agn_surd_project', '/Users/ayan/Programs/SURD/agn_surd_project')
        source = source.replace('/content/drive/MyDrive/SURDS/agn_surd_project.zip', '/Users/ayan/Programs/SURD/agn_surd_project.zip')
        source = source.replace('/content/SURD/utils', '/Users/ayan/Programs/SURD/SURD/utils')
        if 'start_mjd_strict = ' in source:
            source = re.sub(r'start_mjd_strict = .*', 'start_mjd_strict = 47512.0', source)
            source = re.sub(r'end_mjd_strict = .*', 'end_mjd_strict = 49255.0', source)
        cell.source = source

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb', 'w', encoding='utf-8') as f:
    nbformat.write(nb_restr, f)

print("Saved both notebooks based on V8.")
