import nbformat
import sys
import re

def process_notebook(input_path, output_path, mode):
    print(f"Processing {input_path} -> {output_path} (Mode: {mode})")
    with open(input_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    cells_to_keep = []
    
    for cell in nb.cells:
        if cell.cell_type == 'code':
            source = cell.source
            
            # Remove Colab specific stuff
            if 'from google.colab import drive' in source:
                continue # Skip this cell entirely
            if '!unzip -q' in source and 'drive/MyDrive' in source:
                continue
            if '!git clone https://github.com/Computational-Turbulence-Group/SURD.git' in source:
                # We can keep the pip install pymp-pypi, but skip the git clone if local
                source = re.sub(r'!git clone .*', '# SURD repo should be available locally', source)
                
            # Fix paths
            source = source.replace('/content/agn_surd_project', '/Users/ayan/Programs/SURD/agn_surd_project')
            source = source.replace('/content/drive/MyDrive/SURDS/agn_surd_project.zip', '/Users/ayan/Programs/SURD/agn_surd_project.zip')
            source = source.replace('/content/SURD/utils', '/Users/ayan/Programs/SURD/SURD/utils')
            
            # Apply MJD logic
            if 'prepare_combined_data' in source or 'Critical Correction: Restrict data to MJD' in source:
                if mode == 'full':
                    # Remove the MJD restriction
                    source = re.sub(r'start_mjd_strict = 47512\.0', 'start_mjd_strict = -np.inf', source)
                    source = re.sub(r'end_mjd_strict = 49255\.0', 'end_mjd_strict = np.inf', source)
                elif mode == 'restricted':
                    # Ensure it has it (it should from V9)
                    source = re.sub(r'start_mjd_strict = -np.inf', 'start_mjd_strict = 47512.0', source)
                    source = re.sub(r'end_mjd_strict = np.inf', 'end_mjd_strict = 49255.0', source)
                    
                    # Fix ffill/bfill if it exists
                    source = source.replace('series_to_interp.values)', 'series_to_interp.values, left=np.nan, right=np.nan)')
                    
            cell.source = source
        cells_to_keep.append(cell)
        
    nb.cells = cells_to_keep
    
    with open(output_path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
    print(f"Saved {output_path}")

if __name__ == "__main__":
    input_file = '/Users/ayan/Programs/SURD/3_SURDS_AGN_V9.ipynb'
    process_notebook(input_file, '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', 'full')
    process_notebook(input_file, '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb', 'restricted')
