import nbformat
import re

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V8.ipynb', 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

for cell in nb.cells:
    if cell.cell_type == 'code':
        source = cell.source
        
        # Remove colab and update paths
        if 'from google.colab import drive' in source or '!unzip -q' in source:
            cell.source = ""
            continue
        if '!git clone https://github.com/Computational-Turbulence-Group/SURD.git' in source:
            source = re.sub(r'!git clone .*', '# SURD repo should be available locally', source)
            
        source = source.replace('/content/agn_surd_project', '/Users/ayan/Programs/SURD/agn_surd_project')
        source = source.replace('/content/drive/MyDrive/SURDS/agn_surd_project.zip', '/Users/ayan/Programs/SURD/agn_surd_project.zip')
        source = source.replace('/content/SURD/utils', '/Users/ayan/Programs/SURD/SURD/utils')
        
        # We are creating the Full Range version
        if 'prepare_combined_data' in source or 'Critical Correction: Restrict data to MJD' in source:
            source = re.sub(r'start_mjd_strict = 47512\.0', 'start_mjd_strict = -np.inf', source)
            source = re.sub(r'end_mjd_strict = 49255\.0', 'end_mjd_strict = np.inf', source)

        # Patch the run_surd_shuffle_test_extract_metrics to compute all percentiles
        if 'def run_surd_shuffle_test_extract_metrics(' in source or 'def run_surd_block_shuffle_test_extract_metrics(' in source:
            # Add p5, p95, p16, p84, p2_5, p97_5 for both synergy and leak
            # Let's just blindly insert them if we see 'p16_synergy' or 'p5_synergy'
            # Actually, I can just replace the return dict completely!
            # Since the function is complex, let's use regex to find the return dict
            replacement = """        'median_synergy': np.median(all_shuffled_s12, axis=0),
        'p2_5_synergy': np.percentile(all_shuffled_s12, 2.5, axis=0),
        'p5_synergy': np.percentile(all_shuffled_s12, 5, axis=0),
        'p16_synergy': np.percentile(all_shuffled_s12, 16, axis=0),
        'p84_synergy': np.percentile(all_shuffled_s12, 84, axis=0),
        'p95_synergy': np.percentile(all_shuffled_s12, 95, axis=0),
        'p97_5_synergy': np.percentile(all_shuffled_s12, 97.5, axis=0),
        'median_leak': np.median(all_shuffled_leak, axis=0),
        'p2_5_leak': np.percentile(all_shuffled_leak, 2.5, axis=0),
        'p5_leak': np.percentile(all_shuffled_leak, 5, axis=0),
        'p16_leak': np.percentile(all_shuffled_leak, 16, axis=0),
        'p84_leak': np.percentile(all_shuffled_leak, 84, axis=0),
        'p95_leak': np.percentile(all_shuffled_leak, 95, axis=0),
        'p97_5_leak': np.percentile(all_shuffled_leak, 97.5, axis=0)"""
            
            # The dictionary is usually returned directly or assigned
            source = re.sub(r"'median_synergy':.*?np\.median\(all_shuffled_s12, axis=0\).*?(?=\n    \})", replacement, source, flags=re.DOTALL)

        cell.source = source

with open('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Full_Range.ipynb', 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)
print("Saved V10_Full_Range based on V8.")
