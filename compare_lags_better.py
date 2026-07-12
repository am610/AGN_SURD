import json

def extract_lags(nb_file):
    with open(nb_file, 'r', encoding='utf-8') as f:
        nb = json.load(f)
        
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            # Extract Peak Synergy (from the table display)
            if 'df_extended_summary' in "".join(cell['source']):
                for out in cell.get('outputs', []):
                    if out.get('data') and 'text/plain' in out['data']:
                        print(f"=== {nb_file} Extended Summary ===")
                        print("".join(out['data']['text/plain']))
            
            # Extract ICCF
            if 'print(f"ICCF Peak Lag' in "".join(cell['source']):
                for out in cell.get('outputs', []):
                    if out.get('name') == 'stdout':
                        print(f"=== {nb_file} ICCF ===")
                        print("".join(out['text']))

extract_lags('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb')
print("----------------")
extract_lags('/Users/ayan/Programs/SURD/3_SURDS_AGN_V11_Strict_Overlap.ipynb')
