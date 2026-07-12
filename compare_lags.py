import json

def extract_lags(nb_file):
    with open(nb_file, 'r', encoding='utf-8') as f:
        nb = json.load(f)
        
    results = {}
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            if 'prepared_data =' in "".join(cell['source']):
                for out in cell.get('outputs', []):
                    if out.get('name') == 'stdout':
                        lines = "".join(out.get('text', [])).split('\n')
                        for line in lines:
                            if 'Time range: MJD' in line or 'Original integrated Hβ' in line:
                                print(f"{nb_file} Time range: {line}")
            
            # Extract Peak Synergy
            if 'df_extended_summary' in "".join(cell['source']):
                for out in cell.get('outputs', []):
                    if out.get('data') and 'text/html' in out['data']:
                        print(f"{nb_file} Extended Summary (Synergy Lags):")
                        print(out['data']['text/plain'])
            
            # Extract ICCF
            if 'calculate_iccf' in "".join(cell['source']):
                for out in cell.get('outputs', []):
                    if out.get('name') == 'stdout':
                        print(f"{nb_file} ICCF Output:")
                        print("".join(out['text']))
                        
extract_lags('/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb')
print("----------------")
extract_lags('/Users/ayan/Programs/SURD/3_SURDS_AGN_V11_Strict_Overlap.ipynb')
