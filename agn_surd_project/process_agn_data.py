import pandas as pd
import numpy as np
import os
from astropy.io import fits

def load_agnwatch_ngc5548():
    # Continuum 5100A
    c5100_path = "/home/ubuntu/agn_data/ngc5548_agnwatch/c5100.dat"
    df_cont = pd.read_csv(c5100_path, sep='\s+', names=['jd_2440000', 'flux', 'error'], comment='#')
    df_cont['mjd'] = df_cont['jd_2440000'] + 2440000 - 2400000.5
    
    # H-beta
    hb_path = "/home/ubuntu/agn_data/ngc5548_agnwatch/hb.dat"
    df_hb = pd.read_csv(hb_path, sep='\s+', names=['jd_2440000', 'flux', 'error'], comment='#')
    df_hb['mjd'] = df_hb['jd_2440000'] + 2440000 - 2400000.5
    
    return df_cont, df_hb

def load_storm1_ngc5548():
    path = "/home/ubuntu/agn_data/ngc5548_storm1/emission_lines.txt"
    # Fixed-width format starting from line 31
    # Column definitions from the header
    col_names = [
        'thjd', 'f1158', 'e_f1158', 'f1367', 'e_f1367', 'f1469', 'e_f1469', 'f1745', 'e_f1745',
        'f_lya', 'e_f_lya', 'f_nv', 'e_f_nv', 'f_siiv', 'e_f_siiv', 'f_civ', 'e_f_civ', 'f_heii', 'e_f_heii'
    ]
    # Widths are 10 characters each
    widths = [10] * 19
    df = pd.read_fwf(path, widths=widths, names=col_names, skiprows=30)
    # THJD is HJD - 2440000
    df['mjd'] = df['thjd'] + 2440000 - 2400000.5
    return df

def load_storm2_mrk817():
    lc_file = "/home/ubuntu/agn_data/mrk817_storm2/hlsp_storm2_hst_cos_mrk817-go16196_g130m-g160m_v1_lightcurve.fits"
    if not os.path.exists(lc_file):
        for root, dirs, files in os.walk("/home/ubuntu/agn_data/mrk817_storm2"):
            for file in files:
                if "lightcurve.fits" in file:
                    lc_file = os.path.join(root, file)
                    break
    
    data_dict = {}
    with fits.open(lc_file) as hdul:
        for i in range(1, len(hdul)):
            name = hdul[i].name
            data = hdul[i].data
            df = pd.DataFrame(data)
            df.columns = [c.lower() for c in df.columns]
            data_dict[name] = df
            
    return data_dict

if __name__ == "__main__":
    print("Processing NGC 5548 AGN Watch...")
    df_cont, df_hb = load_agnwatch_ngc5548()
    df_cont.to_csv("/home/ubuntu/agn_data/ngc5548_agnwatch_cont_clean.csv", index=False)
    df_hb.to_csv("/home/ubuntu/agn_data/ngc5548_agnwatch_hb_clean.csv", index=False)
    
    print("Processing NGC 5548 STORM 1...")
    df_storm1 = load_storm1_ngc5548()
    df_storm1.to_csv("/home/ubuntu/agn_data/ngc5548_storm1_clean.csv", index=False)
    
    print("Processing Mrk 817 STORM 2...")
    storm2_data = load_storm2_mrk817()
    for name, df in storm2_data.items():
        df.to_csv(f"/home/ubuntu/agn_data/mrk817_storm2_{name.lower()}_clean.csv", index=False)
    
    print("Done.")
