import sys
import numpy as np
import pandas as pd
from scipy.stats import zscore

sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd

# Load raw data
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd_2440000', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd_2440000']
df_hb = pd.read_csv(hb_bins_path)

tmin = max(df_cont['mjd'].min(), df_hb['mjd'].min())
tmax = min(df_cont['mjd'].max(), df_hb['mjd'].max())
dt_final = 1.0
uniform_time_grid = np.arange(tmin, tmax + dt_final, dt_final)

prepared_data = pd.DataFrame({'time': uniform_time_grid})
valid_cont = df_cont.dropna(subset=['flux']).sort_values('mjd')
prepared_data['cont_flux_zscore'] = zscore(np.interp(uniform_time_grid, valid_cont['mjd'], valid_cont['flux']))

for col, new_name in [('blue_wing_flux', 'blue_wing_flux_zscore'), 
                      ('core_flux', 'core_flux_zscore'), 
                      ('red_wing_flux', 'red_wing_flux_zscore')]:
    valid_data = df_hb.dropna(subset=[col]).sort_values('mjd')
    prepared_data[new_name] = zscore(np.interp(uniform_time_grid, valid_data['mjd'], valid_data[col]))

prepared_data = prepared_data.dropna().reset_index(drop=True)
cont_z = prepared_data['cont_flux_zscore'].values
blue_z = prepared_data['blue_wing_flux_zscore'].values
core_z = prepared_data['core_flux_zscore'].values
red_z = prepared_data['red_wing_flux_zscore'].values

def check_lags(target_arr, pred1_arr, pred2_arr, lags, nbins=8):
    results = []
    for lag in lags:
        future_target = target_arr[lag:]
        pred_1 = pred1_arr[:-lag]
        pred_2 = pred2_arr[:-lag]
        
        Y = np.vstack([future_target, pred_1, pred_2])
        hist, _ = np.histogramdd(Y.T, nbins)
        hist = hist / np.sum(hist)
        I_R, I_S, MI, info_leak = surd.surd(hist)
        
        raw_syn = I_S.get((1, 2), 0.0)
        joint_mi = MI.get((1, 2), 1e-14)
        norm_syn = raw_syn / joint_mi if joint_mi > 1e-14 else 0.0
        
        results.append({
            'lag': lag,
            'raw_syn': raw_syn,
            'joint_mi': joint_mi,
            'norm_syn': norm_syn
        })
    return pd.DataFrame(results)

print("Checking Red Wing Hbeta...")
df_red = check_lags(red_z, cont_z, blue_z, np.arange(1, 201))

print("\nTop 5 lags for raw synergy (unnormalized):")
print(df_red.sort_values('raw_syn', ascending=False).head(5))

print("\nTop 5 lags for normalized synergy:")
print(df_red.sort_values('norm_syn', ascending=False).head(5))
