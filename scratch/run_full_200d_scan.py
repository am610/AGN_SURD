import sys
import os
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

def run_collect_2pred(target_arr, pred1_arr, pred2_arr, nlag, nbins):
    future_target = target_arr[nlag:]
    pred_1 = pred1_arr[:-nlag]
    pred_2 = pred2_arr[:-nlag]
    
    Y = np.vstack([future_target, pred_1, pred_2])
    hist, _ = np.histogramdd(Y.T, nbins)
    hist = hist / np.sum(hist)
    I_R, I_S, MI, info_leak = surd.surd(hist)
    
    joint_mi = MI.get((1, 2), 1e-14)
    if joint_mi < 1e-14:
        joint_mi = 1e-14
        
    return {
        "U1": I_R.get((1,), 0.0) / joint_mi,
        "U2": I_R.get((2,), 0.0) / joint_mi,
        "R12": I_R.get((1, 2), 0.0) / joint_mi,
        "S12": I_S.get((1, 2), 0.0) / joint_mi,
        "info_leak": info_leak
    }

lags = np.arange(1, 201)
nbins = 8

syn_core = []
syn_red = []
syn_blue = []

for l in lags:
    syn_core.append(run_collect_2pred(core_z, cont_z, blue_z, l, nbins)['S12'])
    syn_red.append(run_collect_2pred(red_z, cont_z, blue_z, l, nbins)['S12'])
    syn_blue.append(run_collect_2pred(blue_z, cont_z, core_z, l, nbins)['S12'])

print("\n--- 200-DAY SCAN PEAKS (NORMALIZED SYNERGY) ---")
c_peak = np.argmax(syn_core)
r_peak = np.argmax(syn_red)
b_peak = np.argmax(syn_blue)

print(f"Core Target: Peak Lag = {lags[c_peak]} days, Synergy = {syn_core[c_peak]:.4f}")
print(f"Red Target:  Peak Lag = {lags[r_peak]} days, Synergy = {syn_red[r_peak]:.4f}")
print(f"Blue Target: Peak Lag = {lags[b_peak]} days, Synergy = {syn_blue[b_peak]:.4f}")

# Save the curves to check how they look
df_syn = pd.DataFrame({
    'lag': lags,
    'core_synergy': syn_core,
    'red_synergy': syn_red,
    'blue_synergy': syn_blue
})
df_syn.to_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/unconditioned_synergy_200d.csv', index=False)
print("Saved 200d synergy scans to processed/unconditioned_synergy_200d.csv")
