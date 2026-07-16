import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import zscore

# Configure path to SURD utilities
sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd

# ----------------- 1. LOAD RAW DATA -----------------
print("Loading raw data...")
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

# Load Continuum
df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd_2440000', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd_2440000']

# Load Spectroscopic Bins
df_hb = pd.read_csv(hb_bins_path)

# Enforce Strict Overlap Range (MJD 47512 to 49255)
tmin = max(df_cont['mjd'].min(), df_hb['mjd'].min())
tmax = min(df_cont['mjd'].max(), df_hb['mjd'].max())
dt_final = 1.0
uniform_time_grid = np.arange(tmin, tmax + dt_final, dt_final)

# SURD Core Logic (2-predictor decomposition)
def run_collect_2pred(target_arr, pred1_arr, pred2_arr, nlag, nbins):
    future_target = target_arr[nlag:]
    pred_1 = pred1_arr[:-nlag]
    pred_2 = pred2_arr[:-nlag]
    
    Y = np.vstack([future_target, pred_1, pred_2])
    hist, _ = np.histogramdd(Y.T, nbins)
    hist = hist / np.sum(hist)
    I_R, I_S, MI, info_leak = surd.surd(hist)
    
    # Extract normalized components
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

def lag_scan_target3(target_arr, pred1_arr, pred2_arr, lags, nbins=8):
    metrics = {"lag": [], "info_leak": [], "U1": [], "U2": [], "R12": [], "S12": []}
    for lag in lags:
        res = run_collect_2pred(target_arr, pred1_arr, pred2_arr, lag, nbins)
        metrics["lag"].append(lag)
        metrics["info_leak"].append(res["info_leak"])
        metrics["U1"].append(res["U1"])
        metrics["U2"].append(res["U2"])
        metrics["R12"].append(res["R12"])
        metrics["S12"].append(res["S12"])
    return metrics

# ----------------- 2. MONTE CARLO LOOP -----------------
num_iterations = 100
lags = np.arange(1, 121)
nbins = 8

mc_results = {
    "core": {"U1": [], "U2": [], "R12": [], "S12": [], "info_leak": []},
    "red": {"U1": [], "U2": [], "R12": [], "S12": [], "info_leak": []},
    "blue": {"U1": [], "U2": [], "R12": [], "S12": [], "info_leak": []}
}

print(f"Starting {num_iterations} Monte Carlo loops...")

for it in range(num_iterations):
    if (it + 1) % 10 == 0 or it == 0:
        print(f"  Iteration {it + 1} / {num_iterations}...")
        
    # 1. Perturb fluxes by measurement errors
    pert_cont_flux = df_cont['flux'].values + np.random.normal(0, df_cont['err'].values)
    pert_blue_flux = df_hb['blue_wing_flux'].values + np.random.normal(0, df_hb['blue_wing_error'].values)
    pert_core_flux = df_hb['core_flux'].values + np.random.normal(0, df_hb['core_error'].values)
    pert_red_flux = df_hb['red_wing_flux'].values + np.random.normal(0, df_hb['red_wing_error'].values)
    
    # 2. Re-align and interpolate on the 1-day grid
    prepared_data = pd.DataFrame({'time': uniform_time_grid})
    prepared_data['cont_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_cont['mjd'], pert_cont_flux))
    prepared_data['blue_wing_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_hb['mjd'], pert_blue_flux))
    prepared_data['core_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_hb['mjd'], pert_core_flux))
    prepared_data['red_wing_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_hb['mjd'], pert_red_flux))
    prepared_data = prepared_data.dropna().reset_index(drop=True)
    
    cont_z = prepared_data['cont_flux_zscore'].values
    blue_z = prepared_data['blue_wing_flux_zscore'].values
    core_z = prepared_data['core_flux_zscore'].values
    red_z = prepared_data['red_wing_flux_zscore'].values
    
    # 3. Run SURD scans and collect normalized quantities
    # Core Target: Continuum & Blue predicting Core
    m_core = lag_scan_target3(core_z, cont_z, blue_z, lags, nbins=nbins)
    for k in ["U1", "U2", "R12", "S12", "info_leak"]:
        mc_results["core"][k].append(m_core[k])
        
    # Red Target: Continuum & Blue predicting Red
    m_red = lag_scan_target3(red_z, cont_z, blue_z, lags, nbins=nbins)
    for k in ["U1", "U2", "R12", "S12", "info_leak"]:
        mc_results["red"][k].append(m_red[k])
        
    # Blue Target: Continuum & Core predicting Blue
    m_blue = lag_scan_target3(blue_z, cont_z, core_z, lags, nbins=nbins)
    for k in ["U1", "U2", "R12", "S12", "info_leak"]:
        mc_results["blue"][k].append(m_blue[k])

# ----------------- 3. COMPUTE STATS AND SAVE -----------------
print("Computing percentiles and saving results...")

def save_mc_results(lags, mc_dict, name):
    df_data = {'lag': lags}
    for var_name, data_list in mc_dict.items():
        data_arr = np.array(data_list)
        df_data[f'median_{var_name}'] = np.median(data_arr, axis=0)
        df_data[f'p16_{var_name}'] = np.percentile(data_arr, 16, axis=0)
        df_data[f'p84_{var_name}'] = np.percentile(data_arr, 84, axis=0)
        df_data[f'p2_5_{var_name}'] = np.percentile(data_arr, 2.5, axis=0)
        df_data[f'p97_5_{var_name}'] = np.percentile(data_arr, 97.5, axis=0)
        
    df = pd.DataFrame(df_data)
    if 'median_info_leak' in df.columns:
        df = df.rename(columns={
            'median_info_leak': 'median_leak',
            'p16_info_leak': 'p16_leak',
            'p84_info_leak': 'p84_leak',
            'p2_5_info_leak': 'p2_5_leak',
            'p97_5_info_leak': 'p97_5_leak'
        })
    output_dir = "/Users/ayan/Programs/SURD/agn_surd_project/processed/"
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, f"mc_uncertainty_{name}.csv"), index=False)
    print(f"Saved: mc_uncertainty_{name}.csv")

save_mc_results(lags, mc_results["core"], "core")
save_mc_results(lags, mc_results["red"], "red")
save_mc_results(lags, mc_results["blue"], "blue")

print("Monte Carlo Uncertainty computation completed successfully!")
