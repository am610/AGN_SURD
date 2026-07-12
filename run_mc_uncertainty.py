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

# SURD Core Logic
def run_collect(X, nvars, nlag, nbins):
    results = {}
    for i in range(nvars):
        Y = np.vstack([X[i, nlag:], X[:, :-nlag]])
        hist, _ = np.histogramdd(Y.T, nbins)
        I_R, I_S, MI, info_leak = surd.surd(hist)
        results[i] = {"I_R": I_R, "I_S": I_S, "MI": MI, "info_leak": info_leak}
    return results

def lag_scan_target3(X, lags, nbins=8):
    metrics = {"lag": [], "info_leak": [], "S12": []}
    for lag in lags:
        res = run_collect(X=X, nvars=3, nlag=lag, nbins=nbins)[2]
        metrics["lag"].append(lag)
        metrics["info_leak"].append(res["info_leak"])
        metrics["S12"].append(res["I_S"].get((1, 2), np.nan))
    return metrics

# ----------------- 2. MONTE CARLO LOOP -----------------
num_iterations = 100
lags = np.arange(1, 121)
nbins = 8

mc_synergy_core = []
mc_leak_core = []
mc_synergy_red = []
mc_leak_red = []
mc_synergy_blue = []
mc_leak_blue = []

print(f"Starting {num_iterations} Monte Carlo loops...")

for it in range(num_iterations):
    if (it + 1) % 10 == 0 or it == 0:
        print(f"  Iteration {it + 1} / {num_iterations}...")
        
    # 1. Perturb fluxes by measurement errors
    # Continuum
    pert_cont_flux = df_cont['flux'].values + np.random.normal(0, df_cont['err'].values)
    # Spectroscopic bins
    pert_blue_flux = df_hb['blue_wing_flux'].values + np.random.normal(0, df_hb['blue_wing_error'].values)
    pert_core_flux = df_hb['core_flux'].values + np.random.normal(0, df_hb['core_error'].values)
    pert_red_flux = df_hb['red_wing_flux'].values + np.random.normal(0, df_hb['red_wing_error'].values)
    
    # 2. Re-align and interpolate on the 1-day grid
    prepared_data = pd.DataFrame({'time': uniform_time_grid})
    
    # Continuum interpolation
    prepared_data['cont_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_cont['mjd'], pert_cont_flux))
    
    # Wings interpolation
    prepared_data['blue_wing_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_hb['mjd'], pert_blue_flux))
    prepared_data['core_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_hb['mjd'], pert_core_flux))
    prepared_data['red_wing_flux_zscore'] = zscore(np.interp(uniform_time_grid, df_hb['mjd'], pert_red_flux))
    
    prepared_data = prepared_data.dropna().reset_index(drop=True)
    
    # Extract arrays
    cont_z = prepared_data['cont_flux_zscore'].values
    blue_z = prepared_data['blue_wing_flux_zscore'].values
    core_z = prepared_data['core_flux_zscore'].values
    red_z = prepared_data['red_wing_flux_zscore'].values
    
    # 3. Run SURD scans
    # Core Target
    X_core = np.vstack([cont_z, blue_z, core_z])
    m_core = lag_scan_target3(X_core, lags, nbins=nbins)
    mc_synergy_core.append(m_core['S12'])
    mc_leak_core.append(m_core['info_leak'])
    
    # Red Target
    X_red = np.vstack([cont_z, blue_z, red_z])
    m_red = lag_scan_target3(X_red, lags, nbins=nbins)
    mc_synergy_red.append(m_red['S12'])
    mc_leak_red.append(m_red['info_leak'])
    
    # Blue Target
    X_blue = np.vstack([cont_z, core_z, blue_z])
    m_blue = lag_scan_target3(X_blue, lags, nbins=nbins)
    mc_synergy_blue.append(m_blue['S12'])
    mc_leak_blue.append(m_blue['info_leak'])

# ----------------- 3. COMPUTE STATS AND SAVE -----------------
print("Computing percentiles and saving results...")

def save_mc_results(lags, mc_synergy, mc_leak, name):
    mc_synergy = np.array(mc_synergy)
    mc_leak = np.array(mc_leak)
    
    df = pd.DataFrame({
        'lag': lags,
        'median_synergy': np.median(mc_synergy, axis=0),
        'p16_synergy': np.percentile(mc_synergy, 16, axis=0),
        'p84_synergy': np.percentile(mc_synergy, 84, axis=0),
        'p2_5_synergy': np.percentile(mc_synergy, 2.5, axis=0),
        'p97_5_synergy': np.percentile(mc_synergy, 97.5, axis=0),
        'median_leak': np.median(mc_leak, axis=0),
        'p16_leak': np.percentile(mc_leak, 16, axis=0),
        'p84_leak': np.percentile(mc_leak, 84, axis=0),
        'p2_5_leak': np.percentile(mc_leak, 2.5, axis=0),
        'p97_5_leak': np.percentile(mc_leak, 97.5, axis=0)
    })
    
    output_dir = "/Users/ayan/Programs/SURD/agn_surd_project/processed/"
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, f"mc_uncertainty_{name}.csv"), index=False)
    print(f"Saved: mc_uncertainty_{name}.csv")

save_mc_results(lags, mc_synergy_core, mc_leak_core, "core")
save_mc_results(lags, mc_synergy_red, mc_leak_red, "red")
save_mc_results(lags, mc_synergy_blue, mc_leak_blue, "blue")

print("Monte Carlo Uncertainty computation completed successfully!")
