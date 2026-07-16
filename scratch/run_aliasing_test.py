import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import zscore

sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd

# Load observed MJDs to match seasonal gaps
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"
df_hb = pd.read_csv(hb_bins_path)
df_hb = df_hb[(df_hb['mjd'] >= 47512.0) & (df_hb['mjd'] <= 49255.0)].dropna()

cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd']
df_cont = df_cont[(df_cont['mjd'] >= 47512.0) & (df_cont['mjd'] <= 49255.0)].dropna()

obs_mjd_cont = df_cont['mjd'].values
obs_mjd_line = df_hb['mjd'].values
obs_err_cont = df_cont['err'].values
obs_err_line = df_hb['core_error'].values

mjd_grid = np.arange(47512, 49256)
n_grid = len(mjd_grid)

# DRW Simulator
def generate_drw_grid(n_days, tau_param=50.0, sigma_param=0.3):
    flux = np.zeros(n_days)
    mu = 0.0
    flux[0] = np.random.normal(mu, sigma_param)
    coeff = np.exp(-1.0 / tau_param)
    noise_std = sigma_param * np.sqrt(1 - coeff**2)
    for i in range(1, n_days):
        flux[i] = mu + coeff * (flux[i-1] - mu) + np.random.normal(0, noise_std)
    return flux

def run_collect_2pred(target_arr, pred1_arr, pred2_arr, nlag, nbins=6):
    future_target = target_arr[nlag:]
    pred_1 = pred1_arr[:-nlag]
    pred_2 = pred2_arr[:-nlag]
    Y = np.vstack([future_target, pred_1, pred_2])
    hist, _ = np.histogramdd(Y.T, nbins)
    hist = hist / np.sum(hist)
    I_R, I_S, MI, info_leak = surd.surd(hist)
    joint_mi = MI.get((1, 2), 1e-14)
    if joint_mi < 1e-14: joint_mi = 1e-14
    return I_S.get((1, 2), 0.0) / joint_mi

n_realizations = 50
lags = np.arange(1, 201)
nbins = 6

# Store synergy scans for negative control
# Case: Only short-lag (15d) coupling, scanned up to 200d.
all_syn_scans = []

print("Running Monte Carlo seasonal aliasing test (negative control)...")
for r in range(n_realizations):
    if (r + 1) % 10 == 0 or r == 0:
        print(f"  Realization {r + 1} / {n_realizations}...")
    np.random.seed(2000 + r)
    
    # 1. Generate underlying DRW continuous signals
    s1_grid = generate_drw_grid(n_grid, tau_param=50.0, sigma_param=0.3)
    s2_grid = generate_drw_grid(n_grid, tau_param=30.0, sigma_param=0.3)
    
    # Target has ONLY short-lag coupling (15d) with s1 and s2 (linear combination + noise)
    # Zero long-lag coupling exists
    t3_grid = 0.5 * np.roll(s1_grid, 15) + 0.5 * np.roll(s2_grid, 15) + np.random.normal(0, 0.05, n_grid)
    
    # 2. Sample at real epochs with noise
    idx_cont = np.clip(np.searchsorted(mjd_grid, obs_mjd_cont), 0, n_grid-1)
    idx_line = np.clip(np.searchsorted(mjd_grid, obs_mjd_line), 0, n_grid-1)
    
    obs_s1 = s1_grid[idx_cont] + np.random.normal(0, obs_err_cont)
    obs_s2 = s2_grid[idx_line] + np.random.normal(0, obs_err_line)
    obs_t3 = t3_grid[idx_line] + np.random.normal(0, obs_err_line)
    
    # 3. Interpolate back to 1-day grid
    s1_interp = np.interp(mjd_grid, obs_mjd_cont, obs_s1)
    s2_interp = np.interp(mjd_grid, obs_mjd_line, obs_s2)
    t3_interp = np.interp(mjd_grid, obs_mjd_line, obs_t3)
    
    # 4. Standardize
    s1_z = zscore(s1_interp)
    s2_z = zscore(s2_interp)
    t3_z = zscore(t3_interp)
    
    # 5. Run SURD scan up to 200d
    syn_scan = [run_collect_2pred(t3_z, s1_z, s2_z, l, nbins=nbins) for l in lags]
    all_syn_scans.append(syn_scan)

# Compute medians and percentiles
all_syn_scans = np.array(all_syn_scans)
median_syn = np.median(all_syn_scans, axis=0)
p16_syn = np.percentile(all_syn_scans, 16, axis=0)
p84_syn = np.percentile(all_syn_scans, 84, axis=0)
p2_5_syn = np.percentile(all_syn_scans, 2.5, axis=0)
p97_5_syn = np.percentile(all_syn_scans, 97.5, axis=0)

df_alias = pd.DataFrame({
    'lag': lags,
    'median_syn': median_syn,
    'p16_syn': p16_syn,
    'p84_syn': p84_syn,
    'p2_5_syn': p2_5_syn,
    'p97_5_syn': p97_5_syn
})

output_dir = "/Users/ayan/Programs/SURD/agn_surd_project/processed/"
df_alias.to_csv(os.path.join(output_dir, "seasonal_aliasing_null_test.csv"), index=False)
print("Saved seasonal aliasing null test results to processed/seasonal_aliasing_null_test.csv")
