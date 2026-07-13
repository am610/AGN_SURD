import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore

# Configure path to SURD utilities
sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd

# Set premium plotting styles
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.titlesize': 15,
    'font.family': 'sans-serif'
})

print("Loading real NGC 5548 observation epochs...")
# Load observed MJDs and errors to match observational properties
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"
df_hb = pd.read_csv(hb_bins_path)
df_hb = df_hb[(df_hb['mjd'] >= 47512.0) & (df_hb['mjd'] <= 49255.0)].dropna()

cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd']
df_cont = df_cont[(df_cont['mjd'] >= 47512.0) & (df_cont['mjd'] <= 49255.0)].dropna()

# Real dates and errors
obs_mjd_cont = df_cont['mjd'].values
obs_mjd_line = df_hb['mjd'].values
obs_err_cont = df_cont['err'].values
obs_err_line = df_hb['core_error'].values

# Standard 1-day grid
mjd_grid = np.arange(47512, 49256)
n_grid = len(mjd_grid)

# DRW Simulator
def generate_drw_grid(n_days, tau_param=50.0, sigma_param=0.3, seed=None):
    if seed is not None:
        np.random.seed(seed)
    flux = np.zeros(n_days)
    mu = 0.0
    flux[0] = np.random.normal(mu, sigma_param)
    coeff = np.exp(-1.0 / tau_param)
    noise_std = sigma_param * np.sqrt(1 - coeff**2)
    for i in range(1, n_days):
        flux[i] = mu + coeff * (flux[i-1] - mu) + np.random.normal(0, noise_std)
    return flux

# SURD scan helper
def run_collect(X, nvars, nlag, nbins=6):
    results = {}
    for i in range(nvars):
        Y = np.vstack([X[i, nlag:], X[:, :-nlag]])
        hist, _ = np.histogramdd(Y.T, nbins)
        I_R, I_S, MI, info_leak = surd.surd(hist)
        results[i] = {"I_R": I_R, "I_S": I_S, "MI": MI, "info_leak": info_leak}
    return results

def lag_scan_target3(X, lags, nbins=6):
    metrics = {"lag": [], "S12": [], "R12": [], "U1": [], "U2": []}
    for lag in lags:
        res = run_collect(X=X, nvars=3, nlag=lag, nbins=nbins)[2]
        metrics["lag"].append(lag)
        metrics["S12"].append(res["I_S"].get((1, 2), 0.0))
        metrics["R12"].append(res["I_R"].get((1, 2), 0.0))
        metrics["U1"].append(res["I_R"].get((1,), 0.0))
        metrics["U2"].append(res["I_R"].get((2,), 0.0))
    return metrics

# Run realistic simulation over multiple Monte Carlo realizations
n_realizations = 15
lags = np.arange(1, 41)
nbins = 6

# Storage for results
case_results = {1: [], 2: [], 3: [], 4: []}

print("Running Monte Carlo realistic simulations...")
for r in range(n_realizations):
    seed = 1000 + r
    np.random.seed(seed)
    
    # 1. Generate underlying DRW continuous signals
    s1_grid = generate_drw_grid(n_grid, tau_param=50.0, sigma_param=0.3)
    
    # --- CASE 1: Single Driver ---
    s2_grid = np.random.normal(0, 0.3, n_grid)
    t3_grid = np.roll(s1_grid, 15)
    
    # --- CASE 2: Redundant Proxies ---
    s2_grid_c2 = s1_grid + np.random.normal(0, 0.1, n_grid)
    
    # --- CASE 3: Synergistic ---
    s2_grid_c3 = generate_drw_grid(n_grid, tau_param=30.0, sigma_param=0.3)
    # Strictly positive product
    s1_grid_p = s1_grid - s1_grid.min() + 0.1
    s2_grid_p = s2_grid_c3 - s2_grid_c3.min() + 0.1
    t3_grid_c3 = np.roll(s1_grid_p, 15) * np.roll(s2_grid_p, 15)
    
    # --- CASE 4: Two-Zone BLR ---
    s2_grid_c4 = np.roll(s1_grid, 10)
    t3_grid_c4 = 0.5 * np.roll(s1_grid, 10) + 0.5 * np.roll(s2_grid_c4, 20)

    # For each case, simulate observations at real epochs, add observational noise, and interpolate
    for case_num in [1, 2, 3, 4]:
        # Select active signals
        if case_num == 1:
            raw_s1, raw_s2, raw_t3 = s1_grid, s2_grid, t3_grid
        elif case_num == 2:
            raw_s1, raw_s2, raw_t3 = s1_grid, s2_grid_c2, t3_grid
        elif case_num == 3:
            raw_s1, raw_s2, raw_t3 = s1_grid_p, s2_grid_p, t3_grid_c3
        elif case_num == 4:
            raw_s1, raw_s2, raw_t3 = s1_grid, s2_grid_c4, t3_grid_c4
            
        # Sample at real MJD epochs (time mapping)
        idx_cont = np.searchsorted(mjd_grid, obs_mjd_cont)
        idx_line = np.searchsorted(mjd_grid, obs_mjd_line)
        
        # Clip index values to bounds
        idx_cont = np.clip(idx_cont, 0, n_grid-1)
        idx_line = np.clip(idx_line, 0, n_grid-1)
        
        obs_s1 = raw_s1[idx_cont] + np.random.normal(0, obs_err_cont)
        # S2 is observed on continuum epochs for Case 2/4 and line epochs for Case 1/3
        if case_num in [2, 4]:
            obs_s2 = raw_s2[idx_cont] + np.random.normal(0, obs_err_cont)
            mjd_s2 = obs_mjd_cont
        else:
            obs_s2 = raw_s2[idx_line] + np.random.normal(0, obs_err_line)
            mjd_s2 = obs_mjd_line
            
        obs_t3 = raw_t3[idx_line] + np.random.normal(0, obs_err_line)
        
        # Re-interpolate back onto 1-day grid
        s1_interp = np.interp(mjd_grid, obs_mjd_cont, obs_s1)
        s2_interp = np.interp(mjd_grid, mjd_s2, obs_s2)
        t3_interp = np.interp(mjd_grid, obs_mjd_line, obs_t3)
        
        # Standardize and run SURD
        X = np.vstack([zscore(s1_interp), zscore(s2_interp), zscore(t3_interp)])
        metrics = lag_scan_target3(X, lags, nbins=nbins)
        case_results[case_num].append(metrics)

# Compute medians and 1-sigma spreads
final_metrics = {}
for case_num in [1, 2, 3, 4]:
    r_list = case_results[case_num]
    final_metrics[case_num] = {
        'lag': lags,
        'median_S12': np.median([r['S12'] for r in r_list], axis=0),
        'p16_S12': np.percentile([r['S12'] for r in r_list], 16, axis=0),
        'p84_S12': np.percentile([r['S12'] for r in r_list], 84, axis=0),
        'median_R12': np.median([r['R12'] for r in r_list], axis=0),
        'median_U1': np.median([r['U1'] for r in r_list], axis=0),
        'median_U2': np.median([r['U2'] for r in r_list], axis=0)
    }

# Plot Summary Figure (Overwriting Figure 5)
print("Generating Figure 5: Realistic Synthetic Benchmarks with Obs Gaps/Errors...")
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

cases = [
    ('Case 1: Single Driver (True Lag = 15d)', 1, axs[0, 0]),
    ('Case 2: Redundant Proxies (True Lag = 15d)', 2, axs[0, 1]),
    ('Case 3: Synergistic Drivers (True Lag = 15d)', 3, axs[1, 0]),
    ('Case 4: Two-Zone Response (Lags = 10d, 20d)', 4, axs[1, 1])
]

for title, case_num, ax in cases:
    m = final_metrics[case_num]
    ax.plot(m['lag'], m['median_S12'], color='purple', label='Synergy (Median)', linewidth=2.5)
    ax.fill_between(m['lag'], m['p16_S12'], m['p84_S12'], color='purple', alpha=0.25, label='1$\sigma$ Spread')
    ax.plot(m['lag'], m['median_R12'], color='gray', linestyle='--', label='Redundancy', alpha=0.7)
    ax.plot(m['lag'], m['median_U1'], color='blue', linestyle=':', label='Unique (S1)', alpha=0.7)
    ax.plot(m['lag'], m['median_U2'], color='red', linestyle='-.', label='Unique (S2)', alpha=0.7)
    ax.set_xlabel('Lag (days)')
    ax.set_ylabel('Information (bits)')
    ax.set_title(title)
    ax.legend(loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
fig.savefig('overleaf_draft/figure5_synthetic_benchmarks.png', dpi=300)
plt.close(fig)
print("Saved realistic Figure 5 successfully!")
