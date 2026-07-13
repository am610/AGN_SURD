import os
import sys
import numpy as np
import pandas as pd

# Configure path to SURD utilities
sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd
from scipy.stats import zscore

# Load aligned grid data from produce_paper_plots.py environment
# To be fast, we will recreate the grid here
print("Loading data for global max-statistic test...")
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd']
df_hb = pd.read_csv(hb_bins_path)

tmin = max(df_cont['mjd'].min(), df_hb['mjd'].min())
tmax = min(df_cont['mjd'].max(), df_hb['mjd'].max())
mjd_grid = np.arange(tmin, tmax + 1.0, 1.0)

c_z = zscore(np.interp(mjd_grid, df_cont['mjd'], df_cont['flux']))
b_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['blue_wing_flux']))
y_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['core_flux']))
r_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['red_wing_flux']))

lags_120 = np.arange(1, 121)
nbins = 8

def run_collect(X, nvars, nlag, nbins):
    results = {}
    for i in range(nvars):
        Y = np.vstack([X[i, nlag:], X[:, :-nlag]])
        hist, _ = np.histogramdd(Y.T, nbins)
        I_R, I_S, MI, info_leak = surd.surd(hist)
        results[i] = {"I_R": I_R, "I_S": I_S, "MI": MI, "info_leak": info_leak}
    return results

def get_max_synergy_scan(X, lags, nbins=8):
    syn = []
    for lag in lags:
        res = run_collect(X, 3, lag, nbins)[2]
        syn.append(res["I_S"].get((1, 2), 0.0))
    return np.max(syn)

# Real maximum synergy values
print("Computing real max synergies...")
real_max_c = get_max_synergy_scan(np.vstack([c_z, b_z, y_z]), lags_120, nbins=nbins)
real_max_r = get_max_synergy_scan(np.vstack([c_z, b_z, r_z]), lags_120, nbins=nbins)
real_max_b = get_max_synergy_scan(np.vstack([c_z, y_z, b_z]), lags_120, nbins=nbins)

print(f"Real max synergies: Core={real_max_c:.4f}, Red={real_max_r:.4f}, Blue={real_max_b:.4f}")

# Run surrogate realizations
n_surrogates = 30
null_max_c = []
null_max_r = []
null_max_b = []

print(f"Running {n_surrogates} circular shift surrogate scans...")
for s in range(n_surrogates):
    # Circularly shift predictors by random amount between 60 and len(mjd_grid)-60
    shift_c = np.random.randint(60, len(mjd_grid) - 60)
    shift_b = np.random.randint(60, len(mjd_grid) - 60)
    shift_y = np.random.randint(60, len(mjd_grid) - 60)
    shift_r = np.random.randint(60, len(mjd_grid) - 60)
    
    c_sh = np.roll(c_z, shift_c)
    b_sh = np.roll(b_z, shift_b)
    y_sh = np.roll(y_z, shift_y)
    r_sh = np.roll(r_z, shift_r)
    
    # Core target (T3=y, predictors=c, b)
    X_c_sh = np.vstack([c_sh, b_sh, y_z])
    null_max_c.append(get_max_synergy_scan(X_c_sh, lags_120, nbins=nbins))
    
    # Red target (T3=r, predictors=c, b)
    X_r_sh = np.vstack([c_sh, b_sh, r_z])
    null_max_r.append(get_max_synergy_scan(X_r_sh, lags_120, nbins=nbins))
    
    # Blue target (T3=b, predictors=c, y)
    X_b_sh = np.vstack([c_sh, y_sh, b_z])
    null_max_b.append(get_max_synergy_scan(X_b_sh, lags_120, nbins=nbins))
    
    print(f"  Finished surrogate scan {s+1}/{n_surrogates}")

# Calculate global p-values
p_c = np.sum(np.array(null_max_c) >= real_max_c) / n_surrogates
p_r = np.sum(np.array(null_max_r) >= real_max_r) / n_surrogates
p_b = np.sum(np.array(null_max_b) >= real_max_b) / n_surrogates

print("\n--- GLOBAL SIGNIFICANCE TEST RESULTS ---")
print(f"Core Target: real_max={real_max_c:.4f}, median_null_max={np.median(null_max_c):.4f}, p_global={p_c:.4f}")
print(f"Red Target:  real_max={real_max_r:.4f}, median_null_max={np.median(null_max_r):.4f}, p_global={p_r:.4f}")
print(f"Blue Target: real_max={real_max_b:.4f}, median_null_max={np.median(null_max_b):.4f}, p_global={p_b:.4f}")

# Save results to file
pd.DataFrame({
    'target': ['core', 'red', 'blue'],
    'real_max': [real_max_c, real_max_r, real_max_b],
    'median_null_max': [np.median(null_max_c), np.median(null_max_r), np.median(null_max_b)],
    'p_global': [p_c, p_r, p_b]
}).to_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/global_max_statistic_results.csv', index=False)
print("Saved global max-statistic results to CSV.")
