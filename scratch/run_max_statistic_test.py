import os
import sys
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# Configure path to SURD utilities
sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd
from scipy.stats import zscore

print("Loading data for global max-statistic test (4D)...")
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd', 'flux', 'err'])
df_hb = pd.read_csv(hb_bins_path)

tmin = max(df_cont['jd'].min(), df_hb['mjd'].min())
tmax = min(df_cont['jd'].max(), df_hb['mjd'].max())
mjd_grid = np.arange(tmin, tmax + 1.0, 1.0)
len_grid = len(mjd_grid)

c_z = zscore(np.interp(mjd_grid, df_cont['jd'], df_cont['flux']))
b_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['blue_wing_flux']))
y_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['core_flux']))
r_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['red_wing_flux']))

lags_120 = np.arange(1, 121)
nbins = 8

def run_collect_4d(X, target_idx, nlag, nbins):
    Y = np.vstack([X[target_idx, nlag:], X[:, :-nlag]])
    hist, _ = np.histogramdd(Y.T, nbins)
    I_R, I_S, MI, info_leak = surd.surd(hist)
    # The synergy term of interest is for the other two predictors.
    # In Y, row 0 is target. Row 1 is X[0], row 2 is X[1], row 3 is X[2].
    # We want the synergy between the two non-target predictors.
    # Let's map target_idx to find the indices of the other two:
    predictor_indices = [idx + 1 for idx in range(3) if idx != target_idx]
    return I_S.get(tuple(predictor_indices), 0.0)

def get_max_synergy_scan(X, lags, target_idx=2, nbins=8):
    syn = []
    for lag in lags:
        s = run_collect_4d(X, target_idx, lag, nbins)
        syn.append(s)
    return np.max(syn)

def run_single_surrogate(seed_val):
    np.random.seed(seed_val)
    shift_c = np.random.randint(60, len_grid - 60)
    shift_b = np.random.randint(60, len_grid - 60)
    shift_y = np.random.randint(60, len_grid - 60)
    shift_r = np.random.randint(60, len_grid - 60)
    
    c_sh = np.roll(c_z, shift_c)
    b_sh = np.roll(b_z, shift_b)
    y_sh = np.roll(y_z, shift_y)
    r_sh = np.roll(r_z, shift_r)
    
    # Core target (T=y_z, predictors at t = c_sh, b_sh, y_sh)
    # Note: the order of variables in X must match the real scans
    # For Core target, X = stacked [c_z, b_z, y_z]
    X_c_sh = np.vstack([c_sh, b_sh, y_sh])
    max_c = get_max_synergy_scan(X_c_sh, lags_120, target_idx=2, nbins=nbins)
    
    # Red target, X = stacked [c_z, b_z, r_z]
    X_r_sh = np.vstack([c_sh, b_sh, r_sh])
    max_r = get_max_synergy_scan(X_r_sh, lags_120, target_idx=2, nbins=nbins)
    
    # Blue target, X = stacked [c_z, y_z, b_z]
    X_b_sh = np.vstack([c_sh, y_sh, b_sh])
    max_b = get_max_synergy_scan(X_b_sh, lags_120, target_idx=2, nbins=nbins)
    
    return max_c, max_r, max_b

if __name__ == "__main__":
    print("Computing real max synergies...")
    real_max_c = get_max_synergy_scan(np.vstack([c_z, b_z, y_z]), lags_120, target_idx=2, nbins=nbins)
    real_max_r = get_max_synergy_scan(np.vstack([c_z, b_z, r_z]), lags_120, target_idx=2, nbins=nbins)
    real_max_b = get_max_synergy_scan(np.vstack([c_z, y_z, b_z]), lags_120, target_idx=2, nbins=nbins)
    
    print(f"Real max synergies: Core={real_max_c:.4f}, Red={real_max_r:.4f}, Blue={real_max_b:.4f}")
    
    n_surrogates = 500
    print(f"Running {n_surrogates} circular shift surrogate scans in parallel...")
    
    # Use ProcessPoolExecutor to parallelize
    with ProcessPoolExecutor() as executor:
        seeds = range(1000, 1000 + n_surrogates)
        results = list(executor.map(run_single_surrogate, seeds))
        
    null_max_c = [res[0] for res in results]
    null_max_r = [res[1] for res in results]
    null_max_b = [res[2] for res in results]
        
    # Calculate global p-values with standard formula: (1 + count(null >= real)) / (1 + N)
    p_c = (1 + np.sum(np.array(null_max_c) >= real_max_c)) / (1 + n_surrogates)
    p_r = (1 + np.sum(np.array(null_max_r) >= real_max_r)) / (1 + n_surrogates)
    p_b = (1 + np.sum(np.array(null_max_b) >= real_max_b)) / (1 + n_surrogates)
    
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
