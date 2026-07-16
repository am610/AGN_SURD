import os
import sys
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import zscore

# Configure path to SURD utilities
sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd

# Load data
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

lags_200 = np.arange(1, 201)
nbins = 8

def run_collect_2pred_norm(target_arr, pred1_arr, pred2_arr, nlag, nbins):
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
    return I_S.get((1, 2), 0.0) / joint_mi

def get_max_normalized_synergy_scan(target_arr, pred1_arr, pred2_arr, lags, nbins=8):
    syn = []
    for lag in lags:
        s = run_collect_2pred_norm(target_arr, pred1_arr, pred2_arr, lag, nbins)
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
    
    # Core target (T=y_sh, S1=c_sh, S2=b_sh)
    max_c = get_max_normalized_synergy_scan(y_sh, c_sh, b_sh, lags_200, nbins=nbins)
    
    # Red target (T=r_sh, S1=c_sh, S2=b_sh)
    max_r = get_max_normalized_synergy_scan(r_sh, c_sh, b_sh, lags_200, nbins=nbins)
    
    # Blue target (T=b_sh, S1=c_sh, S2=y_sh)
    max_b = get_max_normalized_synergy_scan(b_sh, c_sh, y_sh, lags_200, nbins=nbins)
    
    return max_c, max_r, max_b

if __name__ == "__main__":
    print("Computing real max normalized synergies up to 200 days...")
    real_max_c = get_max_normalized_synergy_scan(y_z, c_z, b_z, lags_200, nbins=nbins)
    real_max_r = get_max_normalized_synergy_scan(r_z, c_z, b_z, lags_200, nbins=nbins)
    real_max_b = get_max_normalized_synergy_scan(b_z, c_z, y_z, lags_200, nbins=nbins)
    
    print(f"Real max normalized synergies: Core={real_max_c:.4f}, Red={real_max_r:.4f}, Blue={real_max_b:.4f}")
    
    n_surrogates = 100
    print(f"Running {n_surrogates} circular shift normalized surrogates in parallel...")
    with ProcessPoolExecutor() as executor:
        seeds = range(4000, 4000 + n_surrogates)
        results = list(executor.map(run_single_surrogate, seeds))
        
    null_max_c = [res[0] for res in results]
    null_max_r = [res[1] for res in results]
    null_max_b = [res[2] for res in results]
        
    p_c = (1 + np.sum(np.array(null_max_c) >= real_max_c)) / (1 + n_surrogates)
    p_r = (1 + np.sum(np.array(null_max_r) >= real_max_r)) / (1 + n_surrogates)
    p_b = (1 + np.sum(np.array(null_max_b) >= real_max_b)) / (1 + n_surrogates)
    
    print("\n--- GLOBAL NORMALIZED SIGNIFICANCE RESULTS (200 DAYS) ---")
    print(f"Core Target: real_max={real_max_c:.4f}, median_null_max={np.median(null_max_c):.4f}, p_global={p_c:.4f}")
    print(f"Red Target:  real_max={real_max_r:.4f}, median_null_max={np.median(null_max_r):.4f}, p_global={p_r:.4f}")
    print(f"Blue Target: real_max={real_max_b:.4f}, median_null_max={np.median(null_max_b):.4f}, p_global={p_b:.4f}")
    
    pd.DataFrame({
        'target': ['core', 'red', 'blue'],
        'real_max': [real_max_c, real_max_r, real_max_b],
        'median_null_max': [np.median(null_max_c), np.median(null_max_r), np.median(null_max_b)],
        'p_global': [p_c, p_r, p_b]
    }).to_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/global_max_statistic_normalized_results.csv', index=False)
    print("Saved global max-statistic results to CSV.")
