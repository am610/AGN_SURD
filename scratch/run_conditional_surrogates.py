import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import zscore
from concurrent.futures import ProcessPoolExecutor

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

len_grid = len(cont_z)

def run_conditional_collect(X, target_idx, predictor_indices, history_idx, nlag, nbins=6):
    future_target = X[target_idx, nlag:]
    pred_1 = X[predictor_indices[0], :-nlag]
    pred_2 = X[predictor_indices[1], :-nlag]
    hist_var = X[history_idx, :-nlag]
    
    data = np.vstack([future_target, pred_1, pred_2, hist_var]).T
    hist_4d, _ = np.histogramdd(data, bins=nbins)
    hist_4d = hist_4d / np.sum(hist_4d)
    
    cond_synergy = 0.0
    
    # Calculate conditional PID as expected value over history bins
    for k in range(nbins):
        p_x3 = np.sum(hist_4d[:, :, :, k])
        if p_x3 > 1e-6:
            hist_3d = hist_4d[:, :, :, k] / p_x3
            try:
                I_R, I_S, MI, info_leak = surd.surd(hist_3d)
                syn_val = I_S.get((1, 2), 0.0)
                cond_synergy += p_x3 * syn_val
            except Exception:
                pass
                
    return cond_synergy

def get_cond_synergy_scan(X, lags, target_idx=2, predictor_indices=[0, 1], history_idx=2, nbins=6):
    syn = []
    for lag in lags:
        s = run_conditional_collect(X, target_idx, predictor_indices, history_idx, lag, nbins)
        syn.append(s)
    return np.array(syn)

def run_single_surrogate_cond(seed_val):
    np.random.seed(seed_val)
    # Generate circular shifts
    shift_c = np.random.randint(60, len_grid - 60)
    shift_b = np.random.randint(60, len_grid - 60)
    shift_y = np.random.randint(60, len_grid - 60)
    
    c_sh = np.roll(cont_z, shift_c)
    b_sh = np.roll(blue_z, shift_b)
    y_sh = np.roll(core_z, shift_y)
    
    X_sh = np.vstack([c_sh, b_sh, y_sh])
    lags_scan = np.arange(1, 121)
    
    # Target is core_z (row 2 in real X, here we shift it as y_sh)
    # predictors: c_sh, b_sh. history: y_sh
    # X_sh = [c_sh, b_sh, y_sh]
    # target_idx = 2, predictors = [0, 1], history = 2
    scan = get_cond_synergy_scan(X_sh, lags_scan, target_idx=2, predictor_indices=[0, 1], history_idx=2, nbins=6)
    return scan

if __name__ == "__main__":
    lags_scan = np.arange(1, 121)
    print("Running real conditional scan...")
    X_real = np.vstack([cont_z, blue_z, core_z])
    real_cond_scan = get_cond_synergy_scan(X_real, lags_scan, target_idx=2, predictor_indices=[0, 1], history_idx=2, nbins=6)
    
    n_surrogates = 100
    print(f"Running {n_surrogates} conditional circular shift surrogates in parallel...")
    with ProcessPoolExecutor() as executor:
        seeds = range(3000, 3000 + n_surrogates)
        results = list(executor.map(run_single_surrogate_cond, seeds))
        
    results = np.array(results)
    
    # Compute percentiles
    median_null = np.median(results, axis=0)
    p2_5_null = np.percentile(results, 2.5, axis=0)
    p97_5_null = np.percentile(results, 97.5, axis=0)
    
    # Global max-statistic significance
    real_max = np.max(real_cond_scan)
    null_maxes = np.max(results, axis=1)
    p_global = (1 + np.sum(null_maxes >= real_max)) / (1 + n_surrogates)
    
    print("\n--- CONDITIONAL SIGNIFICANCE TEST RESULTS ---")
    print(f"Real max conditional synergy: {real_max:.4f}")
    print(f"Median of null maxes:         {np.median(null_maxes):.4f}")
    print(f"Global p-value:               {p_global:.4f}")
    
    df_res = pd.DataFrame({
        'lag': lags_scan,
        'real_cond_synergy': real_cond_scan,
        'median_null_cond': median_null,
        'p2_5_null_cond': p2_5_null,
        'p97_5_null_cond': p97_5_null
    })
    output_dir = "/Users/ayan/Programs/SURD/agn_surd_project/processed/"
    df_res.to_csv(os.path.join(output_dir, "conditional_surrogate_results.csv"), index=False)
    print("Saved conditional surrogate results to processed/conditional_surrogate_results.csv")
