import os
import sys
import numpy as np
import pandas as pd

sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd
from scipy.stats import zscore

print("Loading data for synthetic recovery analysis...")
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"
df_hb = pd.read_csv(hb_bins_path)
df_hb = df_hb[(df_hb['mjd'] >= 47512.0) & (df_hb['mjd'] <= 49255.0)].dropna()
obs_mjd_cont = df_hb['mjd'].values
obs_mjd_line = df_hb['mjd'].values
obs_err_cont = df_hb['core_error'].values * 0.1
obs_err_line = df_hb['core_error'].values * 0.1

mjd_grid = np.arange(47512, 49256)
n_grid = len(mjd_grid)

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

n_realizations = 15
lags = np.arange(1, 41)
nbins = 6

recovered_lags = {1: [], 2: [], 3: [], 4: []}

for r in range(n_realizations):
    seed = 1000 + r
    np.random.seed(seed)
    s1_grid = generate_drw_grid(n_grid, tau_param=50.0, sigma_param=0.3)
    
    # --- Case 1 ---
    s2_grid = np.random.normal(0, 0.3, n_grid)
    t3_grid = np.roll(s1_grid, 15)
    
    # --- Case 2 ---
    s2_grid_c2 = s1_grid + np.random.normal(0, 0.1, n_grid)
    
    # --- Case 3 ---
    s2_grid_c3 = generate_drw_grid(n_grid, tau_param=30.0, sigma_param=0.3)
    s1_grid_p = s1_grid - s1_grid.min() + 0.1
    s2_grid_p = s2_grid_c3 - s2_grid_c3.min() + 0.1
    t3_grid_c3 = np.roll(s1_grid_p, 15) * np.roll(s2_grid_p, 15)
    
    # --- Case 4 ---
    s2_grid_c4 = np.roll(s1_grid, 10)
    t3_grid_c4 = 0.5 * np.roll(s1_grid, 10) + 0.5 * np.roll(s2_grid_c4, 20)

    for case_num in [1, 2, 3, 4]:
        if case_num == 1:
            raw_s1, raw_s2, raw_t3 = s1_grid, s2_grid, t3_grid
        elif case_num == 2:
            raw_s1, raw_s2, raw_t3 = s1_grid, s2_grid_c2, t3_grid
        elif case_num == 3:
            raw_s1, raw_s2, raw_t3 = s1_grid_p, s2_grid_p, t3_grid_c3
        elif case_num == 4:
            raw_s1, raw_s2, raw_t3 = s1_grid, s2_grid_c4, t3_grid_c4
            
        idx_cont = np.clip(np.searchsorted(mjd_grid, obs_mjd_cont), 0, n_grid-1)
        idx_line = np.clip(np.searchsorted(mjd_grid, obs_mjd_line), 0, n_grid-1)
        
        obs_s1 = raw_s1[idx_cont] + np.random.normal(0, obs_err_cont)
        if case_num in [2, 4]:
            obs_s2 = raw_s2[idx_cont] + np.random.normal(0, obs_err_cont)
            mjd_s2 = obs_mjd_cont
        else:
            obs_s2 = raw_s2[idx_line] + np.random.normal(0, obs_err_line)
            mjd_s2 = obs_mjd_line
        obs_t3 = raw_t3[idx_line] + np.random.normal(0, obs_err_line)
        
        s1_interp = np.interp(mjd_grid, obs_mjd_cont, obs_s1)
        s2_interp = np.interp(mjd_grid, mjd_s2, obs_s2)
        t3_interp = np.interp(mjd_grid, obs_mjd_line, obs_t3)
        
        X = np.vstack([zscore(s1_interp), zscore(s2_interp), zscore(t3_interp)])
        metrics = lag_scan_target3(X, lags, nbins=nbins)
        
        if case_num == 1:
            rec_lag = lags[np.argmax(metrics['U1'])]
        elif case_num == 2:
            rec_lag = lags[np.argmax(metrics['R12'])]
        elif case_num == 3:
            rec_lag = lags[np.argmax(metrics['S12'])]
        elif case_num == 4:
            syn = np.array(metrics['S12'])
            peaks = []
            for i in range(1, len(syn)-1):
                if syn[i] > syn[i-1] and syn[i] > syn[i+1]:
                    peaks.append((syn[i], lags[i]))
            peaks = sorted(peaks, key=lambda x: x[0], reverse=True)
            if len(peaks) >= 2:
                rec_lag = sorted([peaks[0][1], peaks[1][1]])
            elif len(peaks) == 1:
                rec_lag = [peaks[0][1], np.nan]
            else:
                rec_lag = [np.nan, np.nan]
                
        recovered_lags[case_num].append(rec_lag)

# Calculate statistics
stats = {}
for c in [1, 2, 3]:
    vals = np.array(recovered_lags[c])
    # Recovery rate defined as being within +/- 5 days of true 15 days
    rec_rate = np.sum(np.abs(vals - 15.0) <= 5.0) / len(vals)
    stats[c] = {
        'mean': np.mean(vals),
        'std': np.std(vals),
        'bias': np.mean(vals) - 15.0,
        'rec_rate': rec_rate
    }

c4_vals = np.array(recovered_lags[4])
p1_vals = c4_vals[:, 0]
p2_vals = c4_vals[:, 1]
stats[4] = {
    'p1_mean': np.nanmean(p1_vals),
    'p1_std': np.nanstd(p1_vals),
    'p2_mean': np.nanmean(p2_vals),
    'p2_std': np.nanstd(p2_vals),
    'bias_p1': np.nanmean(p1_vals) - 10.0,
    'bias_p2': np.nanmean(p2_vals) - 20.0,
}

print("\n--- SYNTHETIC BENCHMARK STATISTICAL PERFORMANCE ---")
for c in [1, 2, 3]:
    print(f"Case {c}: Mean = {stats[c]['mean']:.2f} +/- {stats[c]['std']:.2f} d (Bias: {stats[c]['bias']:.2f} d, Recovery Rate: {stats[c]['rec_rate']:.2f})")
print(f"Case 4: Peak 1 = {stats[4]['p1_mean']:.2f} +/- {stats[4]['p1_std']:.2f} d, Peak 2 = {stats[4]['p2_mean']:.2f} +/- {stats[4]['p2_std']:.2f} d")

# Generate LaTeX table code
print("\nLaTeX Table Code:")
print(f"Case 1: Single Driver & 15.0 & Unique ($U_1$) & ${stats[1]['mean']:.1f} \\pm {stats[1]['std']:.1f}$ & ${stats[1]['bias']:+.1f}$ & {stats[1]['rec_rate']:.2f} \\\\")
print(f"Case 2: Redundant Proxies & 15.0 & Redundancy ($R_{{12}}$) & ${stats[2]['mean']:.1f} \\pm {stats[2]['std']:.1f}$ & ${stats[2]['bias']:+.1f}$ & {stats[2]['rec_rate']:.2f} \\\\")
print(f"Case 3: Synergistic Drivers & 15.0 & Synergy ($S_{{12}}$) & ${stats[3]['mean']:.1f} \\pm {stats[3]['std']:.1f}$ & ${stats[3]['bias']:+.1f}$ & {stats[3]['rec_rate']:.2f} \\\\")
print(f"Case 4: Two-Zone Response & 10.0 / 20.0 & Synergy ($S_{{12}}$) & ${stats[4]['p1_mean']:.1f} \\pm {stats[4]['p1_std']:.1f}$ / ${stats[4]['p2_mean']:.1f} \\pm {stats[4]['p2_std']:.1f}$ & ${stats[4]['bias_p1']:+.1f}$ / ${stats[4]['bias_p2']:+.1f}$ & -- \\\\")
