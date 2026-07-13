import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configure path to SURD utilities
sys.path.append("/Users/ayan/Programs/SURD/SURD/utils")
sys.path.append("/Users/ayan/Programs/SURD/SURD")
import surd

# Load dataset
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"
df_hb = pd.read_csv(hb_bins_path)
df_hb = df_hb[(df_hb['mjd'] >= 47512.0) & (df_hb['mjd'] <= 49255.0)].dropna()

cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd']
df_cont = df_cont[(df_cont['mjd'] >= 47512.0) & (df_cont['mjd'] <= 49255.0)].dropna()

# Grid alignment (1-day resolution)
mjd_grid = np.arange(47512, 49256)
cont_grid = np.interp(mjd_grid, df_cont['mjd'], df_cont['flux'])
blue_grid = np.interp(mjd_grid, df_hb['mjd'], df_hb['blue_wing_flux'])
core_grid = np.interp(mjd_grid, df_hb['mjd'], df_hb['core_flux'])
red_grid = np.interp(mjd_grid, df_hb['mjd'], df_hb['red_wing_flux'])

# Standardize
def zscore(x):
    return (x - np.mean(x)) / np.std(x)

c_z = zscore(cont_grid)
b_z = zscore(blue_grid)
y_z = zscore(core_grid)
r_z = zscore(red_grid)

# Stack for all three targets
# Target indices: Core=2, Red=3, Blue=1
# Predictor combinations for each:
# - Core (idx 2) predicted by Continuum (0) and Blue (1), conditioned on Core (2)
# - Red (idx 3) predicted by Continuum (0) and Blue (1), conditioned on Red (3)
# - Blue (idx 1) predicted by Continuum (0) and Core (2), conditioned on Blue (1)

X = np.vstack([c_z, b_z, y_z, r_z])

def run_conditional_collect(X, target_idx, predictor_indices, history_idx, nlag, nbins=6):
    future_target = X[target_idx, nlag:]
    pred_1 = X[predictor_indices[0], :-nlag]
    pred_2 = X[predictor_indices[1], :-nlag]
    hist_var = X[history_idx, :-nlag]
    
    data = np.vstack([future_target, pred_1, pred_2, hist_var]).T
    hist_4d, _ = np.histogramdd(data, bins=nbins)
    hist_4d = hist_4d / np.sum(hist_4d)
    
    cond_synergy = 0.0
    cond_leak = 0.0
    
    for k in range(nbins):
        p_x3 = np.sum(hist_4d[:, :, :, k])
        if p_x3 > 1e-6:
            hist_3d = hist_4d[:, :, :, k] / p_x3
            try:
                I_R, I_S, MI, info_leak = surd.surd(hist_3d)
                syn_val = I_S.get((1, 2), 0.0)
                cond_synergy += p_x3 * syn_val
                cond_leak += p_x3 * info_leak
            except Exception:
                pass
                
    return cond_synergy, cond_leak

def run_unconditioned_collect(X, target_idx, predictor_indices, nlag, nbins=6):
    future_target = X[target_idx, nlag:]
    pred_1 = X[predictor_indices[0], :-nlag]
    pred_2 = X[predictor_indices[1], :-nlag]
    
    data = np.vstack([future_target, pred_1, pred_2]).T
    hist_3d, _ = np.histogramdd(data, bins=nbins)
    hist_3d = hist_3d / np.sum(hist_3d)
    
    I_R, I_S, MI, info_leak = surd.surd(hist_3d)
    return I_S.get((1, 2), 0.0), info_leak

lags = np.arange(1, 121)

# Scan for Core
core_uncond_syn, core_uncond_leak = [], []
core_cond_syn, core_cond_leak = [], []
for lag in lags:
    us, ul = run_unconditioned_collect(X, 2, [0, 1], lag, nbins=6)
    cs, cl = run_conditional_collect(X, 2, [0, 1], 2, lag, nbins=6)
    core_uncond_syn.append(us)
    core_uncond_leak.append(ul)
    core_cond_syn.append(cs)
    core_cond_leak.append(cl)

print("Core scan completed.")
print(f"Unconditioned Core Max Synergy: {max(core_uncond_syn):.4f} at lag {lags[np.argmax(core_uncond_syn)]} days")
print(f"Conditioned Core Max Synergy: {max(core_cond_syn):.4f} at lag {lags[np.argmax(core_cond_syn)]} days")

# Plot Core Comparison
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(lags, core_uncond_syn, label='Unconditioned Synergy', color='blue', linewidth=2)
plt.plot(lags, core_cond_syn, label='History-Conditioned Synergy', color='red', linewidth=2)
plt.xlabel('Lag (days)')
plt.ylabel('Synergy (bits)')
plt.title('Synergy Comparison: Core Target')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)

plt.subplot(1, 2, 2)
plt.plot(lags, core_uncond_leak, label='Unconditioned Leak', color='blue', linewidth=2)
plt.plot(lags, core_cond_leak, label='History-Conditioned Leak', color='red', linewidth=2)
plt.xlabel('Lag (days)')
plt.ylabel('Information Leak (normalized entropy)')
plt.title('Information Leak Comparison: Core Target')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
os.makedirs("/Users/ayan/Programs/SURD/overleaf_draft", exist_ok=True)
plt.savefig("/Users/ayan/Programs/SURD/overleaf_draft/figure7_history_conditioning.png", dpi=300)
plt.close()
print("Saved comparison figure to overleaf_draft/figure7_history_conditioning.png")
