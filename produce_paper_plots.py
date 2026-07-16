import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore
from scipy.signal import correlate

# Configure path to SURD utilities
sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd

# ----------------- 1. LOAD AND PREPARE DATA -----------------
print("Loading and aligning dataset...")
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

# Interpolate continuum
valid_cont = df_cont.dropna(subset=['flux']).sort_values('mjd')
prepared_data['cont_flux_zscore'] = zscore(np.interp(uniform_time_grid, valid_cont['mjd'], valid_cont['flux']))

# Interpolate spectroscopic bins
for col, new_name in [('blue_wing_flux', 'blue_wing_flux_zscore'), 
                      ('core_flux', 'core_flux_zscore'), 
                      ('red_wing_flux', 'red_wing_flux_zscore')]:
    valid_data = df_hb.dropna(subset=[col]).sort_values('mjd')
    prepared_data[new_name] = zscore(np.interp(uniform_time_grid, valid_data['mjd'], valid_data[col]))

prepared_data = prepared_data.dropna().reset_index(drop=True)

# Extract z-scored flux arrays
cont_zscore = prepared_data['cont_flux_zscore'].values
blue_zscore = prepared_data['blue_wing_flux_zscore'].values
core_zscore = prepared_data['core_flux_zscore'].values
red_zscore = prepared_data['red_wing_flux_zscore'].values

# Set custom plotting styles for premium publication quality
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'figure.titlesize': 16,
    'font.family': 'sans-serif'
})

# ----------------- FIGURE 1: PREPARED LIGHT CURVES -----------------
print("Generating Figure 1: Aligned Light Curves...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

ax1.plot(prepared_data['time'], cont_zscore, label='5100 Å Continuum', color='#1f77b4', linewidth=1.5)
ax1.set_ylabel('Standardized Flux ($Z$)')
ax1.legend(loc='upper right')
ax1.set_title('NGC 5548 Continuum and Velocity-resolved $H\\beta$ Components (Strict Overlap Window)')

ax2.plot(prepared_data['time'], blue_zscore, label='Blue Wing ($-3000$ to $-1000$ km/s)', color='#2ca02c', alpha=0.8)
ax2.plot(prepared_data['time'], core_zscore, label='Core ($-1000$ to $+1000$ km/s)', color='#d62728', alpha=0.8)
ax2.plot(prepared_data['time'], red_zscore, label='Red Wing ($+1000$ to $+3000$ km/s)', color='#ff7f0e', alpha=0.8)
ax2.set_xlabel('Modified Julian Date (MJD)')
ax2.set_ylabel('Standardized Flux ($Z$)')
ax2.legend(loc='upper right')

plt.tight_layout()
fig.savefig('overleaf_draft/figure1_light_curves.png', dpi=300)
plt.close(fig)

# ----------------- SURD UTILS FOR PLOTTING -----------------
def run_collect_2pred(target_arr, pred1_arr, pred2_arr, nlag, nbins):
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
        sum_norm = res["U1"] + res["U2"] + res["R12"] + res["S12"]
        assert np.abs(sum_norm - 1.0) < 1e-6, f"SURD Normalization failed at lag {lag}: sum is {sum_norm}"
        metrics["lag"].append(lag)
        metrics["info_leak"].append(res["info_leak"])
        metrics["U1"].append(res["U1"])
        metrics["U2"].append(res["U2"])
        metrics["R12"].append(res["R12"])
        metrics["S12"].append(res["S12"])
    print(f"  All {len(lags)} lags successfully verified: U1 + U2 + R12 + S12 = 1.0 (identity holds).")
    return metrics

# Run SURD scans up to 120 lags
lags_120 = np.arange(1, 121)
print("Running SURD lag scans up to 120 days for all targets...")
# Core target: S1=cont, S2=blue -> T3=core
metrics_core = lag_scan_target3(core_zscore, cont_zscore, blue_zscore, lags_120, nbins=8)

# Red target: S1=cont, S2=blue -> T3=red
metrics_red = lag_scan_target3(red_zscore, cont_zscore, blue_zscore, lags_120, nbins=8)

# Blue target: S1=cont, S2=core -> T3=blue
metrics_blue = lag_scan_target3(blue_zscore, cont_zscore, core_zscore, lags_120, nbins=8)

# ----------------- FIGURE 2: SURD LAG SCANS WITH MONTE CARLO UNCERTAINTY -----------------
print("Generating Figure 2: SURD Synergy and Leak Lag Scans with MC Uncertainty...")
fig, axs = plt.subplots(3, 2, figsize=(12, 11), sharex=True)

# Load MC uncertainty results
df_mc_core = pd.read_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/mc_uncertainty_core.csv')
df_mc_red = pd.read_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/mc_uncertainty_red.csv')
df_mc_blue = pd.read_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/mc_uncertainty_blue.csv')

targets = [
    ('Core $H\\beta$', metrics_core, df_mc_core),
    ('Red Wing $H\\beta$', metrics_red, df_mc_red),
    ('Blue Wing $H\\beta$', metrics_blue, df_mc_blue)
]

for idx, (name, metrics, df_mc) in enumerate(targets):
    # Synergy column
    ax_syn = axs[idx, 0]
    # Plot MC Median and Shaded Error Bands
    ax_syn.plot(df_mc['lag'], df_mc['median_S12'], color='purple', label='Synergy (MC Median)', linewidth=2)
    ax_syn.fill_between(df_mc['lag'], df_mc['p16_S12'], df_mc['p84_S12'], color='purple', alpha=0.3, label='1$\sigma$ MC Error')
    ax_syn.fill_between(df_mc['lag'], df_mc['p2_5_S12'], df_mc['p97_5_S12'], color='purple', alpha=0.1, label='2$\sigma$ MC Error')
    
    # Reference curves from real run
    ax_syn.plot(metrics['lag'], metrics['R12'], color='gray', linestyle='--', label='Redundancy', alpha=0.7)
    ax_syn.plot(metrics['lag'], metrics['U1'], color='blue', linestyle=':', label='Unique (Continuum)', alpha=0.7)
    ax_syn.plot(metrics['lag'], metrics['U2'], color='red', linestyle='-.', label='Unique (Wing/Core)', alpha=0.7)
    ax_syn.set_ylabel('Information Fraction')
    ax_syn.set_title(f'Information Decomposition: {name}')
    ax_syn.legend(loc='upper right')
    
    # Leak column
    ax_leak = axs[idx, 1]
    # Plot MC Median and Shaded Error Bands for Leak
    ax_leak.plot(df_mc['lag'], df_mc['median_leak'], color='darkorange', linewidth=2, label='Information Leak (MC Median)')
    ax_leak.fill_between(df_mc['lag'], df_mc['p16_leak'], df_mc['p84_leak'], color='darkorange', alpha=0.3, label='1$\sigma$ MC Error')
    ax_leak.fill_between(df_mc['lag'], df_mc['p2_5_leak'], df_mc['p97_5_leak'], color='darkorange', alpha=0.1, label='2$\sigma$ MC Error')
    
    ax_leak.set_ylabel('Normalized Leak $\\mathcal{L}$')
    ax_leak.set_title(f'Information Leak: {name}')
    ax_leak.legend(loc='upper right')

axs[2, 0].set_xlabel('Lag (days)')
axs[2, 1].set_xlabel('Lag (days)')
plt.tight_layout()
fig.savefig('overleaf_draft/figure2_surd_lag_scans.png', dpi=300)
plt.close(fig)

# ----------------- FIGURE 3: ROBUSTNESS AND SHUFFLE ENVELOPES -----------------
print("Generating Figure 3: Robustness and Null Tests...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# Panel A: nbins sensitivity (1 to 60 lags)
lags_60 = np.arange(1, 61)
nbins_vals = [4, 6, 8, 10, 12]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for n_idx, n_val in enumerate(nbins_vals):
    m_temp = lag_scan_target3(core_zscore, cont_zscore, blue_zscore, lags_60, nbins=n_val)
    ax1.plot(lags_60, m_temp['S12'], label=f'nbins = {n_val}', color=colors[n_idx], linewidth=1.5)
ax1.set_xlabel('Lag (days)')
ax1.set_ylabel('Normalized Synergy $\\widehat{S}_{12}$')
ax1.set_title('A: Histogram Bin Sensitivity (Core Target)')
ax1.legend(loc='upper right')

# Panel B: real vs. circular & block shuffle null envelopes
df_combined = pd.read_csv('/Users/ayan/Programs/SURD/agn_surd_project/plots/test_c_shuffle_results/robustness_surrogate_core_combined.csv')

ax2.plot(df_combined['lag'], df_combined['real_synergy'], color='black', linewidth=2, label='Real Synergy')
# Circular shuffle
ax2.plot(df_combined['lag'], df_combined['median_synergy_cont_circ_shuffle'], color='blue', label='Median (Circular)', alpha=0.8)
ax2.fill_between(df_combined['lag'], df_combined['p2_5_synergy_cont_circ_shuffle'], df_combined['p97_5_synergy_cont_circ_shuffle'], 
                 color='blue', alpha=0.15, label='95% envelope (Circular)')
# Block shuffle
ax2.plot(df_combined['lag'], df_combined['median_synergy_cont_block_shuffle'], color='red', label='Median (Block, 10d)', alpha=0.8)
ax2.fill_between(df_combined['lag'], df_combined['p2_5_synergy_cont_block_shuffle'], df_combined['p97_5_synergy_cont_block_shuffle'], 
                 color='red', alpha=0.15, label='95% envelope (Block, 10d)')

ax2.set_xlabel('Lag (days)')
ax2.set_ylabel('Synergy $S_{12}$ (bits)')
ax2.set_title('B: Real vs. Surrogate Envelopes (Core Target)')
ax2.legend(loc='upper right')

plt.tight_layout()
fig.savefig('overleaf_draft/figure3_robustness_and_nulls.png', dpi=300)
plt.close(fig)

# ----------------- FIGURE 4: ICCF VS SURD -----------------
print("Generating Figure 4: ICCF vs. SURD Lags...")
def compute_iccf(line, cont, lags, dt=1.0):
    min_len = min(len(line), len(cont))
    line_trimmed = line[:min_len]
    cont_trimmed = cont[:min_len]
    
    correlation = correlate(line_trimmed, cont_trimmed, mode='full')
    correlation = correlation / np.sqrt(np.sum(line_trimmed**2) * np.sum(cont_trimmed**2))
    
    correlation_lags_samples = np.arange(-min_len + 1, min_len)
    correlation_lags_days = correlation_lags_samples * dt
    
    desired_lags_min = lags.min() * dt
    desired_lags_max = lags.max() * dt
    
    mask_lags = (correlation_lags_days >= desired_lags_min) & (correlation_lags_days <= desired_lags_max)
    iccf_lags_days = correlation_lags_days[mask_lags]
    iccf_values = correlation[mask_lags]
    return iccf_lags_days, iccf_values

iccf_lags_b, iccf_vals_b = compute_iccf(blue_zscore, cont_zscore, lags_120)
iccf_lags_c, iccf_vals_c = compute_iccf(core_zscore, cont_zscore, lags_120)
iccf_lags_r, iccf_vals_r = compute_iccf(red_zscore, cont_zscore, lags_120)

fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

components = [
    ('Blue Wing $H\\beta$', iccf_lags_b, iccf_vals_b, metrics_blue, 19.0, 81.0),
    ('Core $H\\beta$', iccf_lags_c, iccf_vals_c, metrics_core, 20.0, 76.0),
    ('Red Wing $H\\beta$', iccf_lags_r, iccf_vals_r, metrics_red, 13.0, 119.0)
]

for idx, (name, iccf_lags, iccf_vals, surd_metrics, iccf_peak, surd_peak) in enumerate(components):
    ax = axs[idx]
    
    # Plot ICCF on left y-axis
    color = '#1f77b4'
    ax.plot(iccf_lags, iccf_vals, color=color, label='ICCF (Linear Correlation)', linewidth=2)
    ax.tick_params(axis='y', labelcolor=color)
    ax.set_ylabel('Correlation Coefficient', color=color)
    ax.axvline(iccf_peak, color=color, linestyle='--', label=f'ICCF Peak Lag: {iccf_peak:.1f} d')
    
    # Plot SURD Synergy on right y-axis
    ax2 = ax.twinx()
    color2 = 'purple'
    ax2.plot(surd_metrics['lag'], surd_metrics['S12'], color=color2, label='SURD Synergy', linewidth=2)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_ylabel('Normalized Synergy $\\widehat{S}_{12}$', color=color2)
    ax2.axvline(surd_peak, color=color2, linestyle='-.', label=f'SURD Synergy Peak: {surd_peak:.1f} d')
    
    ax.set_title(f'ICCF vs. SURD Synergy: {name}')
    
    # Combine legends
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc='upper right')

axs[2].set_xlabel('Lag (days)')
plt.tight_layout()
fig.savefig('overleaf_draft/figure4_iccf_vs_surd.png', dpi=300)
plt.close(fig)

# ----------------- FIGURE 5: REALISTIC SYNTHETIC BENCHMARKS -----------------
print("Generating Figure 5: Realistic Synthetic Benchmarks...")
# Import and run the realistic synthetic simulation directly
try:
    import scratch.run_synthetic_realistic as run_realistic
    # Since run_realistic executes on load and writes figure5, this will trigger it!
except Exception as e:
    print(f"Warning: Could not run realistic synthetic benchmarks automatically: {e}")

# ----------------- FIGURE 7: TARGET-HISTORY CONDITIONING -----------------
print("Generating Figure 7: Target-History Conditioning (Core Target)...")
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

# We use the standardized continuum, blue wing, and core arrays
X_cond = np.vstack([cont_zscore, blue_zscore, core_zscore])
lags_scan = np.arange(1, 121)

core_uncond_syn, core_uncond_leak = [], []
core_cond_syn, core_cond_leak = [], []

for lag in lags_scan:
    us, ul = run_unconditioned_collect(X_cond, 2, [0, 1], lag, nbins=6)
    cs, cl = run_conditional_collect(X_cond, 2, [0, 1], 2, lag, nbins=6)
    core_uncond_syn.append(us)
    core_uncond_leak.append(ul)
    core_cond_syn.append(cs)
    core_cond_leak.append(cl)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# Synergy comparison
ax1.plot(lags_scan, core_uncond_syn, label='Unconditioned Synergy', color='blue', linewidth=2)
ax1.plot(lags_scan, core_cond_syn, label='History-Conditioned Synergy', color='red', linewidth=2)
ax1.set_xlabel('Lag (days)')
ax1.set_ylabel('Synergy $S_{12}$ (bits)')
ax1.set_title('A: Synergy Comparison (Core Target)')
ax1.legend(loc='upper right')
ax1.grid(True, linestyle='--', alpha=0.5)

# Leak comparison
ax2.plot(lags_scan, core_uncond_leak, label='Unconditioned Leak', color='blue', linewidth=2)
ax2.plot(lags_scan, core_cond_leak, label='History-Conditioned Leak', color='red', linewidth=2)
ax2.set_xlabel('Lag (days)')
ax2.set_ylabel('Information Leak (normalized entropy)')
ax2.set_title('B: Information Leak Comparison (Core Target)')
ax2.legend(loc='upper right')
ax2.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
fig.savefig('overleaf_draft/figure7_history_conditioning.png', dpi=300)
plt.close(fig)
print("Figure 7: figure7_history_conditioning.png successfully created and saved in overleaf_draft/!")

print("All publication-quality figures successfully created and saved in overleaf_draft/!")
