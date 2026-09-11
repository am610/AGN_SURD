import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROCESSED = ROOT / 'agn_surd_project' / 'processed'
OVERLEAF = ROOT / 'overleaf_draft'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/surd-matplotlib')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore

# Configure path to SURD utilities
sys.path.append(str(ROOT / 'SURD'))
sys.path.append(str(ROOT / 'SURD' / 'utils'))
from utils import surd

# ----------------- 1. LOAD AND PREPARE DATA -----------------
print("Loading and aligning dataset...")
cont_path = ROOT / 'agn_surd_project' / 'agn_data' / 'ngc5548_agnwatch' / 'c5100.dat'
hb_bins_path = PROCESSED / 'ngc5548_hb_velocity_bins.csv'

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


def lag_scan_3pred(target_arr, predictors, lags, nbins=8):
    metrics = {"lag": [], "normalized_synergy": [], "synergy_bits": [], "joint_mi": []}
    for lag in lags:
        data = np.column_stack([target_arr[lag:], *(predictor[:-lag] for predictor in predictors)])
        hist, _ = np.histogramdd(data, bins=nbins)
        hist /= hist.sum()
        _, i_s, mi, _ = surd.surd(hist)
        joint_mi = mi.get((1, 2, 3), np.nan)
        synergy = sum(i_s.get(key, 0.0) for key in ((1, 2), (1, 3), (2, 3), (1, 2, 3)))
        metrics["lag"].append(lag)
        metrics["normalized_synergy"].append(synergy / joint_mi if joint_mi > 0 else np.nan)
        metrics["synergy_bits"].append(synergy)
        metrics["joint_mi"].append(joint_mi)
    return metrics

# ----------------- FIGURE 2: ROUND-ROBIN 3-PREDICTOR SURD LAG SCANS -----------------
print("Generating Figure 2: Round-Robin 3-Predictor SURD Synergy and Leak scans...")
lags_200 = np.arange(1, 201)
metrics_core = lag_scan_target3(core_zscore, cont_zscore, blue_zscore, lags_200, nbins=8)
metrics_red = lag_scan_target3(red_zscore, cont_zscore, blue_zscore, lags_200, nbins=8)
metrics_blue = lag_scan_target3(blue_zscore, cont_zscore, core_zscore, lags_200, nbins=8)
df_rr = pd.read_csv(PROCESSED / 'round_robin_results.csv')

fig, axs = plt.subplots(4, 3, figsize=(17, 14), sharex=True)

# List of targets in order: Continuum, Red Wing, Core, Blue Wing
target_list = ['continuum', 'red_wing', 'core', 'blue_wing']
display_names = ['Continuum [reverse]', 'Red Wing H$\\beta$', 'Core H$\\beta$', 'Blue Wing H$\\beta$']

# Color and linestyle maps for Unique terms
unique_styles = {
    'u_continuum': ('blue', ':', 'U_{\\mathrm{cont}}'),
    'u_red_wing': ('green', '--', 'U_{\\mathrm{red}}'),
    'u_core': ('red', '-.', 'U_{\\mathrm{core}}'),
    'u_blue_wing': ('purple', '-.', 'U_{\\mathrm{blue}}')
}

for idx, target in enumerate(target_list):
    df_t = df_rr[df_rr['target'] == target]
    
    # Left column: Information Decomposition
    ax_syn = axs[idx, 0]
    ax_syn.plot(df_t['lag'], df_t['syn'], color='purple', label='Total Synergy $\\widehat{S}$', linewidth=2.5)
    ax_syn.plot(df_t['lag'], df_t['red'], color='gray', label='Total Redundancy $\\widehat{R}$', linewidth=2.5)
    
    # Plot the three uniques
    for col in df_t.columns:
        if col.startswith('u_') and col in unique_styles:
            # Check if this column is not entirely NaN for this target
            if df_t[col].notna().any():
                color, linestyle, label_name = unique_styles[col]
                ax_syn.plot(df_t['lag'], df_t[col], color=color, linestyle=linestyle, label=f'${label_name}$', alpha=0.85, linewidth=1.8)
            
    ax_syn.set_ylabel('Information Fraction')
    ax_syn.set_title(f'Normalized Decomposition: {display_names[idx]}')
    ax_syn.legend(loc='upper right', fontsize=9)
    ax_syn.grid(True, linestyle='--', alpha=0.5)
    
    # Middle column: absolute information and the normalization denominator
    ax_abs = axs[idx, 1]
    ax_abs.plot(df_t['lag'], df_t['joint_mi'], color='#15616d', label='Joint MI (denominator)', linewidth=2.3)
    ax_abs.plot(df_t['lag'], df_t['syn_bits'], color='#78290f', label='Total synergy (bits)', linewidth=2.0)
    ax_abs.set_ylabel('Information (bits)')
    ax_abs.set_title(f'Absolute Information: {display_names[idx]}')
    ax_abs.legend(loc='upper right', fontsize=9)
    ax_abs.grid(True, linestyle='--', alpha=0.5)

    # Right column: Information Leak
    ax_leak = axs[idx, 2]
    ax_leak.plot(df_t['lag'], df_t['leak'], color='darkorange', label='Information Leak $\\mathcal{L}$', linewidth=2.5)
    ax_leak.set_ylabel('Normalized Leak $\\mathcal{L}$')
    ax_leak.set_title(f'Information Leak: {display_names[idx]}')
    ax_leak.legend(loc='upper right', fontsize=9)
    ax_leak.grid(True, linestyle='--', alpha=0.5)
    ax_leak.set_ylim([0, 1])

axs[3, 0].set_xlabel('Lag (days)')
axs[3, 1].set_xlabel('Lag (days)')
axs[3, 2].set_xlabel('Lag (days)')
plt.tight_layout()
fig.savefig('overleaf_draft/figure2_surd_lag_scans.png', dpi=300)
plt.close(fig)

# ----------------- FIGURE 3: ROBUSTNESS AND SHUFFLE ENVELOPES -----------------
print("Generating Figure 3: Robustness and Null Tests...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# Panel A: matched three-predictor bin sensitivity (1 to 60 lags)
lags_60 = np.arange(1, 61)
nbins_vals = [4, 6, 8, 10, 12]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for n_idx, n_val in enumerate(nbins_vals):
    m_temp = lag_scan_3pred(core_zscore, [cont_zscore, red_zscore, blue_zscore], lags_60, nbins=n_val)
    ax1.plot(lags_60, m_temp['normalized_synergy'], label=f'nbins = {n_val}', color=colors[n_idx], linewidth=1.5)
ax1.set_xlabel('Lag (days)')
ax1.set_ylabel('Normalized Total Synergy $\\widehat{S}$')
ax1.set_title('A: Matched 3-predictor Bin Sensitivity (Core Target)')
ax1.legend(loc='upper right')

# Panel B: exact matched circular-shift null envelope
df_matched = pd.read_csv(PROCESSED / 'round_robin_empirical_null_curves.csv')
df_matched = df_matched[df_matched['target'] == 'core']
ax2.plot(df_matched['lag'], df_matched['real_synergy'], color='black', linewidth=2, label='Observed')
ax2.plot(df_matched['lag'], df_matched['null_median'], color='blue', label='Null median', alpha=0.8)
ax2.fill_between(df_matched['lag'], df_matched['null_p2_5'], df_matched['null_p97_5'],
                 color='blue', alpha=0.18, label='95% pointwise envelope')

ax2.set_xlabel('Lag (days)')
ax2.set_ylabel('Normalized Total Synergy $\\widehat{S}$')
ax2.set_title('B: Matched 3-predictor Null, 999 Shifts')
ax2.legend(loc='upper right')

plt.tight_layout()
fig.savefig('overleaf_draft/figure3_robustness_and_nulls.png', dpi=300)
plt.close(fig)

# ----------------- FIGURE 4: ICCF VS SURD -----------------
print("Generating Figure 4: ICCF vs. SURD Lags...")
df_iccf = pd.read_csv(PROCESSED / 'iccf_curves.csv')
df_iccf_summary = pd.read_csv(PROCESSED / 'iccf_summary.csv').set_index('component')

fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

df_rr = pd.read_csv(PROCESSED / 'round_robin_results.csv')
rr_blue = df_rr[df_rr['target'] == 'blue_wing']
rr_core = df_rr[df_rr['target'] == 'core']
rr_red = df_rr[df_rr['target'] == 'red_wing']

components = [
    ('Blue Wing $H\\beta$', 'blue_wing', rr_blue),
    ('Core $H\\beta$', 'core', rr_core),
    ('Red Wing $H\\beta$', 'red_wing', rr_red)
]

for idx, (name, component, surd_df) in enumerate(components):
    ax = axs[idx]
    iccf_component = df_iccf[df_iccf['component'] == component]
    summary = df_iccf_summary.loc[component]
    surd_peak_row = surd_df.loc[surd_df['syn'].idxmax()]
    surd_peak = float(surd_peak_row['lag'])
    
    # Plot ICCF on left y-axis
    color = '#1f77b4'
    ax.plot(iccf_component['lag'], iccf_component['r'], color=color, label='Bidirectional ICCF', linewidth=2)
    ax.tick_params(axis='y', labelcolor=color)
    ax.set_ylabel('Correlation Coefficient', color=color)
    ax.axvspan(summary['centroid_p16'], summary['centroid_p84'], color=color, alpha=0.15)
    ax.axvline(summary['centroid_median'], color=color, linestyle='--',
               label=f"ICCF centroid: {summary['centroid_median']:.1f} d")
    
    # Plot SURD Synergy on right y-axis
    ax2 = ax.twinx()
    color2 = 'purple'
    ax2.plot(surd_df['lag'], surd_df['syn'], color=color2, label='SURD Synergy', linewidth=2)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_ylabel('Normalized Synergy $\\widehat{S}$', color=color2)
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
from run_final_validations import plot_synthetic

synthetic_curves = pd.read_csv(
    'agn_surd_project/processed/synthetic_validation_curves.csv'
)
plot_synthetic(synthetic_curves)

# ----------------- FIGURE 7: TARGET-HISTORY CONDITIONING -----------------
print("Generating Figure 8: Target-History Conditioning (Core Target)...")
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

# Load conditional surrogates (runs up to 120 days)
df_cond_null = pd.read_csv(PROCESSED / 'conditional_surrogate_results.csv')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# Synergy comparison with conditional surrogate envelope
ax1.plot(lags_scan, core_uncond_syn, label='Unconditioned Synergy', color='blue', linewidth=2)
ax1.plot(lags_scan, core_cond_syn, label='History-Conditioned Synergy', color='red', linewidth=2)
# Overlay surrogate envelope
ax1.plot(df_cond_null['lag'], df_cond_null['median_null_cond'], color='orange', label='Median (Shift null)', linestyle='--', alpha=0.8)
ax1.fill_between(df_cond_null['lag'], df_cond_null['p2_5_null_cond'], df_cond_null['p97_5_null_cond'], 
                 color='orange', alpha=0.15, label='95% cond envelope')

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
fig.savefig('overleaf_draft/figure8_history_conditioning.png', dpi=300)
plt.close(fig)
print("Figure 8: figure8_history_conditioning.png successfully created and saved in overleaf_draft/!")

# ----------------- FIGURE 7: MATCHED NEGATIVE-CONTROL ABLATION -----------------
print("Generating Figure 7: Matched Negative-Control Ablation...")
df_ablation = pd.read_csv(PROCESSED / 'matched_negative_control_curves.csv')

fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharex=True)
styles = {
    'complete_daily': ('#2a9d8f', 'Complete daily'),
    'observed_only': ('#e9c46a', 'Observed epochs only'),
    'observed_plus_interpolation': ('#e76f51', 'Observed + interpolation'),
}
for condition, (color, label) in styles.items():
    subset = df_ablation[df_ablation['condition'] == condition]
    axes[0].plot(subset['lag'], subset['normalized_median'], color=color, label=label, linewidth=2)
    axes[1].plot(subset['lag'], subset['synergy_bits_median'], color=color, label=label, linewidth=2)
    axes[2].plot(subset['lag'], subset['joint_mi_median'], color=color, label=label, linewidth=2)
axes[0].set_ylabel('Median normalized total synergy')
axes[1].set_ylabel('Median total synergy (bits)')
axes[2].set_ylabel('Median joint MI (bits)')
for axis, title in zip(axes, ('A: Normalized statistic', 'B: Absolute synergy', 'C: Normalization denominator')):
    axis.set_xlabel('Lag (days)')
    axis.set_title(title)
    axis.grid(True, linestyle='--', alpha=0.5)
axes[0].legend(fontsize=9)

plt.tight_layout()
fig.savefig('overleaf_draft/figure7_seasonal_aliasing.png', dpi=300)
plt.close(fig)
print("Figure 7: figure7_seasonal_aliasing.png successfully created and saved in overleaf_draft/!")

# ----------------- FIGURE 9: CONDITIONAL BINNING SENSITIVITY -----------------
print("Generating Figure 9: Conditional Binning Sensitivity...")
from run_final_validations import plot_binning

binning_curves = pd.read_csv(
    'agn_surd_project/processed/conditional_binning_sensitivity_curves.csv'
)
plot_binning(binning_curves)
print("Figure 9: figure9_conditional_binning_sensitivity.png successfully created and saved in overleaf_draft/!")

print("All publication-quality figures successfully created and saved in overleaf_draft/!")
