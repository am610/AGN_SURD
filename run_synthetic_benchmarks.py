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

# ----------------- 1. DRW SIMULATOR -----------------
def generate_drw(n_days, tau_param=20.0, sigma_param=0.2, dt=1.0, seed=None):
    if seed is not None:
        np.random.seed(seed)
    time = np.arange(0, n_days, dt)
    n = len(time)
    flux = np.zeros(n)
    
    # Ornstein-Uhlenbeck / Damped Random Walk parameters
    mu = 0.0
    flux[0] = np.random.normal(mu, sigma_param)
    
    for i in range(1, n):
        # DRW recurrence relation
        coeff = np.exp(-dt / tau_param)
        noise_std = sigma_param * np.sqrt(1 - coeff**2)
        flux[i] = mu + coeff * (flux[i-1] - mu) + np.random.normal(0, noise_std)
        
    return flux

# ----------------- 2. SURD CALCULATION HELPERS -----------------
def run_collect(X, nvars, nlag, nbins):
    results = {}
    for i in range(nvars):
        Y = np.vstack([X[i, nlag:], X[:, :-nlag]])
        hist, _ = np.histogramdd(Y.T, nbins)
        I_R, I_S, MI, info_leak = surd.surd(hist)
        results[i] = {"I_R": I_R, "I_S": I_S, "MI": MI, "info_leak": info_leak}
    return results

def lag_scan_target3(X, lags, nbins=8):
    metrics = {"lag": [], "info_leak": [], "MI1": [], "MI2": [], "U1": [], "U2": [], "R12": [], "S12": []}
    for lag in lags:
        res = run_collect(X=X, nvars=3, nlag=lag, nbins=nbins)[2]
        metrics["lag"].append(lag)
        metrics["info_leak"].append(res["info_leak"])
        metrics["MI1"].append(res["MI"].get((1,), np.nan))
        metrics["MI2"].append(res["MI"].get((2,), np.nan))
        metrics["U1"].append(res["I_R"].get((1,), np.nan))
        metrics["U2"].append(res["I_R"].get((2,), np.nan))
        metrics["R12"].append(res["I_R"].get((1, 2), np.nan))
        metrics["S12"].append(res["I_S"].get((1, 2), np.nan))
    return metrics

# ----------------- 3. RUN SYNTHETIC BENCHMARKS -----------------
n_days = 2000
dt = 1.0
lags = np.arange(1, 41)
nbins = 6

print("Generating synthetic datasets and running SURD lag scans...")

# --- Case 1: Single Driver ---
# S1 = DRW driver, S2 = Independent White Noise, Target T3 = S1 delayed by 15 days + noise
s1 = generate_drw(n_days, tau_param=50.0, sigma_param=0.3, dt=dt, seed=10)
s2 = np.random.normal(0, 0.3, len(s1))
t3 = np.roll(s1, 15) + np.random.normal(0, 0.05, len(s1))
# Trim edges after roll
trim = 50
X1 = np.vstack([zscore(s1[trim:]), zscore(s2[trim:]), zscore(t3[trim:])])
metrics_case1 = lag_scan_target3(X1, lags, nbins=nbins)

# --- Case 2: Redundant Proxies ---
# S1 = DRW driver, S2 = S1 + White Noise (Proxy), Target T3 = S1 delayed by 15 days + noise
s1 = generate_drw(n_days, tau_param=50.0, sigma_param=0.3, dt=dt, seed=20)
s2 = s1 + np.random.normal(0, 0.1, len(s1))
t3 = np.roll(s1, 15) + np.random.normal(0, 0.05, len(s1))
X2 = np.vstack([zscore(s1[trim:]), zscore(s2[trim:]), zscore(t3[trim:])])
metrics_case2 = lag_scan_target3(X2, lags, nbins=nbins)

# --- Case 3: Synergistic Drivers ---
# S1, S2 = Independent DRW drivers. Target T3 = S1(t-15) * S2(t-15) (nonlinear product) + noise
s1 = generate_drw(n_days, tau_param=50.0, sigma_param=0.3, dt=dt, seed=30)
s2 = generate_drw(n_days, tau_param=30.0, sigma_param=0.3, dt=dt, seed=40)
# Make them strictly positive to avoid canceling signs in product
s1_p = s1 - s1.min() + 0.1
s2_p = s2 - s2.min() + 0.1
t3 = np.roll(s1_p, 15) * np.roll(s2_p, 15) + np.random.normal(0, 0.01, len(s1))
X3 = np.vstack([zscore(s1_p[trim:]), zscore(s2_p[trim:]), zscore(t3[trim:])])
metrics_case3 = lag_scan_target3(X3, lags, nbins=nbins)

# --- Case 4: Two-Zone BLR Response ---
# S1 = Prompt continuum. S2 = Delayed reprocessed continuum (S1 delayed by 10 days).
# Target T3 = responses to both zones: 0.5 * S1(t-10) + 0.5 * S2(t-20)
s1 = generate_drw(n_days, tau_param=50.0, sigma_param=0.3, dt=dt, seed=50)
s2 = np.roll(s1, 10) + np.random.normal(0, 0.05, len(s1))
t3 = 0.5 * np.roll(s1, 10) + 0.5 * np.roll(s2, 20) + np.random.normal(0, 0.05, len(s1))
X4 = np.vstack([zscore(s1[trim:]), zscore(s2[trim:]), zscore(t3[trim:])])
metrics_case4 = lag_scan_target3(X4, lags, nbins=nbins)

# ----------------- 4. PLOT SUMMARY FIGURE -----------------
print("Generating Figure 5: Synthetic Benchmark Scans...")
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

cases = [
    ('Case 1: Single Driver (True Lag = 15d)', metrics_case1, axs[0, 0]),
    ('Case 2: Redundant Proxies (True Lag = 15d)', metrics_case2, axs[0, 1]),
    ('Case 3: Synergistic Drivers (True Lag = 15d)', metrics_case3, axs[1, 0]),
    ('Case 4: Two-Zone Response (Lags = 10d, 20d)', metrics_case4, axs[1, 1])
]

for title, m, ax in cases:
    ax.plot(m['lag'], m['S12'], color='purple', label='Synergy', linewidth=2)
    ax.plot(m['lag'], m['R12'], color='gray', linestyle='--', label='Redundancy', alpha=0.7)
    ax.plot(m['lag'], m['U1'], color='blue', linestyle=':', label='Unique (S1)', alpha=0.7)
    ax.plot(m['lag'], m['U2'], color='red', linestyle='-.', label='Unique (S2)', alpha=0.7)
    ax.set_xlabel('Lag (days)')
    ax.set_ylabel('Information (bits)')
    ax.set_title(title)
    ax.legend(loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
fig.savefig('overleaf_draft/figure5_synthetic_benchmarks.png', dpi=300)
plt.close(fig)

print("Figure 5: figure5_synthetic_benchmarks.png successfully created and saved in overleaf_draft/!")

# Save benchmark results as CSV for verification
output_dir = "/Users/ayan/Programs/SURD/agn_surd_project/processed/"
os.makedirs(output_dir, exist_ok=True)
for name, m in [("case1", metrics_case1), ("case2", metrics_case2), ("case3", metrics_case3), ("case4", metrics_case4)]:
    pd.DataFrame(m).to_csv(os.path.join(output_dir, f"synthetic_{name}_metrics.csv"), index=False)
print("Saved all synthetic benchmark metrics CSV files in processed/ directory.")
