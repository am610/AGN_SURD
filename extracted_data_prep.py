--- Cell 6 ---
import sys
import os

# Configure path to SURD repo utilities
# Change this if your SURD clone is somewhere else.
SURD_UTILS_PATH = '/content/agn_surd_project' # Corrected path to agn_surd_project

if not os.path.exists(SURD_UTILS_PATH):
    raise FileNotFoundError(
        f"Could not find SURD utils directory at {SURD_UTILS_PATH}. "
        "Update SURD_UTILS_PATH to match your local or Colab clone."
    )

sys.path.append(str(SURD_UTILS_PATH))

# The functions 'prepare_agn_timeseries', 'lag_scan_target3', 'plot_lag_metrics',
# 'resample_to_uniform_grid', 'zscore', and 'run_collect' are defined in a previous cell
# and are therefore globally available without needing explicit import from agn_helpers.

# Removed explicit imports for 'surd', 'it_tools', 'analytic_eqs' as they are not found in agn_surd_project.

print("Core SURD modules (surd, it_tools, analytic_eqs) and helper functions (defined in notebook) are ready.")
--- Cell 12 ---
import numpy as np

# Helper functions from cell '0cxiu-N1VNSO' to ensure they are defined in scope
def zscore(x: np.ndarray) -> np.ndarray:
    """
    Standardise an array to zero mean and unit variance.
    """
    x = np.asarray(x, dtype=float)
    return (x - np.mean(x)) / np.std(x)

def resample_to_uniform_grid(times, values, tmin=None, tmax=None, dt=1.0):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)

    if tmin is None:
        tmin = np.min(times)
    if tmax is None:
        tmax = np.max(times)

    grid = np.arange(tmin, tmax + dt, dt)
    vals = np.interp(grid, times, values)
    return grid, vals

def prepare_agn_timeseries(time_opt, flux_opt, time_xray, flux_xray, time_hb, flux_hb, dt=1.0):
    tmin = max(np.min(time_opt), np.min(time_xray), np.min(time_hb))
    tmax = min(np.max(time_opt), np.max(time_xray), np.max(time_hb))

    grid, opt = resample_to_uniform_grid(time_opt, flux_opt, tmin=tmin, tmax=tmax, dt=dt)
    _, xray = resample_to_uniform_grid(time_xray, flux_xray, tmin=tmin, tmax=tmax, dt=dt)
    _, hb = resample_to_uniform_grid(time_hb, flux_hb, tmin=tmin, tmax=tmax, dt=dt)

    opt = zscore(opt)
    xray = zscore(xray)
    hb = zscore(hb)

    X = np.vstack([opt, xray, hb])
    return grid, X


# Extract time and flux data from DataFrames
time_cont = df_ngc5548_cont['mjd'].values
flux_cont = df_ngc5548_cont['flux'].values

time_hb = df_ngc5548_hb['mjd'].values
flux_hb = df_ngc5548_hb['flux'].values

# Define the uniform grid spacing (in days, based on MJD units)
dt_real = 1.0 # 1-day cadence

print(f"Preparing real AGN time series with a uniform cadence of {dt_real} day(s).")

# Prepare the three aligned standardised bands using the helper function:
# Signal 1 = continuum (optical)
# Signal 2 = continuum (placeholder X-ray)
# Signal 3 = Hβ
# For prepare_agn_timeseries, we need three distinct time and flux inputs.
# For now, we will duplicate the continuum data for the 'xray' input.
# This is a temporary setup until actual multi-band/X-ray data is incorporated.

# Note: The prepare_agn_timeseries function expects time_opt, flux_opt, time_xray, flux_xray, time_hb, flux_hb
# We'll feed the continuum data twice for the first two signals.

real_grid, X_real = prepare_agn_timeseries(
    time_cont, flux_cont,       # Signal 1: Continuum
    time_cont, flux_cont,       # Signal 2: Continuum (placeholder for X-ray)
    time_hb, flux_hb,           # Signal 3: Hβ
    dt=dt_real
)

print("\n--- Prepared Real AGN Data Summary ---")
print(f"Shape of X_real (n_signals, n_times): {X_real.shape}")
print(f"Number of signals: {X_real.shape[0]}")
print(f"Number of time samples: {X_real.shape[1]}")
print(f"Cadence (dt): {dt_real} days")
print(f"Total duration of aligned data: {real_grid[-1] - real_grid[0]:.2f} days")
print(f"Time range: MJD {real_grid[0]:.2f} to {real_grid[-1]:.2f}")

# Simple check for interpolation statistics (e.g., density of original points vs. grid)
original_cont_points = len(time_cont)
original_hb_points = len(time_hb)
interpolated_points = X_real.shape[1]

print(f"Original continuum points: {original_cont_points}")
print(f"Original Hβ points: {original_hb_points}")
print(f"Interpolated grid points: {interpolated_points}")

print("\nInterpolation statistics and shapes printed. Proceeding to diagnostic plots.")
--- Cell 15 ---
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from typing import Dict, Any

# Ensure SURD repository is cloned
!git clone https://github.com/Computational-Turbulence-Group/SURD.git /content/SURD

# Configure path to SURD repo utilities
# Change this if your SURD clone is somewhere else.
SURD_UTILS_PATH = Path("/content/SURD/utils") # Corrected path to SURD utils directory

if not SURD_UTILS_PATH.exists():
    raise FileNotFoundError(
        f"Could not find SURD utils directory at {SURD_UTILS_PATH}. "
        "Update SURD_UTILS_PATH to match your local or Colab clone."
    )

sys.path.append(str(SURD_UTILS_PATH))

# Install pymp-pypi dependency for surd
!pip install pymp-pypi

# Core SURD imports
import surd

def zscore(x: np.ndarray) -> np.ndarray:
    """
    Standardise an array to zero mean and unit variance.
    Handles cases where standard deviation is zero and NaN values.
    """
    x = np.asarray(x, dtype=float)
    # Use nanmean and nanstd to ignore NaNs when calculating statistics
    mean_val = np.nanmean(x)
    std_dev = np.nanstd(x)

    if std_dev == 0 or np.isnan(std_dev):
        # If std_dev is zero or NaN (e.g., if x was all NaNs or constant non-NaNs),
        # return an array of zeros.
        z_scored = np.zeros_like(x)
    else:
        z_scored = (x - mean_val) / std_dev

    # Replace any remaining NaNs in the output (corresponding to NaNs in original x) with 0.0
    z_scored[np.isnan(z_scored)] = 0.0
    return z_scored

def run_collect(
    X: np.ndarray,
    nvars: int,
    nlag: int,
    nbins: int,
    axs: np.ndarray | None = None,
    print_results: bool = False,
) -> Dict[int, Dict[str, Any]]:
    """
    Run SURD for each target signal and return a clean dictionary of results.

    Why this helper exists
    ----------------------
    The repo's `surd.run(...)` function loops over all targets and plots all of
    them, but returns only the decomposition of the final target in the loop.
    For analysis, it is much more useful to collect results for every target.

    Parameters
    ----------
    X : np.ndarray
        Array of shape (nvars, ntimes). Each row is one signal.
    nvars : int
        Number of signals.
    nlag : int
        Lag in samples.
    nbins : int
        Number of histogram bins used by np.histogramdd.
    axs : np.ndarray or None
        Optional axes array of shape (nvars, 2) for SURD plotting.
    print_results : bool
        If True, print the SURD decomposition for each target.

    Returns
    -------
    dict
        Dictionary keyed by target index, with values containing:
            - "I_R"
            - "I_S"
            - "MI"
            - "info_leak"
    """
    results = {}

    for i in range(nvars):
        # Future of target i versus present/past of all variables
        Y = np.vstack([X[i, nlag:], X[:, :-nlag]])

        # Multi-dimensional histogram approximation to the joint distribution
        hist, _ = np.histogramdd(Y.T, nbins)

        # SURD decomposition
        I_R, I_S, MI, info_leak = surd.surd(hist)

        results[i] = {
            "I_R": I_R,
            "I_S": I_S,
            "MI": MI,
            "info_leak": info_leak,
        }

        if print_results:
            print(f"\nSURD CAUSALITY FOR TARGET SIGNAL {i + 1}")
            surd.nice_print(I_R, I_S, MI, info_leak)

        if axs is not None:
            surd.plot(I_R, I_S, info_leak, axs[i, :], nvars, threshold=-0.01)

    return results

def lag_scan_target3(X, lags, nbins=8):
    metrics = {
        "lag": [],
        "info_leak": [],
        "MI1": [],
        "MI2": [],
        "U1": [],
        "U2": [],
        "R12": [],
        "S12": [],
    }

    for lag in lags:
        results = run_collect(X=X, nvars=3, nlag=lag, nbins=nbins)
        res = results[2]   # target index 2 = Signal 3 (Hβ)

        metrics["lag"].append(lag)
        metrics["info_leak"].append(res["info_leak"])
        metrics["MI1"].append(res["MI"].get((1,), np.nan))
        metrics["MI2"].append(res["MI"].get((2,), np.nan))
        metrics["U1"].append(res["I_R"].get((1,), np.nan))
        metrics["U2"].append(res["I_R"].get((2,), np.nan))
        metrics["R12"].append(res["I_R"].get((1, 2), np.nan))
        metrics["S12"].append(res["I_S"].get((1, 2), np.nan))

    return metrics

def plot_lag_metrics(metrics, true_tau=None, title="Lag scan for target Signal 3", show_plot=True, save_path=None):
    lags = np.array(metrics["lag"])

    plt.figure(figsize=(10, 6))
    plt.plot(lags, metrics["MI1"], marker="o", label="MI(1)")
    plt.plot(lags, metrics["MI2"], marker="o", label="MI(2)")
    plt.plot(lags, metrics["R12"], marker="o", label="Redundancy(1,2)")
    plt.plot(lags, metrics["S12"], marker="o", label="Synergy(1,2)")
    if true_tau is not None:
        plt.axvline(true_tau, linestyle="--", label="True lag")
    plt.xlabel("Lag")
    plt.ylabel("Information")
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    if save_path: plt.savefig(save_path + "_info.png")
    if show_plot: plt.show()
    plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(lags, metrics["info_leak"], marker="o")
    if true_tau is not None:
        plt.axvline(true_tau, linestyle="--")
    plt.xlabel("Lag")
    plt.ylabel("Information leak")
    plt.title(f"{title}: leak")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    if save_path: plt.savefig(save_path + "_leak.png")
    if show_plot: plt.show()
    plt.close()

# Define lag range in samples (which directly correspond to days since dt_real = 1.0)
lags_real_scan = np.arange(1, 61) # Lags from 1 to 60 days

print(f"Running SURD lag scan for NGC 5548 from {lags_real_scan.min()} to {lags_real_scan.max()} days.")

# Run SURD lag scan on the Hβ target (index 2)
metrics_ngc5548 = lag_scan_target3(X_real, lags_real_scan, nbins=8)

# Plot the results
plot_lag_metrics(
    metrics_ngc5548,
    true_tau=None, # No true lag for real data, will be estimated
    title='SURD Lag Scan: Cont (S1), Blue Wing (S2) -> Core Hβ (Target S3)'
)

print("SURD lag scan completed and plots generated for NGC 5548.")
--- Cell 43 ---
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Ensure the necessary functions are defined (zscore, resample_to_uniform_grid, prepare_agn_timeseries)
# These should be available from cell '0cxiu-N1VNSO'

# Prepare the dataframes for alignment
# df_matching_cont should now be defined from the previous cell if c5100.dat was found and matched

# Check if df_matching_cont exists and is not empty
if 'df_matching_cont' not in locals() or df_matching_cont.empty:
    raise ValueError("df_matching_cont is not defined or is empty. Cannot proceed with data preparation.")

df_cont_aligned = df_matching_cont.rename(columns={'mjd': 'time', 'flux': 'cont_flux', 'error': 'cont_error'})[['time', 'cont_flux', 'cont_error']]

# Rename columns in df_hb_velocity for clarity and merge readiness
df_hb_velocity_aligned = df_hb_velocity.rename(columns={'mjd': 'time'})[
    ['time', 'blue_wing_flux', 'blue_wing_error', 'core_flux', 'core_error', 'red_wing_flux', 'red_wing_error']
]

# --- Critical Correction: Restrict data to MJD 47512–49255 ---
start_mjd_strict = 47512.0
end_mjd_strict = 49255.0

df_cont_aligned = df_cont_aligned[(df_cont_aligned['time'] >= start_mjd_strict) & (df_cont_aligned['time'] <= end_mjd_strict)]
df_hb_velocity_aligned = df_hb_velocity_aligned[(df_hb_velocity_aligned['time'] >= start_mjd_strict) & (df_hb_velocity_aligned['time'] <= end_mjd_strict)]

df_combined = pd.merge(df_cont_aligned, df_hb_velocity_aligned, on='time', how='outer', sort=True)

# Sort by time again just to be sure
df_combined = df_combined.sort_values(by='time').reset_index(drop=True)

# Set a common time grid (e.g., 1-day cadence, same as `dt_real` in previous SURD analysis)
dt_final = 1.0

# Define the full time range for the combined data, strictly using the corrected bounds
min_overall_time = start_mjd_strict
max_overall_time = end_mjd_strict

# Create a uniform time grid for resampling
uniform_time_grid = np.arange(min_overall_time, max_overall_time + dt_final, dt_final)

# Drop any original MJD points outside the strict range before interpolation to avoid issues
df_combined = df_combined[(df_combined['time'] >= min_overall_time) & (df_combined['time'] <= max_overall_time)]

columns_to_interpolate = [
    'cont_flux', 'cont_error',
    'blue_wing_flux', 'blue_wing_error',
    'core_flux', 'core_error',
    'red_wing_flux', 'red_wing_error'
]

prepared_data = pd.DataFrame({'time': uniform_time_grid})

for col in columns_to_interpolate:
    # Only interpolate within the existing data range, no ffill/bfill outside observed bounds
    series_to_interp = df_combined.set_index('time')[col].dropna()
    if not series_to_interp.empty:
        prepared_data[col] = np.interp(uniform_time_grid,
                                     series_to_interp.index,
                                     series_to_interp.values,
                                     left=np.nan, right=np.nan) # Prevent extrapolation
    else:
        prepared_data[col] = np.nan # If no data, all NaNs


# Drop any rows where all *flux* values are NaN (might happen if a series had no data within the common grid)
# We only drop if all four flux types (cont, blue, core, red) are NaN for that time point
flux_cols_to_check = ['cont_flux', 'blue_wing_flux', 'core_flux', 'red_wing_flux']
prepared_data = prepared_data.dropna(subset=flux_cols_to_check, how='all').reset_index(drop=True)

# Check if prepared_data is empty after dropping NaNs
if prepared_data.empty:
    raise ValueError("Prepared data is empty after interpolation and dropping NaNs. Check data ranges and overlap.")

# Z-score standardize the flux columns
# Ensure these columns exist and are not all NaN before z-scoring
for col_name in ['cont_flux', 'blue_wing_flux', 'core_flux', 'red_wing_flux']:
    if col_name in prepared_data.columns and not prepared_data[col_name].isnull().all():
        prepared_data[f'{col_name}_zscore'] = zscore(prepared_data[col_name])
    else:
        prepared_data[f'{col_name}_zscore'] = np.nan # Assign NaN if original was all NaN


print("--- Prepared Consistent Dataset Summary ---")
print(f"Time range: MJD {prepared_data['time'].min():.2f} to {prepared_data['time'].max():.2f}")
print(f"Number of epochs: {len(prepared_data)}")
print(f"Cadence: {dt_final} days")
display(prepared_data.head())

# Visual check of the prepared (interpolated and standardized) light curves
plt.figure(figsize=(16, 9))
plt.plot(prepared_data['time'], prepared_data['cont_flux_zscore'], label='Continuum (Standardized)', color='black', alpha=0.8)
plt.plot(prepared_data['time'], prepared_data['blue_wing_flux_zscore'], label='Blue Wing (Standardized)', color='blue', alpha=0.8)
plt.plot(prepared_data['time'], prepared_data['core_flux_zscore'], label='Core (Standardized)', color='green', alpha=0.8)
plt.plot(prepared_data['time'], prepared_data['red_wing_flux_zscore'], label='Red Wing (Standardized)', color='red', alpha=0.8)

plt.xlabel('Modified Julian Date (MJD)')
plt.ylabel('Standardized Flux')
plt.title('Prepared NGC 5548 Light Curves for SURD (Continuum + Hβ Velocity Bins)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()

print("The combined and standardized dataset is now ready for SURD analysis. It includes the continuum, blue wing, core, and red wing Hβ fluxes on a common time grid.")
print("For SURD, these will be arranged as: continuum (signal 1), blue wing (signal 2), core (signal 3), and red wing (signal 4).")
print("Depending on the specific SURD analysis, you would typically select 3 signals (e.g., cont, blue, core as drivers for red wing target, or cont as driver for blue, core, or red targets individually).")
--- Cell 48 ---
# @title
"""
SURD AGN starter script
=======================

Purpose
-------
This script collects the main steps testing
the SURD repository for an AGN-style reverberation mapping toy problem.

What this script does
---------------------
1. Assumes the SURD repository has already been cloned in Colab or locally.
2. Adds the SURD `utils` directory to Python's import path.
3. Imports the core SURD modules.
4. Builds a synthetic AGN-like time-series data set:
      - hidden UV driver
      - optical continuum proxy
      - X-ray continuum proxy
      - delayed Hβ emission-line response
5. Runs the SURD decomposition for a single lag.
6. Scans across many lags and stores metrics for the final target signal.
7. Provides a helper that returns results for *all* target signals.
8. Shows how to run one of the built-in analytic benchmark systems.

Important notes
---------------
- This script intentionally avoids `transport_map.py`, because that path
  requires the optional `mpart` / MParT dependency.
- The SURD `run()` helper prints and plots results for all targets, but it
  returns only the decomposition for the *last* target processed. For clean
  analysis, this script includes a custom `run_collect()` wrapper.
- The lag parameter `nlag` is in samples, not physical days.

Expected environment
--------------------
This script is written for the same setup that worked in Colab:
    pip install pymp-pypi numpy scipy matplotlib pandas scikit-learn

and then:
    sys.path.append("/content/SURD/utils")

Adapt the SURD_UTILS_PATH below if your cloned repo lives elsewhere.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Any

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# 1. Configure path to SURD repo utilities
# =============================================================================
# Change this if your SURD clone is somewhere else.
SURD_UTILS_PATH = Path("/content/SURD/utils")

if not SURD_UTILS_PATH.exists():
    raise FileNotFoundError(
        f"Could not find SURD utils directory at {SURD_UTILS_PATH}. "
        "Update SURD_UTILS_PATH to match your local or Colab clone."
    )

sys.path.append(str(SURD_UTILS_PATH))

# Core SURD imports
import surd            # noqa: E402
import it_tools        # noqa: E402  # imported for completeness / inspection
import analytic_eqs    # noqa: E402


# =============================================================================
# 2. Utility helpers
# =============================================================================
def zscore(x: np.ndarray) -> np.ndarray:
    """
    Standardise an array to zero mean and unit variance.

    Parameters
    ----------
    x : np.ndarray
        Input array.

    Returns
    -------
    np.ndarray
        Standardised array.
    """
    x = np.asarray(x, dtype=float)
    return (x - np.mean(x)) / np.std(x)


def run_collect(
    X: np.ndarray,
    nvars: int,
    nlag: int,
    nbins: int,
    axs: np.ndarray | None = None,
    print_results: bool = False,
) -> Dict[int, Dict[str, Any]]:
    """
    Run SURD for each target signal and return a clean dictionary of results.

    Why this helper exists
    ----------------------
    The repo's `surd.run(...)` function loops over all targets and plots all of
    them, but returns only the decomposition of the final target in the loop.
    For analysis, it is much more useful to collect results for every target.

    Parameters
    ----------
    X : np.ndarray
        Array of shape (nvars, ntimes). Each row is one signal.
    nvars : int
        Number of signals.
    nlag : int
        Lag in samples.
    nbins : int
        Number of histogram bins used by np.histogramdd.
    axs : np.ndarray or None
        Optional axes array of shape (nvars, 2) for SURD plotting.
    print_results : bool
        If True, print the SURD decomposition for each target.

    Returns
    -------
    dict
        Dictionary keyed by target index, with values containing:
            - "I_R"
            - "I_S"
            - "MI"
            - "info_leak"
    """
    results = {}

    for i in range(nvars):
        # Future of target i versus present/past of all variables
        Y = np.vstack([X[i, nlag:], X[:, :-nlag]])

        # Multi-dimensional histogram approximation to the joint distribution
        hist, _ = np.histogramdd(Y.T, nbins)

        # SURD decomposition
        I_R, I_S, MI, info_leak = surd.surd(hist)

        results[i] = {
            "I_R": I_R,
            "I_S": I_S,
            "MI": MI,
            "info_leak": info_leak,
        }

        if print_results:
            print(f"\nSURD CAUSALITY FOR TARGET SIGNAL {i + 1}")
            surd.nice_print(I_R, I_S, MI, info_leak)

        if axs is not None:
            surd.plot(I_R, I_S, info_leak, axs[i, :], nvars, threshold=-0.01)

    return results


# =============================================================================
# 3. Synthetic AGN-style toy data
# =============================================================================
def make_synthetic_agn(
    N: int = 3000,
    tau: int = 8,
    seed: int = 42,
    random_walk_driver: bool = True,
) -> np.ndarray:
    """
    Build a synthetic AGN-like system with three signals.

    Signal definitions
    ------------------
    1. optical continuum: noisy proxy for hidden UV driver
    2. X-ray continuum:  noisy proxy for hidden UV driver
    3. Hβ emission line: delayed response to the UV driver

    Parameters
    ----------
    N : int
        Number of time samples.
    tau : int
        Reverberation lag in samples.
    seed : int
        Random seed for reproducibility.
    random_walk_driver : bool
        If True, use a random-walk hidden driver (long memory).
        If False, use white noise hidden driver (sharper lag localisation).

    Returns
    -------
    np.ndarray
        Stacked array X of shape (3, N).
    """
    rng = np.random.default_rng(seed)

    # Hidden UV driver
    if random_walk_driver:
        uv = np.cumsum(rng.normal(0, 0.15, N))
    else:
        uv = rng.normal(0, 1.0, N)
    uv = zscore(uv)

    # Two observed continuum proxies
    optical = 0.9 * uv + 0.25 * rng.normal(size=N)
    xray = 0.7 * uv + 0.35 * rng.normal(size=N)

    # Delayed broad-line response
    hbeta = np.zeros(N)
    hbeta[tau:] = 0.8 * uv[:-tau] + 0.25 * rng.normal(size=N - tau)
    hbeta[:tau] = 0.25 * rng.normal(size=tau)

    optical = zscore(optical)
    xray = zscore(xray)
    hbeta = zscore(hbeta)

    return np.vstack([optical, xray, hbeta])


# =============================================================================
# 4. Demonstration: single-lag SURD run
# =============================================================================
def demo_single_lag() -> None:
    """
    Run a single-lag SURD decomposition on the synthetic AGN toy system.
    """
    tau = 8
    X = make_synthetic_agn(N=3000, tau=tau, seed=42, random_walk_driver=True)
    nvars = X.shape[0]

    fig, axs = plt.subplots(nvars, 2, figsize=(12, 4 * nvars), squeeze=False)

    # This is the repo's built-in wrapper.
    # It prints and plots all targets, but returns only the final target result.
    surd.run(X=X, nvars=nvars, nlag=tau, nbins=8, axs=axs)

    plt.tight_layout()
    plt.show()


# =============================================================================
# 5. Demonstration: lag scan for target Signal 3
# =============================================================================
def demo_lag_scan() -> None:
    """
    Scan over lag values and track selected information measures.

    We focus on the final target (Signal 3, Hβ) by using run_collect().
    """
    tau_true = 8
    X = make_synthetic_agn(N=3000, tau=tau_true, seed=42, random_walk_driver=True)
    nvars = X.shape[0]

    lags = np.arange(1, 21)
    mi1, mi2, syn12, leak = [], [], [], []

    for lag in lags:
        results = run_collect(X=X, nvars=nvars, nlag=lag, nbins=8)

        # Target index 2 = Signal 3 in 1-based human counting
        res = results[2]

        # Notes:
        # MI keys refer to variable combinations in the SURD histogram indexing.
        # We keep the same convention used during the Colab exploration.
        mi1.append(res["MI"].get((1,), np.nan))
        mi2.append(res["MI"].get((2,), np.nan))
        syn12.append(res["I_S"].get((1, 2), np.nan))
        leak.append(res["info_leak"])

    plt.figure(figsize=(8, 5))
    plt.plot(lags, mi1, marker="o", label="MI(1)")
    plt.plot(lags, mi2, marker="o", label="MI(2)")
    plt.plot(lags, syn12, marker="o", label="Synergy(1,2)")
    plt.axvline(tau_true, linestyle="--", label="True lag")
    plt.xlabel("Lag (samples)")
    plt.ylabel("Information")
    plt.title("Lag scan for target Signal 3 (Hβ)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(7, 4))
    plt.plot(lags, leak, marker="o")
    plt.axvline(tau_true, linestyle="--")
    plt.xlabel("Lag (samples)")
    plt.ylabel("Information leak")
    plt.title("Leak vs lag for target Signal 3 (Hβ)")
    plt.tight_layout()
    plt.show()

    print("\nInterpretation note:")
    print(
        "If you use a random-walk hidden driver, the information curves often "
        "look broad rather than sharply peaked, because the driver has strong "
        "temporal autocorrelation."
    )
    print(
        "For a cleaner lag-localisation test, rerun with "
        "random_walk_driver=False in make_synthetic_agn()."
    )


# =============================================================================
# 6. Demonstration: compare random-walk vs white-noise hidden driver
# =============================================================================
def demo_driver_comparison() -> None:
    """
    Compare how lag localisation changes when the hidden driver has long memory
    versus no memory.
    """
    tau_true = 8
    lags = np.arange(1, 21)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4), squeeze=False)

    for ax, use_rw, title in zip(
        axs[0],
        [True, False],
        ["Random-walk driver", "White-noise driver"],
        strict=True,
    ):
        X = make_synthetic_agn(
            N=3000,
            tau=tau_true,
            seed=42,
            random_walk_driver=use_rw,
        )

        leak = []
        for lag in lags:
            results = run_collect(X=X, nvars=3, nlag=lag, nbins=8)
            leak.append(results[2]["info_leak"])

        ax.plot(lags, leak, marker="o")
        ax.axvline(tau_true, linestyle="--")
        ax.set_xlabel("Lag (samples)")
        ax.set_ylabel("Information leak")
        ax.set_title(title)

    plt.tight_layout()
    plt.show()


# =============================================================================
# 7. Demonstration: built-in analytic benchmark
# =============================================================================
def demo_analytic_benchmark() -> None:
    """
    Run one of the built-in synthetic benchmark systems from analytic_eqs.

    Available functions discovered during the Colab session:
        - confounder
        - mediator
        - redundant_collider
        - synergistic_collider
    """
    q1, q2, q3 = analytic_eqs.confounder(3000)
    X = np.vstack([q1, q2, q3])

    fig, axs = plt.subplots(3, 2, figsize=(10, 12), squeeze=False)
    surd.run(X=X, nvars=3, nlag=1, nbins=8, axs=axs)
    plt.tight_layout()
    plt.show()


# =============================================================================
# 8. Main execution block
# =============================================================================
if __name__ == "__main__":
    print("Core SURD AGN starter script")
    print("-" * 40)

    # Uncomment whichever demos you want to run.

    # 1) Single-lag toy example
    demo_single_lag()

    # 2) Lag scan for the Hβ target
    demo_lag_scan()

    # 3) Compare long-memory and short-memory hidden drivers
    # demo_driver_comparison()

    # 4) Run a built-in benchmark motif
    # demo_analytic_benchmark()

    print("\nDone.")

--- Cell 49 ---
# @title
# =============================================================================
# 9. Controlled toy cases + real AGN preprocessing helpers
# =============================================================================
#
# This cell extends the original SURD AGN starter workflow with:
#
#   A) Controlled toy cases
#      1. One true driver
#      2. Redundant proxies
#      3. Genuinely synergistic drivers
#
#   B) A reusable lag-scan helper for target Signal 3
#
#   C) A plotting helper for lag-scan outputs
#
#   D) Simple real-AGN preprocessing utilities
#      1. resample to a uniform time grid
#      2. standardise each band
#      3. stack into X = np.vstack([...])
#
#   E) Example execution blocks for both toy data and real AGN data
#
# Assumptions
# -----------
# This script assumes the earlier notebook/script has already defined:
#   - np
#   - plt
#   - zscore(...)
#   - run_collect(...)
#
# In the SURD setup we used:
#   - row 0 = optical continuum
#   - row 1 = X-ray continuum
#   - row 2 = Hβ emission line
#
# and target Signal 3 corresponds to index 2 in Python.
# =============================================================================


# -----------------------------------------------------------------------------
# A1. Controlled toy case: one true driver
# -----------------------------------------------------------------------------
# Physical idea:
#   - Hβ responds only to the optical continuum after a lag tau
#   - X-ray is unrelated noise
#
# Expected SURD behaviour:
#   - optical should show strong unique information near the true lag
#   - X-ray should show very little unique information
#   - redundancy and synergy should be small
# -----------------------------------------------------------------------------
def make_case_one_true_driver(N=3000, tau=8, seed=1):
    rng = np.random.default_rng(seed)

    optical = rng.normal(size=N)
    xray = rng.normal(size=N)   # unrelated noise

    hbeta = np.zeros(N)
    hbeta[tau:] = 0.9 * optical[:-tau] + 0.2 * rng.normal(size=N - tau)
    hbeta[:tau] = 0.2 * rng.normal(size=tau)

    optical = zscore(optical)
    xray = zscore(xray)
    hbeta = zscore(hbeta)

    return np.vstack([optical, xray, hbeta])


# -----------------------------------------------------------------------------
# A2. Controlled toy case: redundant proxies
# -----------------------------------------------------------------------------
# Physical idea:
#   - optical and X-ray are both noisy proxies of a hidden UV driver
#   - Hβ responds to that hidden UV driver after a lag tau
#
# Expected SURD behaviour:
#   - redundancy should be elevated
#   - both observed bands can appear informative
#   - some information leak can remain because the true UV driver is hidden
#
# The random_walk option is useful to study long-memory drivers, which often
# produce broader lag curves instead of sharply localised lag peaks.
# -----------------------------------------------------------------------------
def make_case_redundant_proxies(N=3000, tau=8, seed=2, random_walk=False):
    rng = np.random.default_rng(seed)

    if random_walk:
        uv = np.cumsum(rng.normal(0, 0.15, N))
    else:
        uv = rng.normal(size=N)

    uv = zscore(uv)

    optical = 0.9 * uv + 0.25 * rng.normal(size=N)
    xray = 0.8 * uv + 0.25 * rng.normal(size=N)

    hbeta = np.zeros(N)
    hbeta[tau:] = 0.85 * uv[:-tau] + 0.2 * rng.normal(size=N - tau)
    hbeta[:tau] = 0.2 * rng.normal(size=tau)

    optical = zscore(optical)
    xray = zscore(xray)
    hbeta = zscore(hbeta)

    return np.vstack([optical, xray, hbeta])


# -----------------------------------------------------------------------------
# A3. Controlled toy case: genuinely synergistic drivers
# -----------------------------------------------------------------------------
# Physical idea:
#   - Hβ depends on a nonlinear combination of optical and X-ray
#   - neither optical nor X-ray alone is enough
#   - both together matter
#
# Expected SURD behaviour:
#   - synergy should be noticeably larger
#   - individual-band information may not fully explain the target
# -----------------------------------------------------------------------------
def make_case_synergy(N=3000, tau=8, seed=3):
    rng = np.random.default_rng(seed)

    optical = rng.normal(size=N)
    xray = rng.normal(size=N)

    hbeta = np.zeros(N)
    signal = optical[:-tau] * xray[:-tau]   # nonlinear joint driver
    hbeta[tau:] = signal + 0.25 * rng.normal(size=N - tau)
    hbeta[:tau] = 0.25 * rng.normal(size=tau)

    optical = zscore(optical)
    xray = zscore(xray)
    hbeta = zscore(hbeta)

    return np.vstack([optical, xray, hbeta])


# -----------------------------------------------------------------------------
# B. Reusable lag-scan helper for target Signal 3 (Hβ)
# -----------------------------------------------------------------------------
# This helper scans a list/array of lags and stores selected SURD quantities
# for the target row corresponding to the emission-line signal.
#
# Stored metrics:
#   - info_leak
#   - MI(1), MI(2)
#   - U(1), U(2)  [approx. from I_R singletons as used in repo printout]
#   - R(1,2)
#   - S(1,2)
# -----------------------------------------------------------------------------
def lag_scan_target3(X, lags, nbins=8):
    metrics = {
        "lag": [],
        "info_leak": [],
        "MI1": [],
        "MI2": [],
        "U1": [],
        "U2": [],
        "R12": [],
        "S12": [],
    }

    for lag in lags:
        results = run_collect(X=X, nvars=3, nlag=lag, nbins=nbins)
        res = results[2]   # target index 2 = Signal 3 (Hβ)

        metrics["lag"].append(lag)
        metrics["info_leak"].append(res["info_leak"])
        metrics["MI1"].append(res["MI"].get((1,), np.nan))
        metrics["MI2"].append(res["MI"].get((2,), np.nan))
        metrics["U1"].append(res["I_R"].get((1,), np.nan))
        metrics["U2"].append(res["I_R"].get((2,), np.nan))
        metrics["R12"].append(res["I_R"].get((1, 2), np.nan))
        metrics["S12"].append(res["I_S"].get((1, 2), np.nan))

    return metrics


# -----------------------------------------------------------------------------
# C. Plotting helper for lag-scan results
# -----------------------------------------------------------------------------
# Produces:
#   1. Information curves vs lag
#   2. Information leak vs lag
#
# Optionally mark the true lag with a vertical dashed line.
# -----------------------------------------------------------------------------
def plot_lag_metrics(metrics, true_tau=None, title="Lag scan for target Signal 3", show_plot=True, save_path=None):
    lags = np.array(metrics["lag"])

    plt.figure(figsize=(10, 6))
    plt.plot(lags, metrics["MI1"], marker="o", label="MI(1)")
    plt.plot(lags, metrics["MI2"], marker="o", label="MI(2)")
    plt.plot(lags, metrics["R12"], marker="o", label="Redundancy(1,2)")
    plt.plot(lags, metrics["S12"], marker="o", label="Synergy(1,2)")
    if true_tau is not None:
        plt.axvline(true_tau, linestyle="--", label="True lag")
    plt.xlabel("Lag")
    plt.ylabel("Information")
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6) # Added grid for consistency
    plt.tight_layout()
    if save_path: plt.savefig(save_path + "_info.png")
    if show_plot: plt.show()
    plt.close() # Close figure to free memory

    plt.figure(figsize=(7, 4))
    plt.plot(lags, metrics["info_leak"], marker="o")
    if true_tau is not None:
        plt.axvline(true_tau, linestyle="--")
    plt.xlabel("Lag")
    plt.ylabel("Information leak")
    plt.title(f"{title}: leak")
    plt.grid(True, linestyle='--', alpha=0.6) # Added grid for consistency
    plt.tight_layout()
    if save_path: plt.savefig(save_path + "_leak.png")
    if show_plot: plt.show()
    plt.close() # Close figure to free memory


# -----------------------------------------------------------------------------
# D1. Simple interpolation onto a uniform time grid
# -----------------------------------------------------------------------------
# Why this is needed:
#   The SURD workflow here expects aligned time-series samples.
#   Real AGN light curves are often irregularly sampled, so as a first-pass
#   prototype we interpolate each band onto a common uniform grid.
#
# Important caution:
#   Interpolation can introduce artificial structure. It is acceptable for an
#   initial prototype, but later you may want a more careful uneven-sampling
#   treatment or at least sensitivity tests with different grid spacings.
# -----------------------------------------------------------------------------
def resample_to_uniform_grid(times, values, tmin=None, tmax=None, dt=1.0):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)

    if tmin is None:
        tmin = np.min(times)
    if tmax is None:
        tmax = np.max(times)

    grid = np.arange(tmin, tmax + dt, dt)
    vals = np.interp(grid, times, values)
    return grid, vals


# -----------------------------------------------------------------------------
# D2. Prepare real AGN time series for SURD
# -----------------------------------------------------------------------------
# Inputs:
#   - time_opt, flux_opt   : optical continuum light curve
#   - time_xray, flux_xray : X-ray continuum light curve
#   - time_hb, flux_hb     : Hβ emission-line light curve
#
# Steps:
#   1. restrict to the common overlapping time range
#   2. interpolate each series onto a shared uniform grid
#   3. z-score each band
#   4. stack into X = np.vstack([optical, xray, hbeta])
#
# Returns:
#   - grid : common time grid
#   - X    : stacked standardised data for SURD
# -----------------------------------------------------------------------------
def prepare_agn_timeseries(time_opt, flux_opt, time_xray, flux_xray, time_hb, flux_hb, dt=1.0):
    tmin = max(np.min(time_opt), np.min(time_xray), np.min(time_hb))
    tmax = min(np.max(time_opt), np.max(time_xray), np.max(time_hb))

    grid, opt = resample_to_uniform_grid(time_opt, flux_opt, tmin=tmin, tmax=tmax, dt=dt)
    _, xray = resample_to_uniform_grid(time_xray, flux_xray, tmin=tmin, tmax=tmax, dt=dt)
    _, hb = resample_to_uniform_grid(time_hb, flux_hb, tmin=tmin, tmax=tmax, dt=dt)

    opt = zscore(opt)
    xray = zscore(xray)
    hb = zscore(hb)

    X = np.vstack([opt, xray, hb])
    return grid, X


# -----------------------------------------------------------------------------
# E1. Example: run all three toy benchmark cases
# -----------------------------------------------------------------------------
# These are the first things to run before trusting any real-data result.
# They tell you whether SURD is qualitatively recovering the expected pattern:
#
#   - one true driver      -> strong unique info from optical
#   - redundant proxies    -> elevated redundancy
#   - synergistic drivers  -> elevated synergy
# -----------------------------------------------------------------------------
# lags = np.arange(1, 21)
# tau = 8

# # --- Case 1: one true driver ---
# X_true = make_case_one_true_driver(N=3000, tau=tau, seed=1)
# m_true = lag_scan_target3(X_true, lags, nbins=8)
# plot_lag_metrics(m_true, true_tau=tau, title="Controlled toy case: one true driver")

# # --- Case 2: redundant proxies ---
# X_red = make_case_redundant_proxies(N=3000, tau=tau, seed=2, random_walk=False)
# m_red = lag_scan_target3(X_red, lags, nbins=8)
# plot_lag_metrics(m_red, true_tau=tau, title="Controlled toy case: redundant proxies")

# # --- Case 3: genuinely synergistic drivers ---
# X_syn = make_case_synergy(N=3000, tau=tau, seed=3)
# m_syn = lag_scan_target3(X_syn, lags, nbins=8)
# plot_lag_metrics(m_syn, true_tau=tau, title="Controlled toy case: synergistic drivers")


# -----------------------------------------------------------------------------
# E2. Example: compare the three toy cases at the true lag
# -----------------------------------------------------------------------------
# This gives a quick numerical comparison of the main information terms.
# -----------------------------------------------------------------------------
# true_lag_index = np.where(lags == tau)[0][0]

# print("Summary at the true lag:")
# print("-" * 60)
# print(f"One true driver     : MI1={m_true['MI1'][true_lag_index]:.4f}, "
#       f"MI2={m_true['MI2'][true_lag_index]:.4f}, "
#       f"R12={m_true['R12'][true_lag_index]:.4f}, "
#       f"S12={m_syn['S12'][true_lag_index]:.4f}, "
#       f"Leak={m_true['info_leak'][true_lag_index]:.4f}")

# print(f"Redundant proxies   : MI1={m_red['MI1'][true_lag_index]:.4f}, "
#       f"MI2={m_red['MI2'][true_lag_index]:.4f}, "
#       f"R12={m_red['R12'][true_lag_index]:.4f}, "
#       f"S12={m_red['S12'][true_lag_index]:.4f}, "
#       f"Leak={m_red['info_leak'][true_lag_index]:.4f}")

# print(f"Synergistic drivers : MI1={m_syn['MI1'][true_lag_index]:.4f}, "
#       f"MI2={m_syn['MI2'][true_lag_index]:.4f}, "
#       f"R12={m_syn['R12'][true_lag_index]:.4f}, "
#       f"S12={m_syn['S12'][true_lag_index]:.4f}, "
#       f"Leak={m_syn['info_leak'][true_lag_index]:.4f}")


# -----------------------------------------------------------------------------
# E3. Template for real AGN data
# -----------------------------------------------------------------------------
# Replace the placeholder arrays below with your actual light-curve data.
#
# Example expected variables:
#   time_opt, flux_opt
#   time_xray, flux_xray
#   time_hbeta, flux_hbeta
#
# Then:
#   1. resample to a uniform cadence
#   2. stack into X_real
#   3. scan over lags
#   4. inspect unique/redundant/synergistic information for the Hβ target
#
# To use this section, uncomment it and replace the placeholders.
# -----------------------------------------------------------------------------
"""
# ----- Example placeholders: replace with your real AGN arrays -----
time_opt = np.array([...], dtype=float)
flux_opt = np.array([...], dtype=float)

time_xray = np.array([...], dtype=float)
flux_xray = np.array([...], dtype=float)

time_hbeta = np.array([...], dtype=float)
flux_hbeta = np.array([...], dtype=float)

# Choose a uniform grid spacing in the same units as your time arrays.
# Example: if times are in days, dt = 1.0 means 1-day cadence.
dt = 1.0

# Prepare the three aligned standardised bands
grid, X_real = prepare_agn_timeseries(
    time_opt, flux_opt,
    time_xray, flux_xray,
    time_hbeta, flux_hbeta,
    dt=dt
)

# Visual sanity check of the interpolated / standardised signals
plt.figure(figsize=(12, 5))
plt.plot(grid, X_real[0], label='Optical')
plt.plot(grid, X_real[1], label='X-ray')
plt.plot(grid, X_real[2], label='Hβ')
plt.xlabel('Time')
plt.ylabel('Standardised flux')
plt.title('Prepared AGN time series on common grid')
plt.legend()
plt.tight_layout()
plt.show()

# Define lag range in samples.
# If dt = 2 days and you want to scan to 100 days, use:
# lags = np.arange(1, int(100 / dt) + 1)
lags_real = np.arange(1, 51)

# Run SURD lag scan on the emission-line target
metrics_real = lag_scan_target3(X_real, lags_real, nbins=8)

# Plot the results
plot_lag_metrics(metrics_real, true_tau=None, title='Real AGN lag scan')

# Optional: convert lag samples to physical time
lag_days = dt * np.array(metrics_real['lag'])

# Example extra plot with physical units
plt.figure(figsize=(10, 6))
plt.plot(lag_days, metrics_real["MI1"], marker="o", label="MI(1)")
plt.plot(lag_days, metrics_real["MI2"], marker="o", label="MI(2)")
plt.plot(lag_days, metrics_real["R12"], marker="o", label="Redundancy(1,2)")
plt.plot(lag_days, metrics_real["S12"], marker="o", label="Synergy(1,2)")
plt.xlabel("Lag (days)")
plt.ylabel("Information")
plt.title("Real AGN lag scan in physical units")
plt.legend()
plt.tight_layout()
plt.show()

# Interpretation:
#   - strong U1 / MI1 near some lag: optical carries predictive line information
#   - strong U2 / MI2 near some lag: X-ray carries predictive line information
#   - strong R12: both bands may be redundant proxies for a hidden driver
#   - strong S12: combining both bands helps more than either alone
#   - high leak: substantial line variability remains unexplained by observed bands.
"""
# =============================================================================
# End of extension cell
# =============================================================================

--- Cell 53 ---
# @title
# Choose a uniform grid spacing (e.g., 1 day cadence)
dt = 1.0

# Prepare the three aligned standardised bands using the helper function
grid, X_real_simulated = prepare_agn_timeseries(
    time_opt, flux_opt,
    time_xray, flux_xray,
    time_hbeta, flux_hbeta,
    dt=dt
)

# --- Visual sanity check of the interpolated / standardised signals ---
plt.figure(figsize=(12, 5))
plt.plot(grid, X_real_simulated[0], label='Optical (interpolated & standardized)')
plt.plot(grid, X_real_simulated[1], label='X-ray (interpolated & standardized)')
plt.plot(grid, X_real_simulated[2], label='Hβ (interpolated & standardized)')
plt.xlabel('Time (days)')
plt.ylabel('Standardized flux')
plt.title('Prepared Simulated AGN Time Series on Common Grid')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

--- Cell 60 ---
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Ensure zscore is defined (from cell '_tkPv6RIVMl5')
# Ensure lag_scan_target3 and plot_lag_metrics are defined (from cell '0cxiu-N1VNSO')

def zscore(x: np.ndarray) -> np.ndarray:
    """
    Standardise an array to zero mean and unit variance.
    Handles cases where standard deviation is zero and NaN values.
    """
    x = np.asarray(x, dtype=float)
    # Use nanmean and nanstd to ignore NaNs when calculating statistics
    mean_val = np.nanmean(x)
    std_dev = np.nanstd(x)

    if std_dev == 0 or np.isnan(std_dev):
        # If std_dev is zero or NaN (e.g., if x was all NaNs or constant non-NaNs),
        # return an array of zeros.
        z_scored = np.zeros_like(x)
    else:
        z_scored = (x - mean_val) / std_dev

    # Replace any remaining NaNs in the output (corresponding to NaNs in original x) with 0.0
    z_scored[np.isnan(z_scored)] = 0.0
    return z_scored


def prepare_combined_data(df_matching_cont, df_hb_velocity, dt_cadence):
    """
    Prepares and standardizes the combined continuum and Hβ velocity-resolved data
    onto a uniform time grid with a specified cadence, strictly within the observed Hβ range.

    Args:
        df_matching_cont (pd.DataFrame): DataFrame for the matching continuum.
        df_hb_velocity (pd.DataFrame): DataFrame for velocity-binned Hβ.
        dt_cadence (float): The desired cadence (time step) for the uniform grid.

    Returns:
        pd.DataFrame: A DataFrame containing the prepared and standardized data.
    """
    df_cont_aligned = df_matching_cont.rename(columns={'mjd': 'time', 'flux': 'cont_flux', 'error': 'cont_error'})[['time', 'cont_flux', 'cont_error']]
    df_hb_velocity_aligned = df_hb_velocity.rename(columns={'mjd': 'time'})[
        ['time', 'blue_wing_flux', 'blue_wing_error', 'core_flux', 'core_error', 'red_wing_flux', 'red_wing_error']
    ]

    # Restrict data to the common observed overlap between MJD 47512 and 49255
    start_mjd = 47512.0
    end_mjd = 49255.0

    df_cont_aligned = df_cont_aligned[(df_cont_aligned['time'] >= start_mjd) & (df_cont_aligned['time'] <= end_mjd)]
    df_hb_velocity_aligned = df_hb_velocity_aligned[(df_hb_velocity_aligned['time'] >= start_mjd) & (df_hb_velocity_aligned['time'] <= end_mjd)]

    df_combined = pd.merge(df_cont_aligned, df_hb_velocity_aligned, on='time', how='outer', sort=True)
    df_combined = df_combined.sort_values(by='time').reset_index(drop=True)

    # The uniform time grid should also be within the specified bounds
    uniform_time_grid = np.arange(start_mjd, end_mjd + dt_cadence, dt_cadence)

    columns_to_interpolate = [
        'cont_flux', 'cont_error',
        'blue_wing_flux', 'blue_wing_error',
        'core_flux', 'core_error',
        'red_wing_flux', 'red_wing_error'
    ]

    prepared_data_cadence = pd.DataFrame({'time': uniform_time_grid})

    for col in columns_to_interpolate:
        series_to_interp = df_combined.set_index('time')[col].dropna()
        if not series_to_interp.empty:
            prepared_data_cadence[col] = np.interp(uniform_time_grid,
                                                 series_to_interp.index,
                                                 series_to_interp.values,
                                                 left=np.nan, right=np.nan)
        else:
            prepared_data_cadence[col] = np.nan

    flux_cols_to_check = ['cont_flux', 'blue_wing_flux', 'core_flux', 'red_wing_flux']
    prepared_data_cadence = prepared_data_cadence.dropna(subset=flux_cols_to_check, how='any').reset_index(drop=True)

    if prepared_data_cadence.empty:
        raise ValueError(f"Prepared data is empty after interpolation and dropping NaNs for dt_cadence={dt_cadence}. Check data ranges and overlap.")

    for col_name in flux_cols_to_check:
        if col_name in prepared_data_cadence.columns and not prepared_data_cadence[col_name].isnull().all():
            prepared_data_cadence[f'{col_name}_zscore'] = zscore(prepared_data_cadence[col_name])
        else:
            prepared_data_cadence[f'{col_name}_zscore'] = np.nan # If no data, assign NaN

    return prepared_data_cadence

def run_surd_for_target_and_plot(X_array, lags, nbins, target_name, plot_title_prefix="", true_tau=None, show_plot=True):
    """
    Runs SURD lag scan for a given X_array (drivers + target) and plots the metrics.
    """
    print(f"\n--- Running SURD for {target_name} with nbins={nbins} ---")
    metrics = lag_scan_target3(X_array, lags, nbins=nbins)
    plot_lag_metrics(
        metrics,
        true_tau=true_tau,
        title=f"{plot_title_prefix}SURD Lag Scan: {target_name}",
        show_plot=show_plot
    )
    return metrics

def run_surd_shuffle_test_extract_metrics(
    X_array: np.ndarray,
    lags: np.ndarray,
    nbins: int,
    num_shuffles: int,
    shuffled_signal_index: int
) -> dict:
    """
    Performs SURD lag scans on shuffled data (using circular time shifts) and extracts relevant statistics
    (median, 2.5th, 16th, 84th, 97.5th percentiles) for synergy and information leak.

    Args:
        X_array (np.ndarray): The original stacked array of signals (e.g., [Driver1, Driver2, Target]).
        lags (np.ndarray): Array of lag values to scan.
        nbins (int): Number of histogram bins for SURD.
        num_shuffles (int): Number of shuffle iterations to perform.
        shuffled_signal_index (int): Index of the signal in X_array to shuffle (0 for S1, 1 for S2).

    Returns:
        dict: A dictionary containing median, 2.5th, 16th, 84th, and 97.5th percentiles
              for synergy and information leak across all shuffles, for each lag.
    """
    all_shuffled_s12 = []
    all_shuffled_leak = []

    print(f"Running {num_shuffles} circular time shift iterations for signal {shuffled_signal_index+1}...")

    for i in range(num_shuffles):
        X_shuffled = X_array.copy()
        shift = np.random.randint(X_shuffled.shape[1])
        X_shuffled[shuffled_signal_index] = np.roll(X_shuffled[shuffled_signal_index], shift)

        shuffled_metrics = lag_scan_target3(X_shuffled, lags, nbins=nbins)

        all_shuffled_s12.append(shuffled_metrics['S12'])
        all_shuffled_leak.append(shuffled_metrics['info_leak'])

    all_shuffled_s12 = np.array(all_shuffled_s12)
    all_shuffled_leak = np.array(all_shuffled_leak)

    synergy_stats = {
        'median_synergy': np.median(all_shuffled_s12, axis=0),
        'p2_5_synergy': np.percentile(all_shuffled_s12, 2.5, axis=0),
        'p16_synergy': np.percentile(all_shuffled_s12, 16, axis=0),
        'p84_synergy': np.percentile(all_shuffled_s12, 84, axis=0),
        'p97_5_synergy': np.percentile(all_shuffled_s12, 97.5, axis=0)
    }
    leak_stats = {
        'median_leak': np.median(all_shuffled_leak, axis=0),
        'p2_5_leak': np.percentile(all_shuffled_leak, 2.5, axis=0),
        'p16_leak': np.percentile(all_shuffled_leak, 16, axis=0),
        'p84_leak': np.percentile(all_shuffled_leak, 84, axis=0),
        'p97_5_leak': np.percentile(all_shuffled_leak, 97.5, axis=0)
    }

    results = {
        'lags': lags,
        **synergy_stats,
        **leak_stats
    }
    return results

def run_surd_block_shuffle_test_extract_metrics(
    X_array: np.ndarray,
    lags: np.ndarray,
    nbins: int,
    num_shuffles: int,
    shuffled_signal_index: int,
    block_size_days: float,
    dt_cadence: float
) -> dict:
    """
    Performs SURD lag scans on block-shuffled data and extracts relevant statistics
    (median, 2.5th, 16th, 84th, 97.5th percentiles) for synergy and information leak.
    Block shuffling preserves autocorrelation within blocks.
    """
    all_shuffled_s12 = []
    all_shuffled_leak = []

    n_samples = X_array.shape[1]
    block_size_samples = max(1, int(block_size_days / dt_cadence))
    num_blocks = n_samples // block_size_samples

    print(f"Running {num_shuffles} block shuffle iterations (block_size={block_size_days}d) for signal {shuffled_signal_index+1}...")

    for i in range(num_shuffles):
        X_shuffled = X_array.copy()
        signal_to_shuffle = X_shuffled[shuffled_signal_index]

        blocks = [signal_to_shuffle[j*block_size_samples:(j+1)*block_size_samples] for j in range(num_blocks)]
        remaining = signal_to_shuffle[num_blocks*block_size_samples:]

        np.random.shuffle(blocks)

        shuffled_signal_shuffled_blocks = np.concatenate(blocks)
        shuffled_signal_reconstructed = np.concatenate([shuffled_signal_shuffled_blocks, remaining])

        X_shuffled[shuffled_signal_index] = shuffled_signal_reconstructed

        shuffled_metrics = lag_scan_target3(X_shuffled, lags, nbins=nbins)

        all_shuffled_s12.append(shuffled_metrics['S12'])
        all_shuffled_leak.append(shuffled_metrics['info_leak'])

    all_shuffled_s12 = np.array(all_shuffled_s12)
    all_shuffled_leak = np.array(all_shuffled_leak)

    synergy_stats = {
        'median_synergy': np.median(all_shuffled_s12, axis=0),
        'p2_5_synergy': np.percentile(all_shuffled_s12, 2.5, axis=0),
        'p16_synergy': np.percentile(all_shuffled_s12, 16, axis=0),
        'p84_synergy': np.percentile(all_shuffled_s12, 84, axis=0),
        'p97_5_synergy': np.percentile(all_shuffled_s12, 97.5, axis=0)
    }
    leak_stats = {
        'median_leak': np.median(all_shuffled_leak, axis=0),
        'p2_5_leak': np.percentile(all_shuffled_leak, 2.5, axis=0),
        'p16_leak': np.percentile(all_shuffled_leak, 16, axis=0),
        'p84_leak': np.percentile(all_shuffled_leak, 84, axis=0),
        'p97_5_leak': np.percentile(all_shuffled_leak, 97.5, axis=0)
    }

    results = {
        'lags': lags,
        **synergy_stats,
        **leak_stats
    }
    return results
