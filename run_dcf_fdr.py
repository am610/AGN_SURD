import os
import numpy as np
import pandas as pd
from scipy.stats import zscore

# ----------------- 1. LOAD AND PREPARE DATA (Same as V11) -----------------
print("Loading data...")
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

# Load Continuum
# c5100.dat is space-separated with columns: JD (which is JD-2440000), flux, error
df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd_2440000', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd_2440000'] # actually the notebook maps jd_2440000 directly to mjd

# Load Spectroscopic Bins
df_hb = pd.read_csv(hb_bins_path)
# In df_hb, column 'time' is MJD

# Enforce Strict Overlap Range (MJD 47512 to 49255)
tmin = max(df_cont['mjd'].min(), df_hb['mjd'].min())
tmax = min(df_cont['mjd'].max(), df_hb['mjd'].max())
print(f"Strict overlap MJD range: {tmin} to {tmax}")

# Align on a uniform 1.0-day grid
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
print(f"Prepared uniform dataset size: {prepared_data.shape}")

# ----------------- 2. COMPUTE DISCRETE CORRELATION FUNCTION (DCF) -----------------
# DCF is calculated for irregularly sampled or regularly sampled data.
# Since our prepared data is regularly sampled, the DCF is equivalent to the standard CCF
# but we compute it explicitly to match standard practices.
def compute_dcf(x, y, time, max_lag=120, bin_width=1.0):
    lags = np.arange(1.0, max_lag + bin_width, bin_width)
    dcf_values = []
    
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    var_x = np.var(x)
    var_y = np.var(y)
    
    for lag in lags:
        # We look for pairs (x(t), y(t + lag))
        # Since it is regularly sampled with dt=1.0, a lag of 'lag' days corresponds to 'lag' indices.
        lag_idx = int(round(lag))
        if lag_idx < len(x):
            # Compute correlation for overlapping parts
            cov = np.mean((x[:-lag_idx] - mean_x) * (y[lag_idx:] - mean_y))
            val = cov / np.sqrt(var_x * var_y)
            dcf_values.append(val)
        else:
            dcf_values.append(np.nan)
            
    dcf_values = np.array(dcf_values)
    peak_idx = np.nanargmax(dcf_values)
    peak_lag = lags[peak_idx]
    peak_val = dcf_values[peak_idx]
    return lags, dcf_values, peak_lag, peak_val

print("\n--- Computing DCF (Continuum vs. Hbeta components) ---")
lags, dcf_blue, peak_lag_blue, peak_val_blue = compute_dcf(
    prepared_data['cont_flux_zscore'].values, prepared_data['blue_wing_flux_zscore'].values, prepared_data['time'].values
)
print(f"Blue Wing DCF: Peak Lag = {peak_lag_blue} days, Correlation = {peak_val_blue:.4f}")

_, dcf_core, peak_lag_core, peak_val_core = compute_dcf(
    prepared_data['cont_flux_zscore'].values, prepared_data['core_flux_zscore'].values, prepared_data['time'].values
)
print(f"Core Hbeta DCF: Peak Lag = {peak_lag_core} days, Correlation = {peak_val_core:.4f}")

_, dcf_red, peak_lag_red, peak_val_red = compute_dcf(
    prepared_data['cont_flux_zscore'].values, prepared_data['red_wing_flux_zscore'].values, prepared_data['time'].values
)
print(f"Red Wing DCF: Peak Lag = {peak_lag_red} days, Correlation = {peak_val_red:.4f}")

# ----------------- 3. FALSE DISCOVERY RATE (FDR) CORRECTION -----------------
# We load the shuffle results and compute the p-value at each lag:
# p-value = (number of surrogate runs with synergy >= real synergy) / total surrogate runs
# Since the raw circular and block shuffles were done and saved to CSVs:
# We will read these CSVs and do a statistical analysis.
def apply_fdr_for_component(csv_path, label):
    if not os.path.exists(csv_path):
        print(f"CSV file not found: {csv_path}")
        return
    df = pd.read_csv(csv_path)
    
    # We don't have the individual surrogate runs in the summary CSV (only median, p5, p95).
    # However, to simulate p-value estimation:
    # Under the null hypothesis of no coupling, synergy follows a normal-like distribution.
    # We can estimate the z-score of the real synergy compared to the shuffled median and 95th percentile.
    # For a normal distribution, the 95th percentile is median + 1.645 * sigma.
    # So sigma = (p95 - median) / 1.645
    # Then z-score = (real - median) / sigma
    # p-value = 1 - cdf(z-score)
    # Let's compute this for each lag:
    from scipy.stats import norm
    
    p_values = []
    for idx, row in df.iterrows():
        real = row['real_synergy']
        med = row['median_synergy_cont_shuffle']
        p95 = row['p95_synergy_cont_shuffle']
        
        # Avoid division by zero
        sigma = (p95 - med) / 1.645 if p95 > med else 1e-5
        if sigma <= 0:
            sigma = 1e-5
            
        z = (real - med) / sigma
        p = 1.0 - norm.cdf(z)
        p_values.append(max(p, 1e-15)) # avoid absolute zero
        
    df['p_value'] = p_values
    
    # Benjamini-Hochberg FDR correction
    n = len(df)
    sorted_df = df.sort_values('p_value').copy()
    sorted_df['rank'] = np.arange(1, n + 1)
    sorted_df['fdr_threshold'] = (sorted_df['rank'] / n) * 0.05
    sorted_df['significant'] = sorted_df['p_value'] <= sorted_df['fdr_threshold']
    
    # The largest rank i for which p_i <= (i/n)*alpha is the cutoff.
    sig_df = sorted_df[sorted_df['significant']]
    if not sig_df.empty:
        max_sig_rank = sig_df['rank'].max()
        cutoff_p = sorted_df[sorted_df['rank'] == max_sig_rank]['p_value'].values[0]
        df['significant_fdr'] = df['p_value'] <= cutoff_p
    else:
        df['significant_fdr'] = False
        cutoff_p = 0.0
        
    sig_lags = df[df['significant_fdr']]['lag'].values
    print(f"\n--- FDR Correction (q=0.05) for {label} ---")
    print(f"Total lags scanned: {n}")
    print(f"Significant lags: {len(sig_lags)} / {n}")
    if len(sig_lags) > 0:
        print(f"Significant lag ranges: {sig_lags}")
    else:
        print("No lags are statistically significant after Benjamini-Hochberg FDR correction (q=0.05).")
    return df

core_csv = "/Users/ayan/Programs/SURD/agn_surd_project/plots/test_c_shuffle_results/robustness_surrogate_core.csv"
red_csv = "/Users/ayan/Programs/SURD/agn_surd_project/plots/test_c_shuffle_results/robustness_surrogate_red.csv"
blue_csv = "/Users/ayan/Programs/SURD/agn_surd_project/plots/test_c_shuffle_results/robustness_surrogate_blue.csv"

df_core_fdr = apply_fdr_for_component(core_csv, "Core Hbeta")
df_red_fdr = apply_fdr_for_component(red_csv, "Red Wing Hbeta")
df_blue_fdr = apply_fdr_for_component(blue_csv, "Blue Wing Hbeta")
