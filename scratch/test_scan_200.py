import sys
import numpy as np
import pandas as pd
from scipy.stats import zscore

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

c_z = zscore(np.interp(mjd_grid, df_cont['jd'], df_cont['flux']))
b_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['blue_wing_flux']))
y_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['core_flux']))
r_z = zscore(np.interp(mjd_grid, df_hb['mjd'], df_hb['red_wing_flux']))

def run_collect_2pred(target_arr, pred1_arr, pred2_arr, nlag, nbins=8):
    future_target = target_arr[nlag:]
    pred_1 = pred1_arr[:-nlag]
    pred_2 = pred2_arr[:-nlag]
    Y = np.vstack([future_target, pred_1, pred_2])
    hist, _ = np.histogramdd(Y.T, nbins)
    hist = hist / np.sum(hist)
    I_R, I_S, MI, info_leak = surd.surd(hist)
    joint_mi = MI.get((1, 2), 1e-14)
    if joint_mi < 1e-14: joint_mi = 1e-14
    return I_S.get((1, 2), 0.0) / joint_mi

lags = np.arange(1, 201)
print("Scanning Core...")
syn_core = [run_collect_2pred(y_z, c_z, b_z, l) for l in lags]
print("Scanning Red...")
syn_red = [run_collect_2pred(r_z, c_z, b_z, l) for l in lags]
print("Scanning Blue...")
syn_blue = [run_collect_2pred(b_z, c_z, y_z, l) for l in lags]

print("\n--- Peaks up to 200 days ---")
print(f"Core Peak: Lag={lags[np.argmax(syn_core)]}d, Val={max(syn_core):.4f}")
print(f"Red Peak:  Lag={lags[np.argmax(syn_red)]}d, Val={max(syn_red):.4f}")
print(f"Blue Peak: Lag={lags[np.argmax(syn_blue)]}d, Val={max(syn_blue):.4f}")

# Print synergy values near the boundary for Red
print("\nRed Wing synergy from 100 to 130 days:")
for l in range(100, 131):
    print(f"  Lag {l}d: {syn_red[l-1]:.4f}")
