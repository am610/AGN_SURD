import sys
import numpy as np
import pandas as pd
from scipy.stats import zscore

# Configure path to SURD utilities
sys.path.append('/Users/ayan/Programs/SURD/SURD')
sys.path.append('/Users/ayan/Programs/SURD/SURD/utils')
from utils import surd

# Load aligned data
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd', 'flux', 'err'])
df_hb = pd.read_csv(hb_bins_path)

tmin = max(df_cont['jd'].min(), df_hb['mjd'].min())
tmax = min(df_cont['jd'].max(), df_hb['mjd'].max())

uniform_time_grid = np.arange(tmin, tmax + 1.0, 1.0)
prepared_data = pd.DataFrame({'time': uniform_time_grid})

valid_cont = df_cont.dropna(subset=['flux']).sort_values('jd')
prepared_data['cont_flux_zscore'] = zscore(np.interp(uniform_time_grid, valid_cont['jd'], valid_cont['flux']))

for col, new_name in [('blue_wing_flux', 'blue_wing_flux_zscore'), 
                      ('core_flux', 'core_flux_zscore'), 
                      ('red_wing_flux', 'red_wing_flux_zscore')]:
    valid_data = df_hb.dropna(subset=[col]).sort_values('mjd')
    prepared_data[new_name] = zscore(np.interp(uniform_time_grid, valid_data['mjd'], valid_data[col]))

prepared_data = prepared_data.dropna().reset_index(drop=True)

cont = prepared_data['cont_flux_zscore'].values
blue = prepared_data['blue_wing_flux_zscore'].values
core = prepared_data['core_flux_zscore'].values

# Let's run a test decomposition for lag = 10 days
lag = 10
nbins = 8

# target core at t+lag, predictors: cont, blue, core at t
X = np.vstack([cont, blue, core])
Y_data = np.vstack([core[lag:], X[:, :-lag]])
hist, _ = np.histogramdd(Y_data.T, nbins)

I_R, I_S, MI, info_leak = surd.surd(hist)

print("--- SURD Decomposition Results ---")
print("Target: Core(t+10)")
print("Predictors: Continuum(t), Blue(t), Core(t)")
print(f"Total target entropy H(Y) = {np.sum(-hist.sum(axis=(1,2,3)) * np.log2(hist.sum(axis=(1,2,3)) + 1e-14)):.4f}")
print("Mutual Information MI (raw values):")
for k, v in MI.items():
    print(f"  MI for comb {k} = {v:.4f}")

print("\nRedundancies I_R (keys are predictor combinations):")
sum_r = 0.0
for k, v in I_R.items():
    print(f"  I_R for {k} = {v:.4f}")
    sum_r += v

print("\nSynergies I_S (keys are predictor combinations):")
sum_s = 0.0
for k, v in I_S.items():
    print(f"  I_S for {k} = {v:.4f}")
    sum_s += v

print(f"\nSum of Redundancies (including uniques) = {sum_r:.4f}")
print(f"Sum of Synergies = {sum_s:.4f}")
total_decomposed = sum_r + sum_s
print(f"Total Decomposed Information (Sum(I_R) + Sum(I_S)) = {total_decomposed:.4f}")

# Let's check joint mutual information of the target with all 3 predictors
# Joint entropy H(Y, X1, X2, X3) is entropy of hist
H_joint = -np.sum(hist * np.log2(hist + 1e-14))
# Predictors marginal entropy H(X1, X2, X3)
p_pred = hist.sum(axis=0)
H_pred = -np.sum(p_pred * np.log2(p_pred + 1e-14))
# Target marginal entropy H(Y)
p_target = hist.sum(axis=(1,2,3))
H_target = -np.sum(p_target * np.log2(p_target + 1e-14))

# Joint Mutual Information
MI_joint = H_target + H_pred - H_joint
print(f"Joint Mutual Information I(Y; X1, X2, X3) = {MI_joint:.4f}")

print(f"Difference = Joint MI - Total Decomposed = {MI_joint - total_decomposed:.4f}")
print(f"Leak H(Y|X1, X2, X3) = {H_joint - H_pred:.4f}")
print(f"info_leak returned by surd() = {info_leak:.4f}")
print(f"info_leak calculated as H(Y|X)/H(Y) = {(H_joint - H_pred)/H_target:.4f}")
