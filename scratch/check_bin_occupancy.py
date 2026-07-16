import sys
import numpy as np
import pandas as pd
from scipy.stats import zscore

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

# Let's check bin occupancy at lag = 73 days (the reported peak)
lag = 73
nbins = 6

future_target = core_z[lag:]
pred_1 = cont_z[:-lag]
pred_2 = blue_z[:-lag]
hist_var = core_z[:-lag]

data = np.vstack([future_target, pred_1, pred_2, hist_var]).T
n_samples = len(data)

# 4D histogram
hist_4d, _ = np.histogramdd(data, bins=nbins)

total_bins = nbins ** 4
empty_bins = np.sum(hist_4d == 0)
non_empty_bins = total_bins - empty_bins
fraction_empty = empty_bins / total_bins
occupancy_per_bin = n_samples / total_bins
actual_occupancy = hist_4d[hist_4d > 0]

print(f"Total samples: {n_samples}")
print(f"Total bins in 4D histogram (nbins={nbins}): {total_bins}")
print(f"Empty bins: {empty_bins} ({fraction_empty * 100:.2f}%)")
print(f"Non-empty bins: {non_empty_bins}")
print(f"Average samples per bin (overall): {occupancy_per_bin:.4f}")
print(f"Average samples per non-empty bin: {np.mean(actual_occupancy):.4f}")
print(f"Max samples in a single bin: {np.max(actual_occupancy)}")
print(f"Number of bins with occupancy of exactly 1 sample: {np.sum(actual_occupancy == 1)}")
