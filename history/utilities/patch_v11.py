import nbformat
import json

filename = '/Users/ayan/Programs/SURD/3_SURDS_AGN_V10_Restricted_MJD.ipynb'
out_filename = '/Users/ayan/Programs/SURD/3_SURDS_AGN_V11_Strict_Overlap.ipynb'

with open(filename, 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

for i, cell in enumerate(nb.cells):
    if cell.cell_type == 'code':
        # 1. Patch the loading cell to fix MJD
        if 'df_ngc5548_cont = pd.read_csv(ngc5548_cont_path)' in cell.source:
            # We want to replace the mjd column assignment
            code = cell.source
            code = code.replace("print(f\"Successfully loaded continuum data from: {ngc5548_cont_path}\")",
                                "df_ngc5548_cont['mjd'] = df_ngc5548_cont['jd_2440000']\n    df_ngc5548_hb['mjd'] = df_ngc5548_hb['jd_2440000']\n    print(f\"Successfully loaded continuum data from: {ngc5548_cont_path}\")")
            cell.source = code

        # 2. Patch the global interpolation logic
        if 'df_combined = pd.merge(df_cont_aligned, df_hb_velocity_aligned' in cell.source and 'min_overall_time' in cell.source:
            code = cell.source
            
            old_time_block = """
min_overall_time = df_combined['time'].min()
max_overall_time = df_combined['time'].max()

# Create a uniform time grid for resampling
uniform_time_grid = np.arange(min_overall_time, max_overall_time + dt_final, dt_final)
"""
            new_time_block = """
continuum_mjd_min = df_cont_aligned['time'].min()
continuum_mjd_max = df_cont_aligned['time'].max()
velocity_hbeta_mjd_min = df_hb_velocity_aligned['time'].min()
velocity_hbeta_mjd_max = df_hb_velocity_aligned['time'].max()

tmin = max(continuum_mjd_min, velocity_hbeta_mjd_min)
tmax = min(continuum_mjd_max, velocity_hbeta_mjd_max)

min_overall_time = tmin
max_overall_time = tmax

uniform_time_grid = np.arange(min_overall_time, max_overall_time + dt_final, dt_final)
"""
            code = code.replace(old_time_block.strip(), new_time_block.strip())

            old_interp_block = """
for col in columns_to_interpolate:
    # Forward fill then back fill to handle NaNs at the start/end of individual series before interpolation
    # Ensure the series is sorted by index (time) before ffill/bfill
    filled_series = df_combined.set_index('time')[col].sort_index().fillna(method='ffill').fillna(method='bfill')
    # Only interpolate if the series is not all NaN after filling
    if not filled_series.isnull().all():
        prepared_data[col] = np.interp(uniform_time_grid, filled_series.index, filled_series.values)
    else:
        prepared_data[col] = np.nan
"""
            new_interp_block = """
for col in columns_to_interpolate:
    valid_data = df_combined.set_index('time')[col].dropna().sort_index()
    if not valid_data.empty:
        prepared_data[col] = np.interp(uniform_time_grid, valid_data.index, valid_data.values, left=np.nan, right=np.nan)
    else:
        prepared_data[col] = np.nan

prepared_data = prepared_data.dropna().reset_index(drop=True)

assert prepared_data['time'].min() >= 47512, f"Min time {prepared_data['time'].min()} is below 47512"
assert prepared_data['time'].max() <= 49255, f"Max time {prepared_data['time'].max()} is above 49255"
"""
            code = code.replace(old_interp_block.strip(), new_interp_block.strip())
            cell.source = code

        # 3. Patch the prepare_data_for_cadence function
        if 'def prepare_data_for_cadence' in cell.source:
            code = cell.source
            
            old_time_block2 = """
    min_overall_time = df_combined['time'].min()
    max_overall_time = df_combined['time'].max()
    uniform_time_grid = np.arange(min_overall_time, max_overall_time + dt_cadence, dt_cadence)
"""
            new_time_block2 = """
    continuum_mjd_min = df_cont_aligned['time'].min()
    continuum_mjd_max = df_cont_aligned['time'].max()
    velocity_hbeta_mjd_min = df_hb_velocity_aligned['time'].min()
    velocity_hbeta_mjd_max = df_hb_velocity_aligned['time'].max()
    tmin = max(continuum_mjd_min, velocity_hbeta_mjd_min)
    tmax = min(continuum_mjd_max, velocity_hbeta_mjd_max)
    uniform_time_grid = np.arange(tmin, tmax + dt_cadence, dt_cadence)
"""
            code = code.replace(old_time_block2.strip(), new_time_block2.strip())
            
            old_interp_block2 = """
    for col in columns_to_interpolate:
        filled_series = df_combined.set_index('time')[col].sort_index().ffill().bfill()
        if not filled_series.isnull().all():
            prepared_data_cadence[col] = np.interp(uniform_time_grid, filled_series.index, filled_series.values)
        else:
            prepared_data_cadence[col] = np.nan
"""
            new_interp_block2 = """
    for col in columns_to_interpolate:
        valid_data = df_combined.set_index('time')[col].dropna().sort_index()
        if not valid_data.empty:
            prepared_data_cadence[col] = np.interp(uniform_time_grid, valid_data.index, valid_data.values, left=np.nan, right=np.nan)
        else:
            prepared_data_cadence[col] = np.nan
            
    prepared_data_cadence = prepared_data_cadence.dropna().reset_index(drop=True)
"""
            code = code.replace(old_interp_block2.strip(), new_interp_block2.strip())
            cell.source = code

with open(out_filename, 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)

print("V11 notebook created.")
