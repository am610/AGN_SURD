import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

def load_agn_data(csv_path):
    """Loads a cleaned AGN light curve CSV."""
    return pd.read_csv(csv_path)

def interpolate_light_curve(df, time_col, flux_col, new_times, kind='linear'):
    """Interpolates a light curve onto a new time grid."""
    f = interp1d(df[time_col], df[flux_col], kind=kind, fill_value="extrapolate")
    return f(new_times)

def standardize(data):
    """Standardizes data (zero mean, unit variance)."""
    return (data - np.mean(data)) / np.std(data)

def prepare_surd_input(cont_df, line_df, time_range=None, dt=1.0):
    """
    Prepares continuum and line data for SURD analysis.
    Resamples both to a common grid.
    """
    if time_range is None:
        t_min = max(cont_df.iloc[:, 0].min(), line_df.iloc[:, 0].min())
        t_max = min(cont_df.iloc[:, 0].max(), line_df.iloc[:, 0].max())
    else:
        t_min, t_max = time_range
        
    common_times = np.arange(t_min, t_max, dt)
    
    # Assuming first column is time, second is flux
    c_flux = interpolate_light_curve(cont_df, cont_df.columns[0], cont_df.columns[1], common_times)
    l_flux = interpolate_light_curve(line_df, line_df.columns[0], line_df.columns[1], common_times)
    
    # Standardize
    c_std = standardize(c_flux)
    l_std = standardize(l_flux)
    
    return common_times, c_std, l_std

def stack_for_surd(time_series_list):
    """Stacks multiple standardized time series into a single array X."""
    return np.vstack(time_series_list)
