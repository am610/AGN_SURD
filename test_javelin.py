import sys
import numpy as np
import pandas as pd

# Configure path to JAVELIN
import javelin
from javelin.zylc import LightCurve
from javelin.lcmodel import Cont_Model, Rmap_Model

print("Loading test data...")
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

# Load and process Continuum
df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd_2440000', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd_2440000']
df_cont = df_cont[(df_cont['mjd'] >= 47512.0) & (df_cont['mjd'] <= 49255.0)].dropna()

# Load and process Spectroscopic Bins
df_hb = pd.read_csv(hb_bins_path)
df_hb = df_hb[(df_hb['mjd'] >= 47512.0) & (df_hb['mjd'] <= 49255.0)].dropna()

# Create zylclist for JAVELIN
cont_data = [df_cont['mjd'].values, df_cont['flux'].values, df_cont['err'].values]
core_data = [df_hb['mjd'].values, df_hb['core_flux'].values, df_hb['core_error'].values]

print("Initializing JAVELIN LightCurve...")
lc = LightCurve([cont_data, core_data], names=['continuum', 'core'])

print("Running JAVELIN Continuum Model...")
cont_lc = LightCurve([cont_data], names=['continuum'])
cont_model = Cont_Model(cont_lc)
# Small run for testing
cont_model.do_mcmc(nwalkers=20, nburn=10, nchain=20)
conthpd = cont_model.hpd
print("Continuum HPD:", conthpd)

print("Running JAVELIN Rmap Model...")
rmap = Rmap_Model(lc)
rmap.do_mcmc(conthpd=conthpd, laglimit=[[0.0, 120.0]], nwalkers=20, nburn=10, nchain=20)
print("Rmap HPD:", rmap.hpd)
print("JAVELIN ran successfully!")
