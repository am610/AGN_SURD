import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configure path to JAVELIN
import javelin
from javelin.zylc import LightCurve
from javelin.lcmodel import Cont_Model, Rmap_Model

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

print("Loading raw data...")
cont_path = "/Users/ayan/Programs/SURD/agn_surd_project/agn_data/ngc5548_agnwatch/c5100.dat"
hb_bins_path = "/Users/ayan/Programs/SURD/agn_surd_project/processed/ngc5548_hb_velocity_bins.csv"

# Load Continuum (MJD 47512 to 49255)
df_cont = pd.read_csv(cont_path, sep=r'\s+', header=None, names=['jd_2440000', 'flux', 'err'])
df_cont['mjd'] = df_cont['jd_2440000']
df_cont = df_cont[(df_cont['mjd'] >= 47512.0) & (df_cont['mjd'] <= 49255.0)].dropna()

# Load Spectroscopic Bins (MJD 47512 to 49255)
df_hb = pd.read_csv(hb_bins_path)
df_hb = df_hb[(df_hb['mjd'] >= 47512.0) & (df_hb['mjd'] <= 49255.0)].dropna()

cont_data = [df_cont['mjd'].values, df_cont['flux'].values, df_cont['err'].values]

# ----------------- 1. FIT CONTINUUM DRW -----------------
print("Fitting JAVELIN Continuum Model...")
cont_lc = LightCurve([cont_data], names=['continuum'])
cont_model = Cont_Model(cont_lc)
cont_model.do_mcmc(nwalkers=80, nburn=150, nchain=300)
conthpd = cont_model.hpd
print("Continuum HPD:", conthpd)

# ----------------- 2. FIT VELOCITY WINGS -----------------
nwalkers = 80
nburn = 150
nchain = 300
laglimit = [[0.0, 120.0]]

targets = [
    ('Blue Wing', [df_hb['mjd'].values, df_hb['blue_wing_flux'].values, df_hb['blue_wing_error'].values]),
    ('Core H\\beta', [df_hb['mjd'].values, df_hb['core_flux'].values, df_hb['core_error'].values]),
    ('Red Wing', [df_hb['mjd'].values, df_hb['red_wing_flux'].values, df_hb['red_wing_error'].values])
]

results = {}

for name, t_data in targets:
    print(f"Fitting JAVELIN lag for {name}...")
    lc = LightCurve([cont_data, t_data], names=['continuum', name])
    rmap = Rmap_Model(lc)
    rmap.do_mcmc(conthpd=conthpd, laglimit=laglimit, nwalkers=nwalkers, nburn=nburn, nchain=nchain)
    
    # Save chain results
    # rmap.flatchain shape is (nwalkers * nchain, n_parameters)
    # The parameters are: log(sigma), log(tau), lag, width, scale
    chain = rmap.flatchain
    results[name] = {
        'lag_chain': chain[:, 2],
        'width_chain': chain[:, 3],
        'scale_chain': chain[:, 4]
    }
    
    # Calculate HPD
    hpd = rmap.hpd
    print(f"{name} HPD (lag, width, scale):", hpd)

# ----------------- 3. PLOT AND SAVE POSTERIORS -----------------
print("Generating Figure 6: JAVELIN Lag Posterior Distributions...")
fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

colors = ['#2ca02c', '#d62728', '#ff7f0e']

for idx, (name, res) in enumerate(results.items()):
    ax = axs[idx]
    lag_chain = res['lag_chain']
    
    # Calculate medians and 1-sigma bounds
    median = np.median(lag_chain)
    low = np.percentile(lag_chain, 16)
    high = np.percentile(lag_chain, 84)
    
    # Plot histogram
    ax.hist(lag_chain, bins=60, range=(0, 120), color=colors[idx], alpha=0.7, density=True, label=f'{name} Posterior')
    ax.axvline(median, color='black', linestyle='-', linewidth=2, label=f'Median: {median:.1f} d')
    ax.axvline(low, color='black', linestyle='--', linewidth=1.5, label=f'1$\\sigma$ range: [{low:.1f}, {high:.1f}] d')
    ax.axvline(high, color='black', linestyle='--', linewidth=1.5)
    
    ax.set_ylabel('Probability Density')
    ax.set_title(f'JAVELIN Lag Posterior: {name}')
    ax.legend(loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.5)

axs[2].set_xlabel('Lag (days)')
plt.tight_layout()
fig.savefig('overleaf_draft/figure6_javelin_posteriors.png', dpi=300)
plt.close(fig)

print("Figure 6: figure6_javelin_posteriors.png successfully created and saved in overleaf_draft/!")

# Save MCMC chains to CSV files for verification
output_dir = "/Users/ayan/Programs/SURD/agn_surd_project/processed/"
os.makedirs(output_dir, exist_ok=True)
for name, res in results.items():
    clean_name = name.lower().replace(" ", "_").replace("\\", "")
    pd.DataFrame({
        'lag': res['lag_chain'],
        'width': res['width_chain'],
        'scale': res['scale_chain']
    }).to_csv(os.path.join(output_dir, f"javelin_chain_{clean_name}.csv"), index=False)
print("Saved all JAVELIN MCMC flat chains to CSV files.")
