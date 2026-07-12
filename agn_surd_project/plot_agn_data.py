import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_ngc5548_agnwatch():
    df_cont = pd.read_csv("/home/ubuntu/agn_data/ngc5548_agnwatch_cont_clean.csv")
    df_hb = pd.read_csv("/home/ubuntu/agn_data/ngc5548_agnwatch_hb_clean.csv")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    ax1.errorbar(df_cont['mjd'], df_cont['flux'], yerr=df_cont['error'], fmt='o', markersize=2, label='Continuum 5100A')
    ax1.set_ylabel('Flux (10^-15)')
    ax1.legend()
    ax1.set_title('NGC 5548 AGN Watch')
    
    ax2.errorbar(df_hb['mjd'], df_hb['flux'], yerr=df_hb['error'], fmt='o', color='orange', markersize=2, label='H-beta')
    ax2.set_ylabel('Flux (10^-13)')
    ax2.set_xlabel('MJD')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig("/home/ubuntu/ngc5548_agnwatch_plot.png")
    plt.close()

def plot_mrk817_storm2():
    df_cont = pd.read_csv("/home/ubuntu/agn_data/mrk817_storm2_continuum-1180_clean.csv")
    df_civ = pd.read_csv("/home/ubuntu/agn_data/mrk817_storm2_civ-1590-1638_clean.csv")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    # STORM 2 FITS columns: hjd-2400000, flux, error
    ax1.errorbar(df_cont.iloc[:, 0], df_cont.iloc[:, 1], yerr=df_cont.iloc[:, 2], fmt='o', markersize=2, label='Continuum 1180A')
    ax1.set_ylabel('Flux (10^-13)')
    ax1.legend()
    ax1.set_title('Mrk 817 STORM 2')
    
    ax2.errorbar(df_civ.iloc[:, 0], df_civ.iloc[:, 1], yerr=df_civ.iloc[:, 2], fmt='o', color='red', markersize=2, label='C IV')
    ax2.set_ylabel('Flux (10^-13)')
    ax2.set_xlabel('HJD - 2400000')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig("/home/ubuntu/mrk817_storm2_plot.png")
    plt.close()

if __name__ == "__main__":
    print("Generating plots...")
    plot_ngc5548_agnwatch()
    plot_mrk817_storm2()
    print("Plots saved.")
