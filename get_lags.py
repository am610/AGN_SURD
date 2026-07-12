import pandas as pd
import glob

print("=== V11 STRICT OVERLAP LAGS ===")
for f in glob.glob("/Users/ayan/Programs/SURD/agn_surd_project/plots/test_c_shuffle_results/robustness_surrogate_*.csv"):
    df = pd.read_csv(f)
    print(f"\n{f.split('/')[-1]}")
    # Get the lag with the maximum synergy
    max_idx = df['real_synergy'].idxmax()
    print(f"Max Synergy Lag: {df.loc[max_idx, 'lag']}")
