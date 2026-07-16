import pandas as pd

df = pd.read_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/unconditioned_synergy_200d.csv')
print("Core peak:")
idx_core = df['core_syn'].idxmax()
print(df.loc[idx_core])

print("\nRed peak:")
idx_red = df['red_syn'].idxmax()
print(df.loc[idx_red])

print("\nBlue peak:")
idx_blue = df['blue_syn'].idxmax()
print(df.loc[idx_blue])
