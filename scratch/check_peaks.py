import pandas as pd
df = pd.read_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/unconditioned_synergy_200d.csv')

print("Core Target top 5 synergy values:")
print(df.sort_values('core_synergy', ascending=False).head(5))

print("\nRed Target top 5 synergy values:")
print(df.sort_values('red_synergy', ascending=False).head(5))

print("\nBlue Target top 5 synergy values:")
print(df.sort_values('blue_synergy', ascending=False).head(5))
