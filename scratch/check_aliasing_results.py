import pandas as pd
df = pd.read_csv('/Users/ayan/Programs/SURD/agn_surd_project/processed/seasonal_aliasing_null_test.csv')

print("Seasonal aliasing test results (top 10 highest synergy lags):")
print(df.sort_values('median_syn', ascending=False).head(10))

print("\nSynergy values at lags 15, 75, 110, 160:")
for l in [15, 75, 110, 160]:
    row = df[df['lag'] == l]
    if not row.empty:
        print(f"  Lag {l}d: Median={row['median_syn'].values[0]:.4f}, 95%={row['p97_5_syn'].values[0]:.4f}")
