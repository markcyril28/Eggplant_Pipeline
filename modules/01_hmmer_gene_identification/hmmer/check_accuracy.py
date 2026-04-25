#!/usr/bin/env python3
import sys
import pandas as pd

csv_path = sys.argv[1] if len(sys.argv) > 1 else 'Results.csv'
try:
    df = pd.read_csv(csv_path)
except FileNotFoundError:
    sys.exit(f"Error: file not found: {csv_path}")

if 'accuracy' not in df.columns:
    sys.exit(f"Error: 'accuracy' column not found in {csv_path}. Available columns: {list(df.columns)}")

print('Accuracy column statistics:')
print(f'Total rows: {len(df)}')
print(f'Non-empty accuracy values: {df["accuracy"].notna().sum()}')
print(f'Empty accuracy values: {df["accuracy"].isna().sum()}')
print(f'Accuracy range: {df["accuracy"].min()} to {df["accuracy"].max()}')
print(f'Average accuracy: {df["accuracy"].mean():.3f}')

print('\nSample accuracy values:')
print(df[['gene_id', 'domain_name', 'accuracy']].head(10))