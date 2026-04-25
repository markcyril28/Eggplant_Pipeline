#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv('Results.csv')

print('=== Verification of sorting by domain_name ===')
print('\nFirst 5 rows:')
print(df[['gene_id', 'domain_name']].head())

print('\nDomain name groups:')
for domain in df['domain_name'].unique():
    count = len(df[df['domain_name'] == domain])
    print(f'{domain}: {count} hits')

print('\nRow ranges for each domain:')
for domain in df['domain_name'].unique():
    indices = df[df['domain_name'] == domain].index.tolist()
    print(f'{domain}: rows {min(indices)+2}-{max(indices)+2} (including header)')

print('\nSample from each domain group:')
for domain in df['domain_name'].unique():
    domain_df = df[df['domain_name'] == domain]
    print(f'\n{domain} domain (first 3 genes):')
    print(domain_df[['gene_id', 'domain_name', 'full_seq_evalue']].head(3).to_string(index=False))