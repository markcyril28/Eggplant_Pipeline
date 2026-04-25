#!/usr/bin/env python3
import pandas as pd

# Read the results
df = pd.read_csv('Results.csv')

print('=== HMMER Domain Analysis Results ===')
print(f'Total domain hits: {len(df)}')
print(f'Unique genes identified: {df["gene_id"].nunique()}')

print('\nDomain distribution:')
domain_stats = df.groupby('domain_name').agg({
    'gene_id': 'nunique', 
    'domain_name': 'count'
}).rename(columns={'gene_id': 'unique_genes', 'domain_name': 'total_hits'})
print(domain_stats)

print('\nTop genes by domain count:')
gene_stats = df.groupby('gene_id').agg({
    'domain_name': lambda x: len(set(x)), 
    'gene_id': 'count'
}).rename(columns={'domain_name': 'unique_domains', 'gene_id': 'total_hits'}).sort_values('total_hits', ascending=False)
print(gene_stats.head(10))

print('\nGenes with multiple domain types:')
multi_domain_genes = gene_stats[gene_stats['unique_domains'] > 1].sort_values('unique_domains', ascending=False)
if len(multi_domain_genes) > 0:
    print(multi_domain_genes)
    # Show which domains each multi-domain gene has
    for gene_id in multi_domain_genes.index:
        domains = df[df['gene_id'] == gene_id]['domain_name'].unique()
        print(f"  {gene_id}: {', '.join(domains)}")
else:
    print("No genes found with multiple domain types")

print('\nDomain type summary:')
for domain in df['domain_name'].unique():
    domain_df = df[df['domain_name'] == domain]
    print(f"\n{domain} domain:")
    print(f"  - Total hits: {len(domain_df)}")
    print(f"  - Unique genes: {domain_df['gene_id'].nunique()}")
    print(f"  - Best E-value: {domain_df['full_seq_evalue'].min()}")
    print(f"  - Genes: {', '.join(domain_df['gene_id'].unique()[:5])}{'...' if domain_df['gene_id'].nunique() > 5 else ''}")