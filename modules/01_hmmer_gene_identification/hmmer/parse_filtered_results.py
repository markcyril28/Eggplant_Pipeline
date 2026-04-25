#!/usr/bin/env python3
"""
Script to parse HMMER filtered domain results and create a consolidated CSV file
"""

import argparse
import csv
import re
import sys
from pathlib import Path

def parse_domtbl_line(line):
    """Parse a single line from a HMMER domain table"""
    # Split the line and handle the description field which might contain spaces
    parts = line.strip().split()
    
    # The description starts after the accuracy score (22nd field)
    main_fields = parts[:22]
    description = ' '.join(parts[22:]) if len(parts) > 22 else ''
    
    # Extract gene name and annotation from description
    gene_match = re.search(r'gene=([^\s]+)', description)
    gene_name = gene_match.group(1) if gene_match else ''
    
    name_match = re.search(r'Name:"([^"]*)"', description)
    annotation = name_match.group(1) if name_match else ''
    
    return {
        'gene_id': main_fields[0],
        'protein_length': main_fields[2],  # tlen (target length)
        'domain_name': main_fields[3],     # query name (SSXT, WRC, QLQ)
        'pfam_id': main_fields[4],         # query accession (PF05030.17, etc.)
        'model_length': main_fields[5],    # qlen (query length)
        'full_seq_evalue': main_fields[6],
        'full_seq_score': main_fields[7],
        'full_seq_bias': main_fields[8],
        'domain_num': main_fields[9],
        'total_domains': main_fields[10],
        'c_evalue': main_fields[11],
        'i_evalue': main_fields[12],
        'domain_score': main_fields[13],
        'domain_bias': main_fields[14],
        'hmm_from': main_fields[15],
        'hmm_to': main_fields[16],
        'ali_from': main_fields[17],
        'ali_to': main_fields[18],
        'env_from': main_fields[19],
        'env_to': main_fields[20],
        'accuracy': main_fields[21] if len(main_fields) > 21 else '',
        'gene_name': gene_name,
        'annotation': annotation
    }

def process_filtered_files(base_dir: Path, output_file: Path):
    """Process all filtered domain table files and create consolidated CSV"""
    
    # Domain files to process
    domain_files = {
        'PF05030_SSXT': base_dir / 'PF05030_SSXT' / 'PF05030_SSXT_hits_filtered.domtbl',
        'PF08879_WRC': base_dir / 'PF08879_WRC' / 'PF08879_WRC_hits_filtered.domtbl',
        'PF08880_QLQ': base_dir / 'PF08880_QLQ' / 'PF08880_QLQ_hits_filtered.domtbl'
    }
    
    all_results = []
    
    # Process each domain file
    for domain_type, file_path in domain_files.items():
        if file_path.exists():
            print(f"Processing {domain_type}...")
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):  # Skip empty lines and comments
                        try:
                            result = parse_domtbl_line(line)
                            all_results.append(result)
                        except Exception as e:
                            print(f"Error parsing line in {domain_type}: {line}")
                            print(f"Error: {e}")
        else:
            print(f"Warning: {file_path} not found")
    
    # Sort results by gene_id first, then by full_seq_evalue (best E-values first)
    all_results.sort(key=lambda x: (x['gene_id'], float(x['full_seq_evalue'])))
    
    # Write to CSV
    if all_results:
        fieldnames = [
            'gene_id', 'gene_name', 'annotation', 'protein_length',
            'domain_name', 'pfam_id', 'model_length',
            'full_seq_evalue', 'full_seq_score', 'full_seq_bias',
            'domain_num', 'total_domains',
            'c_evalue', 'i_evalue', 'domain_score', 'domain_bias',
            'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to',
            'accuracy'
        ]
        
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)
        
        print(f"\nSuccessfully wrote {len(all_results)} results to {output_file}")
        
        # Print summary statistics
        domain_counts = {}
        for result in all_results:
            domain = result['domain_name']
            domain_counts[domain] = domain_counts.get(domain, 0) + 1
        
        print("\nSummary by domain:")
        for domain, count in sorted(domain_counts.items()):
            print(f"  {domain}: {count} hits")
            
        # Count unique genes
        unique_genes = set(result['gene_id'] for result in all_results)
        print(f"\nTotal unique genes: {len(unique_genes)}")
        
    else:
        print("No results found to write to CSV")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse HMMER filtered domain results into a consolidated CSV"
    )
    parser.add_argument(
        "base_dir",
        nargs="?",
        default=".",
        help="Directory containing per-domain subdirs with *_hits_filtered.domtbl files (default: cwd)"
    )
    parser.add_argument(
        "-o", "--output",
        default="Results.csv",
        help="Output CSV path (default: Results.csv in cwd)"
    )
    args = parser.parse_args()
    process_filtered_files(Path(args.base_dir), Path(args.output))