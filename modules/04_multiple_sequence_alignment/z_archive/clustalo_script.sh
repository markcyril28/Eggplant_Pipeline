#!/bin/bash
set -euo pipefail
# NOTE: This is a legacy script. The active pipeline uses d_msa_alignment.sh + modules/04_multiple_sequence_alignment/align.sh

# ===================== IMPORTANT VARIABLES =====================
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REFSEQS_DIR="$PIPELINE_DIR/1_RefSeqs/a_Smel_RefSeqs"

Eggplant_V3_Chromosomes_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3/Eggplant_V3_Chromosomes.fa"
Eggplant_V3_CDS_putative_function_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3/Eggplant_V3.CDS.putative_function.fa"
assembly_chr_FASTA_DIR="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3"
# ===============================================================


# Databases (Date Accessed: September 22, 2024)
Eggplant_V3_Chromosomes_DB="$assembly_chr_FASTA_DIR/Eggplant_V3_Chromosomes_db/Eggplant_V3_Chromosomes"
Eggplant_V3_CDS_putative_function_DB="$assembly_chr_FASTA_DIR/Eggplant_V3_CDS_putative_function_db/Eggplant_V3_CDS_putative_function"

# DMPs from Arabidopsis thaliana, Oryza sativa, and Zea mays (Date Accessed: September 22-26, 2024)

# AtDMPs from (https://www.arabidopsis.org/).
AtDMPs_v1_fasta_FILE="$PIPELINE_DIR/2_INPUTS/DMP/query_fasta/Arabidopsis_thaliana/v1_keyword_searched/AtDMPs_nuc_seq_merged_v1.fasta"
AtDMPs_fasta_FILE="$PIPELINE_DIR/2_INPUTS/DMP/query_fasta/Arabidopsis_thaliana/v2_based_from_paper/AtDMPs_nuc_seq_merged_v2.fasta"

# OsDMPs from the query DMP sequences use in the paper of <paper SWEETPOTATO>, and sequences lifted from annotated sequences obtained from (https://rice.uga.edu/).
OsDMPs_fasta_FILE="$PIPELINE_DIR/2_INPUTS/DMP/query_fasta/Oryza_sativa/OsDMPs_nuc_seq_merged.fasta"

# ZmDMPs from a simple keyword search of "DMP" in this website, (https://www.maizegdb.org/gene_center/gene).
ZmDMPs_v1_fasta_FILE="$PIPELINE_DIR/2_INPUTS/DMP/query_fasta/Zea_mays/v1/ZmDMPs_nuc_seq_merged_v1.fasta"
ZmDMPs_v2_fasta_FILE="$PIPELINE_DIR/2_INPUTS/DMP/query_fasta/Zea_mays/v2/ZmDMPs_nuc_seq_merged_v2.fasta"

# List and Parameters
Eggplant_db_LIST=($Eggplant_V3_Chromosomes_DB $Eggplant_V3_CDS_putative_function_DB)
DMP_fasta_LIST=($OsDMPs_fasta_FILE $AtDMPs_fasta_FILE $ZmDMPs_v1_fasta_FILE $ZmDMPs_v2_fasta_FILE)

clustalo \
	-i $AtDMPs_v1_fasta_FILE -o "04_MSA/AtDMPs_v1_clustalo.fa" -v --force
