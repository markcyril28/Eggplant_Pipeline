#!/bin/bash
set -euo pipefail

# List of RefSeq IDs
NCBI_RefSeq_IDS_LIST=(
    "XM_008677403.2"
    "NM_001371746.1"
    "NM_001150103.2"
    "CM000784.4"
    "XM_008677254.4"
    "NM_001150831.1"
    "XM_008660852.4"
    "NM_001151895.2"
    "XM_020544574.3"
    "XM_008659040.1"
    "XM_008658511.4"
    "XM_008657814.4"
    "XM_008658777.3"
    "NM_001151404.2"
)

gene_name_LIST=(
	"probable_WRKY_transcription_factor_63" 
	"ZmDMP1" 
	"ZmDMP2" 
	"ZmDMP2" 
	"ZmDMP3_7_weird" 
	"ZmDMP3_7_weird" 
	"ZmDMP3_7_weird" 
	"ZmDMP4" 
	"ZmDMP4" 
	"ZmDMP4_7_weird" 
	"ZmDMP6" 
	"ZmDMP6_7_weird" 
	"ZmDMP7" 
	"ZmDMP8_7_weird"
)


# Get the length of the arrays (assuming they are the same length)
length=${#NCBI_RefSeq_IDS_LIST[@]}

# Loop through the indexes
for ((i=0; i<length; i++)); do
    refseq_id="${NCBI_RefSeq_IDS_LIST[i]}"
    gene_name="${gene_name_LIST[i]}"
    echo "Processing RefSeq ID: $refseq_id with Gene Name: $gene_name"
    wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=$refseq_id&db=nuccore&report=fasta&retmode=text" -O "${gene_name}_$refseq_id.fasta"
done




: << 'SAMPLE'

# Loop through each RefSeq ID
for ID in "${REFSEQ_IDS[@]}"; do
    echo "Downloading FASTA for $ID..."
    wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=$ID&db=nuccore&report=fasta&retmode=text" -O "${ID}.fasta"
done

# Set the Internal Field Separator to a comma
IFS=,

# Read the CSV file line by line
while read column1 column2; do
    # Skip the header
    if [[ "$column1" != "column_1" ]]; then
        echo "Column 1: $column1, Column 2: $column2"
        #echo "Downloading FASTA for $column1..."
        #wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=$column1&db=nuccore&report=fasta&retmode=text" -O "${column1}_${column2}.fasta"
    fi
done < List_of_ZmDMPs.csv

SAMPLE