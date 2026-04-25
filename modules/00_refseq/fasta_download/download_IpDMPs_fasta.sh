#!/bin/bash
set -euo pipefail

: << 'OLD_WAY'
# List of RefSeq IDs
NCBI_RefSeq_IDS_LIST=(
    "NM_001371746.1"
    "XM_020544574.3"
    "XM_008677403.2"
    "XM_008677254.4"
    "XM_008660852.4"
    "XM_008659040.1"
    "XM_008658777.3"
    "XM_008658511.4"
    "XM_008657814.4"
    "XM_008647463.2"
    "NM_001151895.2"
    "NM_001151404.2"
    "NM_001150831.1"
    "NM_001150103.2"
)

Name_in_Paper_LIST=(
    "ZmDMP1"
    "ZmDMP13"
    "ZmDMP3"
    "ZmDMP2"
    "ZmDMP11"
    "ZmDMP7"
    "ZmDMP10"
    "ZmDMP9"
    "ZmDMP8"
    "ZmDMP4"
    "ZmDMP14"
    "ZmDMP12"
    "ZmDMP5"
    "ZmDMP15"
)

GenPept_name_LIST=(
    "DMP1-like"
    "DMP4"
    "probable_WRKY_transcription_factor_63"
    "DMP3"
    "DMP3"
    "DMP4"
    "DMP7"
    "DMP6"
    "DMP6"
    "DMP6"
    "DMP4-like"
    "DMP8"
    "DMP3-like"
    "DMP2"
)

GenBank_name_LIST=(
    "uncharacterized_LOC100277531"
    "ZmDMP4_LOC103639910"
    "probable_WRKY_transcription_factor_63"
    "ZmDMP7_LOC103651591"
    "ZmDMP7_LOC103637805"
    "ZmDMP7_LOC103636687"
    "ZmDMP7_LOC103636425"
    "ZmDMP6"
    "ZmDMP7_LOC103635334"
    "ZmDMP7_LOC103627158"
    "ZmDMP4"
    "ZmDMP7_LOC100277972"
    "ZmDMP7_LOC100277191"
    "ZmDMP2"
)


# Get the length of the arrays (assuming they are the same length)
length=${#NCBI_RefSeq_IDS_LIST[@]}

# Loop through the indexes
for ((i=0; i<length; i++)); do
    echo $i
    refseq_id="${NCBI_RefSeq_IDS_LIST[i]}"
    ZmP_name="${Name_in_Paper_LIST[i]}"
    GP_name="${GenPept_name_LIST[i]}"
    GB_name="${GenBank_name_LIST[i]}"
    echo "Processing RefSeq ID: $refseq_id with Names: Paper Name: $ZmP_name, GenPept Name: $GP_name, GenBank: $GB_name "
    wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=$refseq_id&db=nuccore&report=fasta&retmode=text" -O "${refseq_id}_ZmP_${ZmP_name}_GP_${GP_name}_GB_${GB_name}.fasta"
done

OLD_WAY

# Set the field separator to comma for CSV file
IFS=','

# Loop through the CSV file line by line
while read -r column1 column2 column3; do
    # Do something with the data, e.g., print it
    echo "Column 1: $column1"
    echo "Column 2: $column2"
    echo "Column 3: $column3"
    echo "-------------"
done < input_file.csv
