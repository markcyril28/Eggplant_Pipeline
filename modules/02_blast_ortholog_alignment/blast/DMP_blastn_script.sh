#!/bin/bash
set -euo pipefail
# NOTE: This is a legacy script. The active pipeline uses 2_BLAST_Alignment.sh + modules/01_identification/blastn.sh

# ===================== IMPORTANT VARIABLES =====================
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REFSEQS_DIR="$PIPELINE_DIR/1_RefSeqs/a_Smel_RefSeqs"

Eggplant_V3_assembly_chr_FASTA_DIR="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3"
Eggplant_V3_Chromosomes_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3/Eggplant_V3_Chromosomes.fa"
Eggplant_V3_CDS_putative_function_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3/Eggplant_V3.CDS.putative_function.fa"

Eggplant_V4_assembly_chr_FASTA_DIR="$REFSEQS_DIR/Solanum_melongena_V4_pangenome/Eggplant_V4"
Eggplant_V4_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_V4_pangenome/Eggplant_V4/Eggplant_V4.fa"
Eggplant_V4_Chromosomes_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_V4_pangenome/Eggplant_V4/"
Eggplant_V4_chr0_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_V4_pangenome/Eggplant_V4/Eggplant_V4.chromosomes.fa"

Eggplant_V4_1_fasta_FASTA_DIR="$REFSEQS_DIR/Solanum_melongena_v4.1"
Eggplant_V4_1_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1.fa"
Eggplant_V4_1_transcript_fasta_FILE="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1_transcripts.function.fa"
# ===============================================================

#mkdir 01_Identification

# Installation of NCBI-BLAST
#sudo apt install ncbi-blast+ -y

: << 'MAKEBLASTDB'

mkdir -p "$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_Chromosomes_db"
makeblastdb -in $Eggplant_V3_Chromosomes_fasta_FILE \
	-dbtype nucl \
	-out "$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_Chromosomes_db/Eggplant_V3_Chromosomes"

mkdir -p 01_Identification/DMP_query_fasta_file

mkdir -p "$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_CDS_putative_function_db"
makeblastdb -in $Eggplant_V3_CDS_putative_function_fasta_FILE \
	-dbtype nucl \
	-out "$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_CDS_putative_function_db/Eggplant_V3_CDS_putative_function"

Eggplant_V4_assembly_chr_FASTA_DIR="$REFSEQS_DIR/Solanum_melongena_V4_pangenome/Eggplant_V4"
for fasta in $Eggplant_V4_assembly_chr_FASTA_DIR/*; do
	fasta_BASENAME=$(basename $fasta .fa) 
	mkdir -p "$Eggplant_V4_assembly_chr_FASTA_DIR/${fasta_BASENAME}_db"
	makeblastdb -in $fasta \
		-dbtype nucl \
		-out "$Eggplant_V4_assembly_chr_FASTA_DIR/${fasta_BASENAME}_db/$fasta_BASENAME"
done

for fasta in $Eggplant_V4_1_fasta_FASTA_DIR/*; do
	fasta_BASENAME=$(basename $fasta .fa) 
	mkdir -p "$Eggplant_V4_1_fasta_FASTA_DIR/${fasta_BASENAME}_db"
	makeblastdb -in $fasta \
		-dbtype nucl \
		-out "$Eggplant_V4_1_fasta_FASTA_DIR/${fasta_BASENAME}_db/$fasta_BASENAME"
done

for fasta in $Eggplant_V4_1_fasta_FASTA_DIR/*; do
	fasta_BASENAME=$(basename $fasta .fa) 
	mkdir -p "$Eggplant_V4_1_fasta_FASTA_DIR/${fasta_BASENAME}_db"
	makeblastdb -in $fasta \
		-dbtype nucl \
		-out "$Eggplant_V4_1_fasta_FASTA_DIR/${fasta_BASENAME}_db/$fasta_BASENAME"
done

MAKEBLASTDB


: << 'SINGLE_RUN'
blastn -db "$assembly_chr_FASTA_DIR/Eggplant_V3_Chromosomes_db/Eggplant_V3_Chromosomes" \
	-query $AtDMPs_fasta_FILE \
	-evalue 1e-10 -word_size 11 > 01_Identification/blastn_V3_Chromosomes_against_AtDMPs_result.txt

SINGLE_RUN


# -------------------------------------------------------------------------------------------------------------------------------
# Databases (Date Accessed: September 22, 2024)
Eggplant_V3_Chromosomes_DB="$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_Chromosomes_db/Eggplant_V3_Chromosomes"
Eggplant_V3_CDS_putative_function_DB="$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_CDS_putative_function_db/Eggplant_V3_CDS_putative_function"

Eggplant_V4_DB="$Eggplant_V4_assembly_chr_FASTA_DIR/Eggplant_V4_db/Eggplant_V4"
Eggplant_V4_Chromosomes_DB="$Eggplant_V4_assembly_chr_FASTA_DIR/Eggplant_V4.chromosomes_db/Eggplant_V4.chromosomes"
Eggplant_V4_chr0="$Eggplant_V4_assembly_chr_FASTA_DIR/Eggplant_V4_chr0_db/Eggplant_V4_chr0"

Eggplant_V4_1_DB="$Eggplant_V4_1_fasta_FASTA_DIR/Eggplant_V4.1_db/Eggplant_V4.1"
Eggplant_V4_1_transcript_DB="$Eggplant_V4_1_fasta_FASTA_DIR/Eggplant_V4.1_transcripts.function_db/Eggplant_V4.1_transcripts.function"

# -------------------------------------------------------------------------------------------------------------------------------
# DMPs from Arabidopsis thaliana, Oryza sativa, and Zea mays (Date Accessed: September 22-26, 2024)

# Keyword searched of AtDMPs from (https://www.arabidopsis.org/).
AtDMPs_v1_fasta_FILE="01_Identification/DMP_query_fasta_file/Arabidopsis_thaliana/v1_keyword_searched/AtDMPs_nuc_seq_merged_v1.fasta"
# List of AtDMPs based on SweetPotato_paper, and nucleotide fasta sequences lifted from blastn datasets of Aradipopsis thalian from (https://www.arabidopsis.org/download/list?dir=Sequences%2FTAIR10_blastsets)
AtDMPs_v2_fasta_FILE="01_Identification/DMP_query_fasta_file/Arabidopsis_thaliana/v2_based_from_paper/AtDMPs_nuc_seq_merged_v2.fasta"

# OsDMPs from the query DMP sequences use in the paper of <paper SWEETPOTATO>, and sequences lifted from annotated sequences obtained from (https://rice.uga.edu/).
OsDMPs_fasta_FILE="01_Identification/DMP_query_fasta_file/Oryza_sativa/OsDMPs_nuc_seq_merged.fasta"

# ZmDMPs from a simple keyword search of "DMP" in this website, (https://www.maizegdb.org/gene_center/gene).
ZmDMPs_v1_fasta_FILE="01_Identification/DMP_query_fasta_file/Zea_mays/v1/ZmDMPs_nuc_seq_merged_v1.fasta"
# ZmDMPs from SweetPotato_paper, and searched/blastp aligned in NCBI. 
ZmDMPs_v2_fasta_FILE="01_Identification/DMP_query_fasta_file/Zea_mays/v2_based_on_paper_blastp/ZmDMPs_nuc_seq_merged_v2.fa"
# ZmDMPs_v2 fastas aligned using cds assembly fasta
ZmDMPs_v3_fasta_FILE=""
# Zea_mays_Zm-B73-REFERENCE-NAM-5.0_ncbi_dataset
ZmDMPs_v4_fasta_FILE="01_Identification/DMP_query_fasta_file/Zea_mays/v4_based_on_genome/ZmDMPs_nuc_seq_merged_v4.fasta"

# Ipomoea batatas
IpDMPs_fasta_FILE="01_Identification/DMP_query_fasta_file/Ipomea_batatas/IpDMPs_merged_fasta.fa"

# Capsicum_annuum_UCD10Xv1.1_ncbi_dataset
CaDMPs_fasta_FILE="01_Identification/DMP_query_fasta_file/Capsicum_annuum_Pepper/CaDMPs_nuc_seq_merged.fasta"

# Glycine max (Wm82.a2.v1); used from this paper (10.3389/fpls.2023.1216082); lifted from (); interestingly mentioned in this paper (Mutation of GmDMP genes triggers haploid induction in soybean)
GmDMP_fasta_FILE_v1="01_Identification/DMP_query_fasta_file/Glycine_max_Pepper/GmDMPs_nuc_seq_merged_v1.fasta"


# Other_DMPs from SweetPotato_paper
OtherDMPs_fasta_FILE="01_Identification/DMP_query_fasta_file/Other_DMPs/OtherDMPs_merged_fasta.fa"
# -------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------
# Selected DMPs Lifted from the 2024 Paper. 
# GROUP_1
#Setaria italica. Genome Assembly Setaria_italica_v2.0. Lifted from (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000263155.2/).
SiDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_1_Setaria_italica/SiDMPs_merged_fasta.fa"
# Physcomitrella patens. Pddmp lifted from Phypa V3 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002425.4/).
Pddmp_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_1_Physcomitrella_patens/PpDMPs_merged_fasta.fa"

# GROUP_2
# Citrus sinensis. Genome Assembly DVS_A1.0 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_022201045.2/).
CsDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_2_Citrus_sinensis/CsDMPs_merged_fasta.fa"
# Ananas comosus. ASM154086v1 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001540865.1/).
AcDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_2_Ananas_comosus/AcDMPs_merged_fasta.fa"

# GROUP_3
# Musa acuminata. Genome Assembly Cavendish_Baxijiao_AAA (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036884655.1/).
MaDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_3_Musa_acuminata/MaDMPs_merged_fasta.fa"
# Glycine max. Genome Assembly Glycine_max_v4.0 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004515.6/). 
GmDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_3_Glycine_max/GmDMPs_merged_fasta.fa"

# GROUP_4
# Solanum lycopersicum. SlDMP lifted from Genome Assembly SLM_r2.1 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036512215.1/).
SlDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_4_Solanum_lycopersicum/SlDMPs_merged_fasta.fa"
# Populus trichocarpa. PtDMP lifted from (Genome Assembly P.trichocarpa_v4.1). (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002775.5/).
PtDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_4_Populus_trichocarpa/PtDMPs_merged_fasta.fa"

# GROUP_5
# Gossypium raimondii. Genome Assembly ASM2569854v1 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_025698545.1/).
GrDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_5_Gossypium_raimondii/GrDMPs_merged_fasta.fa"
# Brachypodium distachyon. Genome Assembly Brachypodium_distachyon_v3.0 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005505.3/).
BdDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_5_Brachypodium_distachyon/BdDMPs_merged_fasta.fa"
# Sorghum bicolor. Genome Assembly Sorghum_bicolor_NCBIv3 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003195.3/).
SbDMP_fasta_FILE="01_Identification/DMP_query_fasta_file/2024_Paper_DMPs/GROUP_5_Sorghum_bicolor/SbDMPs_merged_fasta.fa"

# -------------------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------------------------
# Lists and Parameters
Eggplant_db_LIST=(
	#"$Eggplant_V3_Chromosomes_DB"
	#"$Eggplant_V3_CDS_putative_function_DB"
	#"$Eggplant_V4_DB"
	#"$Eggplant_V4_Chromosomes_DB"
	"$Eggplant_V4_1_DB"
	#"$Eggplant_V4_1_transcript_DB"
)

DMP_fasta_LIST=(
	"$AtDMPs_v1_fasta_FILE"
	"$AtDMPs_v2_fasta_FILE"
	"$OsDMPs_fasta_FILE"
	#"$ZmDMPs_v1_fasta_FILE"
	#"$ZmDMPs_v2_fasta_FILE"
	"$ZmDMPs_v4_fasta_FILE"
	"$IpDMPs_fasta_FILE"
	"$CaDMPs_fasta_FILE"
	"$GmDMP_fasta_FILE_v1"
	"$OtherDMPs_fasta_FILE"
	"$SiDMP_fasta_FILE"
	"$Pddmp_fasta_FILE"
	"$CsDMP_fasta_FILE"
	#"$AcDMP_fasta_FILE"
	"$MaDMP_fasta_FILE"
	"$GmDMP_fasta_FILE"
	"$SlDMP_fasta_FILE"
	"$PtDMP_fasta_FILE"
	"$GrDMP_fasta_FILE"
	"$BdDMP_fasta_FILE"
	#"$SbDMP_fasta_FILE"
)

# -------------------------------------------------------------------------------------------------------------------------------

# NOTE: Lower e-value -> more significant match
# 0.05 	- Default in NCBI Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi); 
# 1e-10 - Default in-website blast in solgenomics (https://solgenomics.net/tools/blast/?db_id=320), and in paper<>.
#e_value_LIST=(0.05 1e-10)
e_value_LIST=(1e-10)

# NOTE: 
# - smaller word size, more sensitive, but increases false positive; 
# - longer word size, less sensitive, but require specificity (longer exact matches)
# 11 - Default in solgenomics (https://solgenomics.net/tools/blast/?db_id=320).
# 28 - Default in NCBI Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
#word_size_LIST=(11 20 28)
word_size_LIST=(11)


: << 'DEBUGGER'
for Eggplant_db in ${Eggplant_db_LIST[@]}; do
	echo $Eggplant_db
done

for DMP_fasta in ${DMP_fasta_LIST[@]}; do
	echo $DMP_fasta
done
DEBUGGER

#: << 'BLASTN_FOR_LOOP'
rm -rf 01_Identification/DMP_blastn_results/* 2>/dev/null || true

for DB in "${Eggplant_db_LIST[@]}"; do
	for DMP_fasta_FILE in "${DMP_fasta_LIST[@]}"; do
		for e_value in "${e_value_LIST[@]}"; do
			for word_size in "${word_size_LIST[@]}"; do
				DB_BASENAME=$(basename "$DB")
				DMP_fasta_BASENAME=$(basename "$DMP_fasta_FILE" .fasta)

				#mkdir -p "01_Identification/DMP_blastn_results/ev_${e_value}_ws_${word_size}"
				mkdir -p "01_Identification/DMP_blastn_results/ev_${e_value}_ws_${word_size}_${DB_BASENAME}_as_Subject_DB_$(date +%F_%H-%M)"
				blastn_result_csv_VAR="01_Identification/DMP_blastn_results/ev_${e_value}_ws_${word_size}_${DB_BASENAME}_as_Subject_DB_$(date +%F_%H-%M)/blastn_${DB_BASENAME}_VS_${DMP_fasta_BASENAME}_${e_value}_${word_size}.csv"
				blastn_result_txt_VAR="01_Identification/DMP_blastn_results/ev_${e_value}_ws_${word_size}_${DB_BASENAME}_as_Subject_DB_$(date +%F_%H-%M)/blastn_${DB_BASENAME}_VS_${DMP_fasta_BASENAME}_${e_value}_${word_size}.txt"

				echo "$DMP_fasta_BASENAME is being aligned with $DB_BASENAME."
				echo "Subject ID,Query ID,E-value,Percent Identity,Alignment Length,Query Length,Target Length,Query Coverage,Mismatches,Gaps,Bit Score,Max Score,Total Score" > "$blastn_result_csv_VAR"
				blastn -db "$DB" \
					-query "$DMP_fasta_FILE" \
					-evalue "$e_value" \
					-word_size "$word_size" \
					-outfmt "10 sseqid qseqid evalue pident length qlen slen qcovs mismatch gaps bitscore score score" >> "$blastn_result_csv_VAR"

				blastn -db "$DB" \
					-query "$DMP_fasta_FILE" \
					-evalue "$e_value" \
					-word_size "$word_size" > "$blastn_result_txt_VAR"

				echo -e "Alignment of $DMP_fasta_BASENAME with $DB_BASENAME. SUCCESS.\n"
			done
		done
	done
done

#BLASTN_FOR_LOOP


#: << 'MERGE_RESULT_INTO_SINGLE_CSV'

# Set the directory to search
directory="01_Identification/DMP_blastn_results"

# Output file
mkdir -p "01_Identification/DMP_curated_blastn_results"
output_file="01_Identification/DMP_curated_blastn_results/merged_output_run_$(date +%F_%H-%M).csv"

# Temporary file to store the header
header_file="01_Identification/header.tmp"

# Remove output file if it already exists
[ -f "$output_file" ] && rm "$output_file"

# Find all CSV files and process them
find "$directory" -type f -name "*.csv" | while read -r file; do
    # If the output file is empty, copy the header
    if [ ! -f "$header_file" ]; then
        head -n 1 "$file" > "$header_file"
        cat "$header_file" > "$output_file"
    fi

    # Append the file contents excluding the header
    tail -n +2 "$file" >> "$output_file"
done

# Remove the temporary header file
rm -f "$header_file"

echo "CSV files merged into $output_file."

#MERGE_RESULT_INTO_SINGLE_CSV
