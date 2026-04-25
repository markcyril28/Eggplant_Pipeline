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

mkdir -p 01_Identification/query_fasta_file

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

# Databases (Date Accessed: September 22, 2024)
Eggplant_V3_Chromosomes_DB="$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_Chromosomes_db/Eggplant_V3_Chromosomes"
Eggplant_V3_CDS_putative_function_DB="$Eggplant_V3_assembly_chr_FASTA_DIR/Eggplant_V3_CDS_putative_function_db/Eggplant_V3_CDS_putative_function"

Eggplant_V4_DB="$Eggplant_V4_assembly_chr_FASTA_DIR/Eggplant_V4_db/Eggplant_V4"
Eggplant_V4_Chromosomes_DB="$Eggplant_V4_assembly_chr_FASTA_DIR/Eggplant_V4.chromosomes_db/Eggplant_V4.chromosomes"
Eggplant_V4_chr0="$Eggplant_V4_assembly_chr_FASTA_DIR/Eggplant_V4_chr0_db/Eggplant_V4_chr0"

Eggplant_V4_1_DB="$Eggplant_V4_1_fasta_FASTA_DIR/Eggplant_V4.1_db/Eggplant_V4.1"
Eggplant_V4_1_transcript_DB="$Eggplant_V4_1_fasta_FASTA_DIR/Eggplant_V4.1_transcripts.function_db/Eggplant_V4.1_transcripts.function"


# DMPs from Arabidopsis thaliana, Oryza sativa, and Zea mays (Date Accessed: September 22-26, 2024)

# Keyword searched of AtDMPs from (https://www.arabidopsis.org/).
AtDMPs_v1_fasta_FILE="01_Identification/query_fasta_file/Arabidopsis_thaliana/v1_keyword_searched/AtDMPs_nuc_seq_merged_v1.fasta"
# List of AtDMPs based on SweetPotato_paper, and nucleotide fasta sequences lifted from blastn datasets of Aradipopsis thalian from (https://www.arabidopsis.org/download/list?dir=Sequences%2FTAIR10_blastsets)
AtDMPs_v2_fasta_FILE="01_Identification/query_fasta_file/Arabidopsis_thaliana/v2_based_from_paper/AtDMPs_nuc_seq_merged_v2.fasta"

# OsDMPs from the query DMP sequences use in the paper of <paper SWEETPOTATO>, and sequences lifted from annotated sequences obtained from (https://rice.uga.edu/).
OsDMPs_fasta_FILE="01_Identification/query_fasta_file/Oryza_sativa/OsDMPs_nuc_seq_merged.fasta"

# ZmDMPs from a simple keyword search of "DMP" in this website, (https://www.maizegdb.org/gene_center/gene).
ZmDMPs_v1_fasta_FILE="01_Identification/query_fasta_file/Zea_mays/v1/ZmDMPs_nuc_seq_merged_v1.fasta"
# ZmDMPs from SweetPotato_paper, and searched/blastp aligned in NCBI. 
ZmDMPs_v2_fasta_FILE="01_Identification/query_fasta_file/Zea_mays/v2_based_on_paper_blastp/ZmDMPs_nuc_seq_merged_v2.fa"
# ZmDMPs_v2 fastas aligned using cds assembly fasta
ZmDMPs_v3_fasta_FILE=""

# Other_DMPs from SweetPotato_paper
OtherDMPs_fasta_FILE="01_Identification/query_fasta_file/Other_DMPs/OtherDMPs_merged_fasta.fa"


# Group 1 
# aquiegia coreulea (2024 Paper)
AqDMP_v1_fasta_FILE=""
# Physcomitrella patens (from 2024 Paper)




# Lists and Parameters
Eggplant_db_LIST=(
	#"$Eggplant_V3_Chromosomes_DB"
	#"$Eggplant_V3_CDS_putative_function_DB"
	#"$Eggplant_V4_DB"
	#"$Eggplant_V4_Chromosomes_DB"
	"$Eggplant_V4_1_DB"
	"$Eggplant_V4_1_transcript_DB"
)

DMP_fasta_LIST=(
	"$AtDMPs_v1_fasta_FILE"
	"$AtDMPs_v2_fasta_FILE"
	"$OsDMPs_fasta_FILE"
	"$ZmDMPs_v1_fasta_FILE"
	"$ZmDMPs_v2_fasta_FILE"
	"$OtherDMPs_fasta_FILE"
)

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
#rm -r 01_Identification/blastn_results/*

for DB in "${Eggplant_db_LIST[@]}"; do
	for DMP_fasta_FILE in "${DMP_fasta_LIST[@]}"; do
		for e_value in "${e_value_LIST[@]}"; do
			for word_size in "${word_size_LIST[@]}"; do
				DB_BASENAME=$(basename "$DB")
				DMP_fasta_BASENAME=$(basename "$DMP_fasta_FILE" .fasta)

				mkdir -p "01_Identification/blastn_results/ev_${e_value}_ws_${word_size}"
				mkdir -p "01_Identification/blastn_results/ev_${e_value}_ws_${word_size}/${DB_BASENAME}_as_Subject_DB"
				blastn_result_csv_VAR="01_Identification/blastn_results/ev_${e_value}_ws_${word_size}/${DB_BASENAME}_as_Subject_DB/blastn_${DB_BASENAME}_VS_${DMP_fasta_BASENAME}_${e_value}_${word_size}.csv"
				blastn_result_txt_VAR="01_Identification/blastn_results/ev_${e_value}_ws_${word_size}/${DB_BASENAME}_as_Subject_DB/blastn_${DB_BASENAME}_VS_${DMP_fasta_BASENAME}_${e_value}_${word_size}.txt"

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

				echo -e "Alignment of $DMP_fasta_BASENAME with $DB_BASENAME is a success.\n"
			done
		done
	done
done

#BLASTN_FOR_LOOP
