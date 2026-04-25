#!/bin/bash
set -euo pipefail
# NOTE: This is a legacy script. The active pipeline uses 2_BLAST_Alignment.sh + modules/01_identification/blastn.sh

# ===================== IMPORTANT VARIABLES =====================
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
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
	-query $AtPLAs_fasta_FILE \
	-evalue 1e-10 -word_size 11 > 01_Identification/blastn_V3_Chromosomes_against_AtPLAs_result.txt

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

IbPLA_fasta_FILE="01_Identification/query_fasta_file/Ipomea_batatas/IpPLAs_merged_fasta.fa"

# -------------------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------------------------
# Lists and Parameters
Eggplant_db_LIST=(
	#"$Eggplant_V3_Chromosomes_DB"
	#"$Eggplant_V3_CDS_putative_function_DB"
	#"$Eggplant_V4_DB"
	#"$Eggplant_V4_Chromosomes_DB"
	#"$Eggplant_V4_1_DB"
	"$Eggplant_V4_1_transcript_DB"
)

PLA_fasta_LIST=(
	"$IbPLA_fasta_FILE"
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

for PLA_fasta in ${PLA_fasta_LIST[@]}; do
	echo $PLA_fasta
done
DEBUGGER

#: << 'BLASTN_FOR_LOOP'
rm -rf 01_Identification/PLA_blastn_results/* 2>/dev/null || true

for DB in "${Eggplant_db_LIST[@]}"; do
	for PLA_fasta_FILE in "${PLA_fasta_LIST[@]}"; do
		for e_value in "${e_value_LIST[@]}"; do
			for word_size in "${word_size_LIST[@]}"; do
				DB_BASENAME=$(basename "$DB")
				PLA_fasta_BASENAME=$(basename "$PLA_fasta_FILE" .fasta)

				#mkdir -p "01_Identification/PLA_blastn_results/ev_${e_value}_ws_${word_size}"
				mkdir -p "01_Identification/PLA_blastn_results/ev_${e_value}_ws_${word_size}_${DB_BASENAME}_as_Subject_DB_$(date +%F_%H-%M)"
				blastn_result_csv_VAR="01_Identification/PLA_blastn_results/ev_${e_value}_ws_${word_size}_${DB_BASENAME}_as_Subject_DB_$(date +%F_%H-%M)/blastn_${DB_BASENAME}_VS_${PLA_fasta_BASENAME}_${e_value}_${word_size}.csv"
				blastn_result_txt_VAR="01_Identification/PLA_blastn_results/ev_${e_value}_ws_${word_size}_${DB_BASENAME}_as_Subject_DB_$(date +%F_%H-%M)/blastn_${DB_BASENAME}_VS_${PLA_fasta_BASENAME}_${e_value}_${word_size}.txt"

				echo "$PLA_fasta_BASENAME is being aligned with $DB_BASENAME."
				echo "Subject ID,Query ID,E-value,Percent Identity,Alignment Length,Query Length,Target Length,Query Coverage,Mismatches,Gaps,Bit Score,Max Score,Total Score" > "$blastn_result_csv_VAR"
				blastn -db "$DB" \
					-query "$PLA_fasta_FILE" \
					-evalue "$e_value" \
					-word_size "$word_size" \
					-outfmt "10 sseqid qseqid evalue pident length qlen slen qcovs mismatch gaps bitscore score score" >> "$blastn_result_csv_VAR"

				blastn -db "$DB" \
					-query "$PLA_fasta_FILE" \
					-evalue "$e_value" \
					-word_size "$word_size" > "$blastn_result_txt_VAR"

				echo -e "Alignment of $PLA_fasta_BASENAME with $DB_BASENAME. SUCCESS.\n"
			done
		done
	done
done

#BLASTN_FOR_LOOP


#: << 'MERGE_RESULT_INTO_SINGLE_CSV'

# Set the directory to search
directory="01_Identification/PLA_blastn_results"

# Output file
mkdir -p "01_Identification/PLA_curated_blastn_results"
output_file="01_Identification/PLA_curated_blastn_results/merged_output_run_$(date +%F_%H-%M).csv"

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
