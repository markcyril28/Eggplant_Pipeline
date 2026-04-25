#!/bin/bash
set -euo pipefail

: <<'DESCRIPTIONS'
Write the bash script for this. 

When writing the bash script for this, 
    Limit catching errors. 
    Keep the script as simple and concise as possible.
    Use comments to explain each step of the script. Ensure the script is easy to read and understand.
    Avoid using complex bash features or constructs.
    Work within the current directory.


# Instruction 1:
Copy the folder '000_Full_Datasets_GRF-GIF_curated' and all its contents 
    to a new folder named '000_Full_Datasets_GRF-GIF_curated_renamed'. 

    For each fasta file in the subfolders of '000_Full_Datasets_GRF-GIF_curated_renamed',
        copy the fasta file name to the end of every header in that fasta file. 

# Instruction 2:
For each the fasta in the subfolders of '000_Full_Datasets_GRF-GIF_curated_renamed'
    renamed the fasta file accordingly, 
        According to species (abbreviated), GRF or GIF, and type of sequence. 
            e.g. 
                Arabidopsis_thaliana/Araport11_cds_20240409_GIF_curated.fasta -> Arabidopsis_thaliana/AtGIF_cds.fasta
                Arabidopsis_thaliana/Araport11_cds_20240409_growth-regulating_factor_curated.fasta -> Arabidopsis_thaliana/AtGRF_cds.fasta
                Brachypodium_distachyon/GCF_000005505.3_Brachypodium_distachyon_v3.0_cds_from_genomic.fna_-interacting factor_curated.fasta -> Brachypodium_distachyon/BdGIF_cds.fasta
                Brachypodium_distachyon/GCF_000005505.3_Brachypodium_distachyon_v3.0_cds_from_genomic.fna_growth-regulating_factor_curated.fasta -> Brachypodium_distachyon/BdGRF_cds.fasta
            Reference (for species abbreviations): 
                Arabidopsis thaliana -> At
                Brachypodium distachyon -> Bd
                Oryza sativa -> Os
                Zea mays -> Zm
                Sorghum bicolor -> Sb
                Setaria italica -> Si
                Triticum aestivum -> Ta
                Hordeum vulgare -> Hv
            Reference (for GRF or GIF): 
                GRF -> GRF
                GIF -> GIF
                growth-regulating factor -> GRF
                interacting factor -> GIF
            Reference (for type of sequence):
                cds -> cds
                protein -> prot
                peptide -> pep
                mrna -> mrna
        Skip renaming if the fasta file is already named in the correct format.

DESCRIPTIONS

# ============================================
# Simple renamer for GRF/GIF FASTA datasets
# ============================================

# 1) Duplicate the dataset folder
cp -r 000_Full_Datasets_GRF-GIF_curated 000_Full_Datasets_GRF-GIF_curated_renamed

# 2) Append the file name to every FASTA header (>...) in each file
find 000_Full_Datasets_GRF-GIF_curated_renamed -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.faa" \) | while read -r file; do
  fname="$(basename "$file")"
  awk -v name="$fname" '/^>/{print $0," Extracted from:",name; next} {print}' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
done

# 3) Rename files to <SpeciesAbbr><GRF|GIF>_<seqtype>.fasta
#    Detect species, family, and sequence type.

# Species abbreviations
declare -A SP=(
  ["Arabidopsis_thaliana"]="At"
  ["Brachypodium_distachyon"]="Bd"
  ["Citrus_sinensis"]="Cs"
  ["Glycine_max"]="Gm"
  ["Gossypium_hirsutum"]="Gh"
  ["Lactuca_sativa"]="Ls"
  ["Oryza_sativa"]="Os"
  ["Physcomitrella_patens"]="Pp"
  ["Setaria_italica"]="Si"
  ["Solanum_lycopersicum"]="Sl"
  ["Triticum_aestivum"]="Ta"
  ["Vitis_vinifera"]="Vv"
  ["Zea_mays"]="Zm"
)

# Sequence type mapping
declare -A TMAP=( ["cds"]="cds" ["protein"]="prot" ["peptide"]="pep" ["mrna"]="mrna" )

find 000_Full_Datasets_GRF-GIF_curated_renamed -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.faa" \) | while read -r file; do
  dir="$(dirname "$file")"
  base="$(basename "$file")"
  species="$(basename "$dir")"

  # --- Detect species abbreviation ---
  abbr="${SP[$species]}"

  # Try stripping extra suffix (e.g., Oryza_sativa_osa → Oryza_sativa)
  if [[ -z "$abbr" ]]; then
    stripped="${species%_*}"  # remove last underscore section
    abbr="${SP[$stripped]}"
  fi

  # Fallback: detect from file name prefix
  if [[ -z "$abbr" ]]; then
    case "$base" in
      osa*|Osa*) abbr="Os" ;;   # Oryza sativa
      zm*|Zm*)   abbr="Zm" ;;   # Zea mays
      bd*|Bd*)   abbr="Bd" ;;   # Brachypodium distachyon
      sl*|Sl*)   abbr="Sl" ;;   # Solanum lycopersicum
      gm*|Gm*)   abbr="Gm" ;;   # Glycine max
      gh*|Gh*)   abbr="Gh" ;;   # Gossypium hirsutum
      *) abbr="" ;;
    esac
  fi

  # --- Detect GRF or GIF family ---
  family=""
  if [[ "$base" =~ (growth[-_ ]?regulating[-_ ]?factor|GRF) ]]; then
    family="GRF"
  elif [[ "$base" =~ (interacting[-_ ]?factor|GIF) ]]; then
    family="GIF"
  fi

  # --- Detect sequence type ---
  seq=""
  for k in cds protein peptide mrna; do
    [[ "$base" == *"$k"* ]] && seq="${TMAP[$k]}" && break
  done

  # --- Construct new filename ---
  if [[ -n "$abbr" && -n "$family" && -n "$seq" ]]; then
    # Always normalize extension to .fasta for consistency
    newpath="${dir}/${abbr}${family}_${seq}.fasta"
    [[ "$file" != "$newpath" ]] && mv "$file" "$newpath"
  else
    echo "Skipping (incomplete detection): $file"
  fi
done

echo "Done."
