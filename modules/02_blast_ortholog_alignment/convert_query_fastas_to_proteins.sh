#!/bin/bash
# Convert configured BLAST nucleotide query FASTAs into protein FASTAs.
# Reads pairwise paths from [ortholog_blast].query_fastas and query_protein_fastas.

set -euo pipefail

GENE_GROUP="DMP"
THREADS=1
OVERWRITE=false
CONFIG_FILE=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
MODULES="$PIPELINE_DIR/modules"
TOML_PARSER="$MODULES/utils/parse_toml.py"
TRANSLATE_MODULE="$MODULES/08_protein_structure_analysis/translate.sh"

usage() {
	cat <<'EOF'
Usage:
  bash modules/01_identification/convert_query_fastas_to_proteins.sh [options]

Options:
  --gene-group <name>   Gene group name for split config lookup (default: DMP)
  --config <path>       Explicit TOML file path (overrides --gene-group)
  --threads <int>       Threads passed to translator module (default: 1)
  --overwrite           Rebuild protein FASTAs even if output exists
  --help                Show this help text
EOF
}

while [[ $# -gt 0 ]]; do
	case "$1" in
		--gene-group)
			GENE_GROUP="$2"
			shift 2
			;;
		--config)
			CONFIG_FILE="$2"
			shift 2
			;;
		--threads)
			THREADS="$2"
			shift 2
			;;
		--overwrite)
			OVERWRITE=true
			shift
			;;
		--help|-h)
			usage
			exit 0
			;;
		*)
			echo "ERROR: Unknown argument: $1" >&2
			usage
			exit 1
			;;
	esac
done

if [[ -z "$CONFIG_FILE" ]]; then
	CANDIDATE="$PIPELINE_DIR/config/$GENE_GROUP/02_blast_ortholog_alignment_${GENE_GROUP}.toml"
	if [[ -f "$CANDIDATE" ]]; then
		CONFIG_FILE="$CANDIDATE"
	else
		CONFIG_FILE="$PIPELINE_DIR/config/$GENE_GROUP.toml"
	fi
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
	echo "ERROR: Config file not found: $CONFIG_FILE" >&2
	exit 1
fi

if [[ ! -f "$TRANSLATE_MODULE" ]]; then
	echo "ERROR: Translation module not found: $TRANSLATE_MODULE" >&2
	exit 1
fi

HAS_TRANSEQ=true
if ! command -v transeq >/dev/null 2>&1; then
	HAS_TRANSEQ=false
fi

get_toml() {
	python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"
}

read_toml_list() {
	local section="$1"
	local key="$2"
	local raw

	raw=$(get_toml "$section" "$key" 2>/dev/null || true)
	[[ -z "${raw// }" ]] && return 0

	# parse_toml.py emits one item per line; mapfile preserves paths with spaces
	local items=()
	mapfile -t items <<< "$raw"
	printf '%s\n' "${items[@]}"
}

resolve_path() {
	local p="$1"
	if [[ "$p" == /* || "$p" =~ ^[A-Za-z]:[\\/] ]]; then
		printf '%s\n' "$p"
	else
		printf '%s\n' "$PIPELINE_DIR/$p"
	fi
}

mapfile -t NUC_FASTAS < <(read_toml_list ortholog_blast query_fastas)
mapfile -t PROT_FASTAS < <(read_toml_list ortholog_blast query_protein_fastas)

if [[ ${#NUC_FASTAS[@]} -eq 0 ]]; then
	echo "ERROR: No entries found for ortholog_blast.query_fastas in $CONFIG_FILE" >&2
	exit 1
fi

if [[ ${#PROT_FASTAS[@]} -eq 0 ]]; then
	echo "ERROR: No entries found for ortholog_blast.query_protein_fastas in $CONFIG_FILE" >&2
	exit 1
fi

if [[ ${#NUC_FASTAS[@]} -ne ${#PROT_FASTAS[@]} ]]; then
	echo "ERROR: query_fastas count (${#NUC_FASTAS[@]}) != query_protein_fastas count (${#PROT_FASTAS[@]})" >&2
	exit 1
fi

echo "Config: $CONFIG_FILE"
echo "Pairs to process: ${#NUC_FASTAS[@]}"
echo "Threads per translation: $THREADS"
echo "Overwrite existing outputs: $OVERWRITE"
echo "Translator backend: $([[ "$HAS_TRANSEQ" == "true" ]] && echo "transeq" || echo "python_fallback")"

converted=0
skipped=0
missing=0

for i in "${!NUC_FASTAS[@]}"; do
	nuc_rel="${NUC_FASTAS[$i]}"
	prot_rel="${PROT_FASTAS[$i]}"

	nuc_abs="$(resolve_path "$nuc_rel")"
	prot_abs="$(resolve_path "$prot_rel")"

	if [[ ! -f "$nuc_abs" ]]; then
		echo "WARN: Missing nucleotide FASTA: $nuc_abs"
		((missing+=1))
		continue
	fi

	if [[ -f "$prot_abs" && "$OVERWRITE" != "true" ]]; then
		echo "SKIP: Exists: $prot_abs"
		((skipped+=1))
		continue
	fi

	echo "RUN : $(basename "$nuc_abs") -> $(basename "$prot_abs")"
	if [[ "$HAS_TRANSEQ" == "true" ]]; then
		bash "$TRANSLATE_MODULE" --input "$nuc_abs" --output "$prot_abs" --threads "$THREADS"
	else
		mkdir -p "$(dirname "$prot_abs")"
		python3 - "$nuc_abs" "$prot_abs" <<'PY'
import sys

in_path, out_path = sys.argv[1], sys.argv[2]

CODON_TABLE = {
	"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
	"CTT":"L","CTC":"L","CTA":"L","CTG":"L",
	"ATT":"I","ATC":"I","ATA":"I","ATG":"M",
	"GTT":"V","GTC":"V","GTA":"V","GTG":"V",
	"TCT":"S","TCC":"S","TCA":"S","TCG":"S",
	"CCT":"P","CCC":"P","CCA":"P","CCG":"P",
	"ACT":"T","ACC":"T","ACA":"T","ACG":"T",
	"GCT":"A","GCC":"A","GCA":"A","GCG":"A",
	"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
	"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
	"AAT":"N","AAC":"N","AAA":"K","AAG":"K",
	"GAT":"D","GAC":"D","GAA":"E","GAG":"E",
	"TGT":"C","TGC":"C","TGA":"*","TGG":"W",
	"CGT":"R","CGC":"R","CGA":"R","CGG":"R",
	"AGT":"S","AGC":"S","AGA":"R","AGG":"R",
	"GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate_frame0(seq: str) -> str:
	seq = seq.upper().replace("U", "T")
	aa = []
	end = len(seq) - (len(seq) % 3)
	for i in range(0, end, 3):
		codon = seq[i:i+3]
		if any(base not in "ACGT" for base in codon):
			aa.append("X")
		else:
			aa.append(CODON_TABLE.get(codon, "X"))
	return "".join(aa).rstrip("*")

def write_wrapped(handle, seq: str, width: int = 60) -> None:
	for i in range(0, len(seq), width):
		handle.write(seq[i:i+width] + "\n")

count = 0
with open(in_path, "r", encoding="utf-8") as fin, open(out_path, "w", encoding="utf-8") as fout:
	header = None
	seq_parts = []

	for raw in fin:
		line = raw.strip()
		if not line:
			continue
		if line.startswith(">"):
			if header is not None:
				protein = translate_frame0("".join(seq_parts))
				fout.write(header + "\n")
				write_wrapped(fout, protein)
				count += 1
			header = line
			seq_parts = []
		else:
			seq_parts.append(line)

	if header is not None:
		protein = translate_frame0("".join(seq_parts))
		fout.write(header + "\n")
		write_wrapped(fout, protein)
		count += 1

print(f"Translated {count} sequences: {in_path} -> {out_path}")
PY
	fi
	((converted+=1))
done

echo "Done. converted=$converted skipped=$skipped missing_inputs=$missing"

