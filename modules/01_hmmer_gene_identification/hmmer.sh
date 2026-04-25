#!/bin/bash
# Module: HMMER Profile-HMM Gene Identification (Dual Protein + Nucleotide)
# Usage: bash hmmer.sh --hmm <profile.hmm> --proteins <fasta> --transcripts <fasta> --evalue <val> --outdir <dir> [--threads N] [--alignment <seed_alignment>] [--skip-build] [--search-mode prot|nucl|both] [--config <config.toml>]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../logging/logging_utils.sh"

# Default values for standalone mode
CPU=12
E_VALUE="1e-10"
OUTPUT_DIR="."
HMM_PROFILE=""
PROTEINS_FASTA=""
TRANSCRIPTS_FASTA=""
ALIGNMENT_FILE=""
CDHIT_THRESHOLD=0.9
SKIP_BUILD=false
CONFIG_FILE=""
SEARCH_MODE="both"         # prot, nucl, or both
THREADS_FROM_CLI=false     # True when --threads is explicitly passed (overrides config)
USE_GA=true                # Use Pfam gathering thresholds
NHMMER_EVALUE="1e-5"       # Nucleotide search E-value
NHMMER_AUTO_BUILD=true     # Auto-build nucleotide HMM from protein hits
NHMMER_MIN_HITS=3          # Min transcript hits to build nucleotide HMM
CDHIT_AL=0.8               # Alignment coverage for longer sequence
CDHIT_AS=0.8               # Alignment coverage for shorter sequence
CDHIT_WORD_SIZE=5          # Word size for cd-hit protein (5 for >= 0.7 identity)
CDHIT_MEMORY=16000         # Memory limit in MB for cd-hit
CDHIT_CLUSTER_ALGO=1       # 0 = fast, 1 = accurate (slower)
CDHIT_EST_WORD_SIZE=8      # Word size for cd-hit-est
CDHIT_EST_MIN_LENGTH=30    # Minimum nucleotide sequence length
CDHIT_EST_MEMORY=16000     # Memory limit in MB for cd-hit-est
CDHIT_EST_CLUSTER_ALGO=1   # 0 = fast, 1 = accurate (slower)
CDHIT_EST_STRAND=0         # 0 = both strands, 1 = + strand only

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)              CONFIG_FILE="$2"; shift 2 ;;
        --hmm)                 HMM_PROFILE="$2"; shift 2 ;;
        --proteins)            PROTEINS_FASTA="$2"; shift 2 ;;
        --transcripts)         TRANSCRIPTS_FASTA="$2"; shift 2 ;;
        --evalue)              E_VALUE="$2"; shift 2 ;;
        --outdir)              OUTPUT_DIR="$2"; shift 2 ;;
        --threads)             CPU="$2"; THREADS_FROM_CLI=true; shift 2 ;;
        --alignment)           ALIGNMENT_FILE="$2"; shift 2 ;;
        --cdhit-threshold)     CDHIT_THRESHOLD="$2"; shift 2 ;;
        --skip-build)          SKIP_BUILD=true; shift ;;
        --search-mode)         SEARCH_MODE="$2"; shift 2 ;;
        --use-ga)              USE_GA="$2"; shift 2 ;;
        --nhmmer-evalue)       NHMMER_EVALUE="$2"; shift 2 ;;
        --nhmmer-auto-build)   NHMMER_AUTO_BUILD="$2"; shift 2 ;;
        --nhmmer-min-hits)     NHMMER_MIN_HITS="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# Load configuration from TOML if provided
if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    TOML_PARSER="$SCRIPT_DIR/../utils/parse_toml.py"
    # cfg <"subpath key"> <default>  — reads under [identification]; word-splits $1 for nested paths
    cfg() { # cfg "<subpath> key" default
        local val
        # shellcheck disable=SC2086
        val=$(python3 "$TOML_PARSER" "$CONFIG_FILE" identification $1 2>/dev/null) || true
        echo "${val:-$2}"
    }

    SEARCH_MODE=$(cfg "db_type"                          "$SEARCH_MODE")
    # Only read thread count from config if --threads was NOT passed by orchestrator.
    # The orchestrator calculates THREADS_PER_GENOME = total_cores / concurrent_genomes;
    # overriding it with the per-tool optimal would cause massive CPU over-subscription.
    if [[ "$THREADS_FROM_CLI" == "false" ]]; then
        CPU=$(cfg "hmmer_params threads" "$CPU")
    fi
    E_VALUE=$(cfg      "hmmer_params e_value"             "$E_VALUE")
    USE_GA=$(cfg       "hmmer_params use_ga"              "$USE_GA")

    CDHIT_THRESHOLD=$(cfg    "cdhit_params identity"      "$CDHIT_THRESHOLD")
    CDHIT_AL=$(cfg           "cdhit_params aL"            "$CDHIT_AL")
    CDHIT_AS=$(cfg           "cdhit_params aS"            "$CDHIT_AS")
    CDHIT_WORD_SIZE=$(cfg    "cdhit_params word_size"     "$CDHIT_WORD_SIZE")
    CDHIT_MEMORY=$(cfg       "cdhit_params memory"        "$CDHIT_MEMORY")
    CDHIT_CLUSTER_ALGO=$(cfg "cdhit_params cluster_algo"  "$CDHIT_CLUSTER_ALGO")

    CDHIT_EST_WORD_SIZE=$(cfg    "cdhit_est_params word_size"   "$CDHIT_EST_WORD_SIZE")
    CDHIT_EST_MIN_LENGTH=$(cfg   "cdhit_est_params min_length"  "$CDHIT_EST_MIN_LENGTH")
    CDHIT_EST_MEMORY=$(cfg       "cdhit_est_params memory"      "$CDHIT_EST_MEMORY")
    CDHIT_EST_CLUSTER_ALGO=$(cfg "cdhit_est_params cluster_algo" "$CDHIT_EST_CLUSTER_ALGO")
    CDHIT_EST_STRAND=$(cfg       "cdhit_est_params strand"      "$CDHIT_EST_STRAND")

    NHMMER_EVALUE=$(cfg       "nhmmer_params e_value"           "$NHMMER_EVALUE")
    NHMMER_AUTO_BUILD=$(cfg   "nhmmer_params auto_build"        "$NHMMER_AUTO_BUILD")
    NHMMER_MIN_HITS=$(cfg     "nhmmer_params min_hits_for_build" "$NHMMER_MIN_HITS")
fi

[[ -z "$HMM_PROFILE" ]]    && { log_error "Missing --hmm"; exit 1; }
[[ -z "$PROTEINS_FASTA" ]] && { log_error "Missing --proteins"; exit 1; }

PROFILE_BASE=$(basename "$HMM_PROFILE" .hmm)
HMMER_DIR="$OUTPUT_DIR/a_HMMER_RESULTS/$PROFILE_BASE"
RAW_DIR="$OUTPUT_DIR/b_RAW_EXTRACTED"
CDHIT_DIR="$OUTPUT_DIR/c_CD_HIT_Reduced/$PROFILE_BASE"
NHMMER_DIR="$OUTPUT_DIR/a_NHMMER_RESULTS/$PROFILE_BASE"
mkdir -p "$HMMER_DIR" "$RAW_DIR" "$CDHIT_DIR"

# Build/press HMM (skipped when orchestrator already did this)
if [[ "$SKIP_BUILD" == false ]]; then
    if [[ -n "$ALIGNMENT_FILE" && -f "$ALIGNMENT_FILE" ]]; then
        log_info "Building HMM profile for $PROFILE_BASE..."
        hmmbuild --cpu "$CPU" "$HMM_PROFILE" "$ALIGNMENT_FILE"
    fi
    log_info "Pressing HMM profile for $PROFILE_BASE..."
    hmmpress -f "$HMM_PROFILE"
fi

# ========================= PHASE 1: PROTEIN SEARCH =========================
PROT_HIT_IDS=""
if [[ "$SEARCH_MODE" == "prot" || "$SEARCH_MODE" == "both" ]]; then
    log_info "Running hmmsearch (protein) for $PROFILE_BASE (e-value <= $E_VALUE)..."

    # Build hmmsearch flags
    HMMSEARCH_FLAGS=(--cpu "$CPU")
    if [[ "$USE_GA" == "true" ]]; then
        HMMSEARCH_FLAGS+=(--cut_ga)
    fi
    HMMSEARCH_FLAGS+=(--seed 42 --noali --notextw)
    HMMSEARCH_FLAGS+=(--tblout "$HMMER_DIR/${PROFILE_BASE}_hits.tbl")
    HMMSEARCH_FLAGS+=(--domtblout "$HMMER_DIR/${PROFILE_BASE}_dom.tbl")

    hmmsearch "${HMMSEARCH_FLAGS[@]}" \
        "$HMM_PROFILE" \
        "$PROTEINS_FASTA" > "$HMMER_DIR/hmmsearch_${PROFILE_BASE}.log"

    # Filter by e-value (domain E-value, col 7 in domtblout)
    awk -v ev="$E_VALUE" '$1 !~ /^#/ && $7+0 <= ev+0 {print $0}' \
        "$HMMER_DIR/${PROFILE_BASE}_dom.tbl" > "$HMMER_DIR/${PROFILE_BASE}_hits_filtered.domtbl"

    HMM_RESULT="$HMMER_DIR/${PROFILE_BASE}_hits_filtered.domtbl"
    if [[ -s "$HMM_RESULT" ]]; then
        PROT_HIT_IDS="$HMMER_DIR/${PROFILE_BASE}_prot_hit_ids.txt"
        awk '{print $1}' "$HMM_RESULT" | sort -u > "$PROT_HIT_IDS"
        log_info "Protein search: $(wc -l < "$PROT_HIT_IDS") unique hits"
    else
        log_warn "Protein search: 0 hits for $PROFILE_BASE"
    fi
fi

# ======================== PHASE 2: NUCLEOTIDE SEARCH ========================
NUCL_HIT_IDS=""
if [[ ("$SEARCH_MODE" == "nucl" || "$SEARCH_MODE" == "both") && -n "$TRANSCRIPTS_FASTA" && -f "$TRANSCRIPTS_FASTA" ]]; then
    mkdir -p "$NHMMER_DIR"

    # Determine nucleotide HMM source
    NUCL_HMM=""

    # Option A: Auto-build nucleotide HMM from protein-search transcript hits
    if [[ "$NHMMER_AUTO_BUILD" == "true" && -n "$PROT_HIT_IDS" && -s "$PROT_HIT_IDS" ]]; then
        HIT_COUNT=$(wc -l < "$PROT_HIT_IDS")
        if (( HIT_COUNT >= NHMMER_MIN_HITS )); then
            log_info "Auto-building nucleotide HMM from $HIT_COUNT protein-hit transcripts..."

            # Extract transcripts for protein hits
            SEED_TRANSCRIPTS="$NHMMER_DIR/${PROFILE_BASE}_seed_transcripts.fa"
            awk '
                NR==FNR {ids[$1]=1; next}
                /^>/ {
                    if (p) print ""
                    p=0
                    split($0, a, " "); seqid=substr(a[1], 2)
                    if (seqid in ids) p=1
                    if (p) print
                    next
                }
                p {gsub(/\.$/, ""); print}
                END { if (p) print "" }
            ' "$PROT_HIT_IDS" "$TRANSCRIPTS_FASTA" > "$SEED_TRANSCRIPTS"

            if [[ -s "$SEED_TRANSCRIPTS" ]]; then
                # Align seed transcripts with MAFFT (fast, accurate for profile building)
                SEED_ALN="$NHMMER_DIR/${PROFILE_BASE}_seed_aln.fa"
                if command -v mafft &>/dev/null; then
                    mafft --auto --thread "$CPU" "$SEED_TRANSCRIPTS" > "$SEED_ALN" 2>/dev/null
                else
                    # Fallback: use unaligned sequences (hmmbuild handles this)
                    cp "$SEED_TRANSCRIPTS" "$SEED_ALN"
                fi

                # Build nucleotide HMM
                NUCL_HMM="$NHMMER_DIR/${PROFILE_BASE}_nucl.hmm"
                hmmbuild --cpu "$CPU" --dna "$NUCL_HMM" "$SEED_ALN"
                log_info "Nucleotide HMM built: $NUCL_HMM"
            else
                log_warn "No matching transcripts for nucleotide HMM build"
            fi
        else
            log_warn "Only $HIT_COUNT protein hits (need $NHMMER_MIN_HITS); skipping nucleotide HMM build"
        fi
    fi

    # Run nhmmer if we have a nucleotide HMM
    if [[ -n "$NUCL_HMM" && -f "$NUCL_HMM" ]]; then
        log_info "Running nhmmer (nucleotide) for $PROFILE_BASE (e-value <= $NHMMER_EVALUE)..."

        nhmmer --cpu "$CPU" \
            --seed 42 --noali --notextw \
            --tblout "$NHMMER_DIR/${PROFILE_BASE}_nhmmer_hits.tbl" \
            "$NUCL_HMM" \
            "$TRANSCRIPTS_FASTA" > "$NHMMER_DIR/nhmmer_${PROFILE_BASE}.log"

        # Filter nhmmer tblout: col 1 = target, col 13 = E-value
        awk -v ev="$NHMMER_EVALUE" '$1 !~ /^#/ && $13+0 <= ev+0 {print $1}' \
            "$NHMMER_DIR/${PROFILE_BASE}_nhmmer_hits.tbl" | sort -u \
            > "$NHMMER_DIR/${PROFILE_BASE}_nucl_hit_ids.txt"

        if [[ -s "$NHMMER_DIR/${PROFILE_BASE}_nucl_hit_ids.txt" ]]; then
            NUCL_HIT_IDS="$NHMMER_DIR/${PROFILE_BASE}_nucl_hit_ids.txt"
            log_info "Nucleotide search: $(wc -l < "$NUCL_HIT_IDS") unique hits"
        else
            log_warn "Nucleotide search: 0 hits for $PROFILE_BASE"
        fi
    fi
fi

# ========================= PHASE 3: MERGE HIT IDs ==========================
MERGED_IDS="$HMMER_DIR/${PROFILE_BASE}_hit_ids.txt"
if [[ -n "$PROT_HIT_IDS" && -s "$PROT_HIT_IDS" && -n "$NUCL_HIT_IDS" && -s "$NUCL_HIT_IDS" ]]; then
    # Union of protein and nucleotide hits
    sort -u "$PROT_HIT_IDS" "$NUCL_HIT_IDS" > "$MERGED_IDS"
    PROT_ONLY=$(comm -23 "$PROT_HIT_IDS" "$NUCL_HIT_IDS" | wc -l)
    NUCL_ONLY=$(comm -13 "$PROT_HIT_IDS" "$NUCL_HIT_IDS" | wc -l)
    BOTH_HITS=$(comm -12 "$PROT_HIT_IDS" "$NUCL_HIT_IDS" | wc -l)
    log_info "Merged hits: $(wc -l < "$MERGED_IDS") total (prot-only=$PROT_ONLY, nucl-only=$NUCL_ONLY, both=$BOTH_HITS)"
elif [[ -n "$PROT_HIT_IDS" && -s "$PROT_HIT_IDS" ]]; then
    cp "$PROT_HIT_IDS" "$MERGED_IDS"
    log_info "Using protein-only hits: $(wc -l < "$MERGED_IDS")"
elif [[ -n "$NUCL_HIT_IDS" && -s "$NUCL_HIT_IDS" ]]; then
    cp "$NUCL_HIT_IDS" "$MERGED_IDS"
    log_info "Using nucleotide-only hits: $(wc -l < "$MERGED_IDS")"
else
    log_warn "No hits found for $PROFILE_BASE in any search mode"
    log_info "HMMER workflow complete for $PROFILE_BASE (no hits)."
    exit 0
fi

HIT_COUNT=$(wc -l < "$MERGED_IDS")
log_info "Extracting $HIT_COUNT unique hits for $PROFILE_BASE (bulk)..."

# ======================= PHASE 4: EXTRACT SEQUENCES =========================
# Bulk extract proteins
awk '
    NR==FNR {ids[$1]=1; next}
    /^>/ {
        if (p) print ""
        p=0
        split($0, a, " "); seqid=substr(a[1], 2)
        if (seqid in ids) p=1
        if (p) print
        next
    }
    p {gsub(/\.$/, ""); print}
    END { if (p) print "" }
' "$MERGED_IDS" "$PROTEINS_FASTA" > "$RAW_DIR/${PROFILE_BASE}_matched_proteins.fa"

# Bulk extract transcripts
if [[ -n "$TRANSCRIPTS_FASTA" && -f "$TRANSCRIPTS_FASTA" ]]; then
    awk '
        NR==FNR {ids[$1]=1; next}
        /^>/ {
            if (p) print ""
            p=0
            split($0, a, " "); seqid=substr(a[1], 2)
            if (seqid in ids) p=1
            if (p) print
            next
        }
        p {gsub(/\.$/, ""); print}
        END { if (p) print "" }
    ' "$MERGED_IDS" "$TRANSCRIPTS_FASTA" > "$RAW_DIR/${PROFILE_BASE}_matched_transcripts.fa"
fi

# =================== PHASE 5: REDUNDANCY REDUCTION =========================
log_info "CD-HIT redundancy reduction (protein threshold=$CDHIT_THRESHOLD)..."

# Split threads between protein and nucleotide CD-HIT when both run in parallel
if [[ -f "$RAW_DIR/${PROFILE_BASE}_matched_transcripts.fa" && -s "$RAW_DIR/${PROFILE_BASE}_matched_transcripts.fa" ]]; then
    CDHIT_THREADS=$(( CPU / 2 ))
    (( CDHIT_THREADS < 1 )) && CDHIT_THREADS=1
else
    CDHIT_THREADS=$CPU
fi

# CD-HIT flags for protein
CDHIT_PROT_FLAGS=(-c "$CDHIT_THRESHOLD" -n "$CDHIT_WORD_SIZE" -T "$CDHIT_THREADS" -M "$CDHIT_MEMORY" -g "$CDHIT_CLUSTER_ALGO")
[[ "$CDHIT_AL" != "" ]] && CDHIT_PROT_FLAGS+=(-aL "$CDHIT_AL")
[[ "$CDHIT_AS" != "" ]] && CDHIT_PROT_FLAGS+=(-aS "$CDHIT_AS")

cd-hit -i "$RAW_DIR/${PROFILE_BASE}_matched_proteins.fa" \
    -o "$CDHIT_DIR/${PROFILE_BASE}_proteins_cdhit.fa" \
    "${CDHIT_PROT_FLAGS[@]}" &
CDHIT_PID_PROT=$!

# CD-HIT-EST for nucleotide sequences
if [[ -f "$RAW_DIR/${PROFILE_BASE}_matched_transcripts.fa" && -s "$RAW_DIR/${PROFILE_BASE}_matched_transcripts.fa" ]]; then
    log_info "CD-HIT-EST redundancy reduction (nucleotide threshold=$CDHIT_THRESHOLD)..."
    CDHIT_EST_FLAGS=(-c "$CDHIT_THRESHOLD" -n "$CDHIT_EST_WORD_SIZE" -T "$CDHIT_THREADS" -M "$CDHIT_EST_MEMORY" -g "$CDHIT_EST_CLUSTER_ALGO" -l "$CDHIT_EST_MIN_LENGTH" -r "$CDHIT_EST_STRAND")
    [[ "$CDHIT_AL" != "" ]] && CDHIT_EST_FLAGS+=(-aL "$CDHIT_AL")
    [[ "$CDHIT_AS" != "" ]] && CDHIT_EST_FLAGS+=(-aS "$CDHIT_AS")

    cd-hit-est -i "$RAW_DIR/${PROFILE_BASE}_matched_transcripts.fa" \
        -o "$CDHIT_DIR/${PROFILE_BASE}_transcripts_cdhit.fa" \
        "${CDHIT_EST_FLAGS[@]}" &
    CDHIT_PID_TRANS=$!
fi

wait "$CDHIT_PID_PROT" || { log_error "CD-HIT (proteins) failed"; exit 1; }
[[ -n "${CDHIT_PID_TRANS:-}" ]] && { wait "$CDHIT_PID_TRANS" || { log_error "CD-HIT-EST (transcripts) failed"; exit 1; }; }

log_info "HMMER workflow complete for $PROFILE_BASE."
