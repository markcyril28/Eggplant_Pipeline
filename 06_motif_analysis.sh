#!/bin/bash
# ============================================================================
# Program 6: Motif Analysis
# ============================================================================
# Comment in/out the gene groups below, then run:
#   bash f_motif_analysis.sh
# ============================================================================

set -euo pipefail

# ===============================================================
# GENE_GROUPS, CPU, MAX_PARALLEL, OVERWRITE, and OPERATIONS are loaded from
# 06_motif_analysisCONFIG.toml  [pipeline]  section
# (gene-group overrides in config/{GROUP}/06_motif_analysis_{GROUP}.toml).
# Edit gene_groups in the TOML to control which groups run.

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

should_run() { [[ " ${OPERATIONS[@]} " =~ " $1 " ]]; }

load_pipeline_config() {
    MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
    CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || nproc 2>/dev/null || echo "12")
    MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || nproc 2>/dev/null || echo "4")
    OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "true")

    local ops_str
    ops_str=$(get_toml pipeline operations 2>/dev/null || true)
    if [[ -n "$ops_str" ]]; then
        # parse_toml.py outputs one item per line; mapfile -t splits on newlines
        mapfile -t OPERATIONS <<< "$ops_str"
    else
        OPERATIONS=("RUN_READY_DATABASE" "RUN_MEME_SUITE_Analysis")
    fi
}

resolve_plantcare_raw_dir() {
    local base_dir="$1"
    local configured_rel configured_abs fallback

    configured_rel=$(get_toml plantcare raw_results_dir 2>/dev/null || true)
    if [[ -n "$configured_rel" ]]; then
        configured_abs="$PIPELINE_DIR/$configured_rel"
        if [[ -d "$configured_abs" ]]; then
            echo "$configured_abs"
            return 0
        fi
    fi

    fallback=$(find "$base_dir/$GENE_SEQUENCE_OUTPUT_DIR" "$base_dir/05_Gene_Sequence_Analysis" -mindepth 2 -maxdepth 4 \
        -type d \( -name "$PLANTCARE_RAW_STAGE_DIR" -o -name "PlantCARE_Results" \) 2>/dev/null | head -n 1 || true)

    if [[ -n "$fallback" ]]; then
        echo "$fallback"
        return 0
    fi

    echo ""
}

TEMP_FILES=()
cleanup_temp() {
    for f in "${TEMP_FILES[@]+"${TEMP_FILES[@]}"}"; do rm -f "$f"; done
    safe_teardown_logging
}
trap cleanup_temp EXIT

GENE_SEQUENCE_OUTPUT_DIR="06_Motif_Analysis"
GTF_STAGE_DIR="01_GTF_Extraction"
FASTA_STAGE_DIR="02_FASTA_with_upstream_and_downstream"
PLANTCARE_STAGE_DIR="03_PlantCARE_Analysis"
PLANTCARE_RAW_STAGE_DIR="03_PlantCARE_Results"
MEME_STAGE_DIR="04_MEME_Analysis"

# ---------------------------------------------------------------------------
# Auto-discover phylo-ordered alignment file for the motif locations diagram.
# Search convention (first hit wins):
#   1. III_RESULT/<GROUP>/10_Secondary_Structure_Analysis/**/*<base>*phylo_ordered*.aln
#   2. III_RESULT/<GROUP>/05_Phylogenetics/**/*<base>*phylo_ordered*.aln
#   3. III_RESULT/<GROUP>/04_MSA/**/*<base>*AMINO_ACID*.aln          (alignment order fallback)
# <base> is the FASTA stem with alphabet tokens stripped, so AA + NT inputs
# resolve to the same base and share the AA-tree topology.
# Returns empty string when no candidate is found.
# ---------------------------------------------------------------------------
auto_find_phylo_order() {
    local group="$1" fasta_path="$2"
    local stem base
    stem=$(basename "$fasta_path")
    stem="${stem%.fasta}"
    stem="${stem%.fa}"
    base="$stem"
    base="${base/_AMINO_ACID_Sequence/}"
    base="${base/_NUCLEOTIDE_Sequence/}"
    base="${base%_amino_acid}"
    base="${base%_nucleotide_fixed}"
    base="${base%_nucleotide}"
    base="${base%_protein}"

    local -a search_dirs=(
        "$PIPELINE_DIR/III_RESULT/$group/10_Secondary_Structure_Analysis"
        "$PIPELINE_DIR/III_RESULT/$group/05_Phylogenetics"
    )
    local dir match
    for dir in "${search_dirs[@]}"; do
        [[ -d "$dir" ]] || continue
        match=$(find "$dir" -type f -name "*${base}*phylo_ordered*.aln" 2>/dev/null | sort | head -n 1)
        [[ -n "$match" ]] && { echo "$match"; return 0; }
    done

    # Fallback: any AMINO_ACID alignment from the MSA stage (not strictly phylo-ordered)
    local msa_dir="$PIPELINE_DIR/III_RESULT/$group/04_MSA"
    if [[ -d "$msa_dir" ]]; then
        match=$(find "$msa_dir" -type f -name "*${base}*AMINO_ACID*.aln" 2>/dev/null | sort | head -n 1)
        [[ -n "$match" ]] && { echo "$match"; return 0; }
    fi

    echo ""
}

resolve_top_folder() {
    local genome_path="$1"

    case "$genome_path" in
        *Solanum_melongena_v4.1*|*smel_v4_1*|*Eggplant*)
            echo "Solanum_melongena_v4.1"
            ;;
        *GPE001970*|*unito*)
            echo "GPE001970_SMEL5"
            ;;
        *)
            echo "shared"
            ;;
    esac
}

resolve_group_config() {
    local group="$1"
    local config_dir="$PIPELINE_DIR/config/${group}"
    local MERGE_TOML="$MODULES/utils/merge_toml.py"

    if [[ -d "$config_dir" ]]; then
        CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${group}_geneseq_cfg_XXXXXX.toml")
        TEMP_FILES+=("$CONFIG_FILE")

        # Required: 00_common_<GROUP>.toml provides [general] + [output_dirs] +
        # [reference]. Fail loudly if it's missing so the user sees which
        # file path is wrong, instead of a cryptic downstream "key not found".
        local common_file="$config_dir/00_common_${group}.toml"
        if [[ ! -f "$common_file" ]]; then
            echo "ERROR: required config not found: $common_file" >&2
            echo "       (expected naming: <stage>_<GROUP>.toml under config/<GROUP>/)" >&2
            ls -la "$config_dir" >&2
            exit 1
        fi

        # Deep-merge: shared defaults first, gene-group overrides win on conflict.
        # merge_toml.py logs skipped/missing files to stderr.
        python3 "$MERGE_TOML" \
            "$PIPELINE_DIR/06_motif_analysisCONFIG.toml" \
            "$common_file" \
            "$config_dir/06_motif_analysis_${group}.toml" \
            "$config_dir/07_plantcare_analysis_${group}.toml" \
            "$config_dir/01_hmmer_gene_identification_${group}.toml" \
            > "$CONFIG_FILE"
    else
        CONFIG_FILE="$PIPELINE_DIR/config/${group}.toml"
    fi
}

resolve_gene_sequence_dir() {
    local base_dir="$1"
    local configured_dir

    configured_dir=$(get_toml output_dirs gene_sequence 2>/dev/null || true)
    if [[ -n "$configured_dir" && "$configured_dir" != "$GENE_SEQUENCE_OUTPUT_DIR" ]]; then
        log_warn "Config output_dirs.gene_sequence is '$configured_dir'; using '$GENE_SEQUENCE_OUTPUT_DIR'"
    fi

    echo "$base_dir/$GENE_SEQUENCE_OUTPUT_DIR"
}

migrate_gene_sequence_dir_if_needed() {
    local base_dir="$1"
    local target_dir="$base_dir/$GENE_SEQUENCE_OUTPUT_DIR"
    local legacy_dir

    for legacy_dir in "05_Gene_Sequence" "05_Gene_Sequence_Anaysis" "05_Gene_Sequence_Analysis"; do
        local src="$base_dir/$legacy_dir"
        if [[ -d "$src" && "$src" != "$target_dir" ]]; then
            if [[ -d "$target_dir" ]]; then
                log_warn "Both '$src' and '$target_dir' exist; keeping existing target"
            else
                log_info "Renaming output folder: $legacy_dir -> $GENE_SEQUENCE_OUTPUT_DIR"
                mv "$src" "$target_dir"
            fi
        fi
    done
}

run_sequence_extraction() {
    local group="$1"
    local base_dir="$2"
    local upstream downstream

    should_run "Extraction_of_sequence_with_upstream_and_downstream" || return 0

    log_step "Gene Sequence Extraction: $group"

    upstream=$(get_toml gene_sequence upstream)
    downstream=$(get_toml gene_sequence downstream)
    local gene_seq_dir
    gene_seq_dir=$(resolve_gene_sequence_dir "$base_dir")
    mkdir -p "$gene_seq_dir"

    # Try legacy single gtf_file key (monolithic configs: GRF_GIF, PLA)
    local gtf_file_rel
    gtf_file_rel=$(get_toml gene_sequence gtf_file 2>/dev/null || true)

    if [[ -n "$gtf_file_rel" ]]; then
        local genome gtf_file top_folder seq_dir
        genome="$PIPELINE_DIR/$(get_toml reference eggplant_v4_1_genome)"
        gtf_file="$PIPELINE_DIR/$gtf_file_rel"
        top_folder=$(resolve_top_folder "$genome")
        seq_dir="$gene_seq_dir/$top_folder/$PLANTCARE_STAGE_DIR/$FASTA_STAGE_DIR"
        mkdir -p "$seq_dir"

        bash "$MODULES/06_motif_analysis/extract_sequences.sh" \
            --genome "$genome" \
            --gtf "$gtf_file" \
            --outdir "$seq_dir" \
            --upstream "$upstream" \
            --downstream "$downstream" \
            --threads "$CPU"
    else
        # Multi-genome mode: use GTFs generated by run_gtf_extraction
        local -a gtf_labels gtf_genomes
        mapfile -t gtf_labels < <(get_toml gene_sequence gtf_extraction labels)
        mapfile -t gtf_genomes < <(get_toml gene_sequence gtf_extraction genomes)

        for i in "${!gtf_labels[@]}"; do
            local genome_label="${gtf_labels[$i]}"
            local genome="$PIPELINE_DIR/${gtf_genomes[$i]}"
            local gtf_file="$gene_seq_dir/$genome_label/$PLANTCARE_STAGE_DIR/$GTF_STAGE_DIR/${genome_label}.gtf"
            local seq_dir="$gene_seq_dir/$genome_label/$PLANTCARE_STAGE_DIR/$FASTA_STAGE_DIR"

            if [[ ! -f "$gtf_file" ]]; then
                log_warn "No GTF found for $genome_label at $gtf_file"
                log_warn "Skipping extraction"
                continue
            fi

            log_info "Extracting sequences for $genome_label"
            bash "$MODULES/06_motif_analysis/extract_sequences.sh" \
                --genome "$genome" \
                --gtf "$gtf_file" \
                --outdir "$seq_dir" \
                --upstream "$upstream" \
                --downstream "$downstream" \
                --threads "$CPU"
        done
    fi

    log_step "Gene sequence extraction complete: $group"
}

run_plantcare_analysis() {
    local group="$1"
    local base_dir="$2"
    local raw_dir mod_dir dpi pc_out_rel pc_dir genes funcs
    local -a cmd

    # Check if any PlantCARE sub-step is enabled
    local any_plantcare=false
    for _op in "PlantCARE_post_processing" "PlantCARE_matrix_generation" \
               "PlantCARE_heatmap_generation" "PlantCARE_combined_heatmap" \
               "PlantCARE_groups_heatmap"; do
        should_run "$_op" && { any_plantcare=true; break; }
    done
    $any_plantcare || return 0

    log_step "PlantCARE Promoter Analysis: $group"

    # Build --steps from enabled sub-operations
    local -a step_list=()
    should_run "PlantCARE_post_processing"    && step_list+=("post-process")
    should_run "PlantCARE_matrix_generation"  && step_list+=("matrix")
    should_run "PlantCARE_heatmap_generation" && step_list+=("heatmap")
    should_run "PlantCARE_combined_heatmap"   && step_list+=("combined-heatmap")
    should_run "PlantCARE_groups_heatmap"     && step_list+=("groups-heatmap")
    local steps_str
    steps_str=$(IFS=','; echo "${step_list[*]:-post-process,matrix,heatmap}")

    raw_dir=$(resolve_plantcare_raw_dir "$base_dir")
    if [[ -z "$raw_dir" ]]; then
        log_error "PlantCARE raw results directory not found for $group"
        return 1
    fi

    if ! find "$raw_dir" -type f \( -name "*.tab" -o -name "*.tar.gz" \) | grep -q .; then
        log_error "No PlantCARE .tab or .tar.gz files found in: $raw_dir"
        return 1
    fi

    log_info "PlantCARE raw directory: $raw_dir"
    mod_dir="$MODULES/07_plantcare_analysis"
    dpi=$(get_toml plantcare dpi 2>/dev/null || echo "900")

    pc_out_rel=$(get_toml plantcare output_dir 2>/dev/null || true)
    if [[ -n "$pc_out_rel" ]]; then
        pc_dir="$PIPELINE_DIR/$pc_out_rel"
    else
        if [[ "$(basename "$raw_dir")" == "$PLANTCARE_RAW_STAGE_DIR" ]]; then
            pc_dir="$(dirname "$raw_dir")"
        else
            pc_dir="$(dirname "$raw_dir")/$PLANTCARE_STAGE_DIR"
        fi
    fi
    mkdir -p "$pc_dir"

    genes=$(get_toml plantcare selected_genes 2>/dev/null || true)
    funcs=$(get_toml plantcare selected_functions 2>/dev/null || true)

    cmd=(
        bash "$MODULES/07_plantcare_analysis/run_pipeline.sh"
        --raw-dir "$raw_dir"
        --outdir "$pc_dir"
        --module-dir "$mod_dir"
        --steps "$steps_str"
        --dpi "$dpi"
        --threads "$CPU"
    )

    [[ "$OVERWRITE" == true ]] && cmd+=(--overwrite)
    [[ -n "$genes" ]] && cmd+=(--genes "$(echo "$genes" | tr '\n' ',' | sed 's/,$//')")
    [[ -n "$funcs" ]] && cmd+=(--functions "$(echo "$funcs" | tr '\n' ',' | sed 's/,$//')")

    "${cmd[@]}"

    log_step "PlantCARE analysis complete: $group"
}

run_ready_database() {
    local group="$1"

    should_run "RUN_READY_DATABASE" || return 0

    log_step "MEME Database Setup: $group"

    local databases_dir
    databases_dir="$PIPELINE_DIR/$(get_toml meme databases_dir \
        2>/dev/null || echo "II_INPUTS/meme_motif_databases")"

    local overwrite_arg=""
    [[ "$OVERWRITE" == true ]] && overwrite_arg="--overwrite"

    bash "$MODULES/06_motif_analysis/meme_suite/setup_meme_databases.sh" \
        --outdir "$databases_dir" \
        $overwrite_arg

    log_step "MEME database setup complete: $group"
}

run_meme_analysis() {
    local group="$1"
    local base_dir="$2"

    should_run "RUN_MEME_SUITE_Analysis" || return 0

    log_step "MEME Motif Analysis: $group"

    local meme_enabled
    meme_enabled=$(get_toml meme enabled 2>/dev/null || echo "true")
    [[ "$meme_enabled" == "false" ]] && { log_info "MEME disabled in config — skipping"; return 0; }

    local gene_seq_dir
    gene_seq_dir=$(resolve_gene_sequence_dir "$base_dir")

    # Read all MEME parameters from TOML config
    local databases_dir steps nmotifs minw maxw mod
    local optimal_threads time_limit markov_order tomtom_dbs fimo_dbs
    databases_dir="$PIPELINE_DIR/$(get_toml meme databases_dir \
        2>/dev/null || echo "II_INPUTS/meme_motif_databases")"
    steps=$(get_toml meme steps 2>/dev/null || echo "meme,tomtom,fimo")
    nmotifs=$(get_toml meme nmotifs 2>/dev/null || echo "10")
    minw=$(get_toml meme minw 2>/dev/null || echo "6")
    maxw=$(get_toml meme maxw 2>/dev/null || echo "50")
    mod=$(get_toml meme mod 2>/dev/null || echo "anr")
    optimal_threads=$(get_toml meme optimal_threads 2>/dev/null || echo "$CPU")
    time_limit=$(get_toml meme time_limit 2>/dev/null || echo "300")
    markov_order=$(get_toml meme markov_order 2>/dev/null || echo "0")
    tomtom_dbs=$(get_toml meme tomtom_databases 2>/dev/null || \
        echo "JASPAR/JASPAR2024_CORE_plants_non-redundant_v2.meme,ARABD/ArabidopsisDAPv1.meme")
    fimo_dbs=$(get_toml meme fimo_databases 2>/dev/null || \
        echo "JASPAR/JASPAR2024_CORE_plants_non-redundant_v2.meme,ARABD/ArabidopsisDAPv1.meme")

    # JPEG / color scheme parameters
    local color_scheme jpeg_dpi jpeg_columns jpeg_quality jpeg_logo_padding motif_location_palette motif_location_bg motif_location_font_scale
    color_scheme=$(get_toml meme color_scheme 2>/dev/null || echo "default")
    jpeg_dpi=$(get_toml meme jpeg_dpi 2>/dev/null || echo "300")
    jpeg_columns=$(get_toml meme jpeg_columns 2>/dev/null || echo "5")
    jpeg_quality=$(get_toml meme jpeg_quality 2>/dev/null || echo "92")
    jpeg_logo_padding=$(get_toml meme jpeg_logo_padding 2>/dev/null || echo "10")
    motif_location_palette=$(get_toml meme motif_location_palette 2>/dev/null || echo "wong")
    motif_location_bg=$(get_toml meme motif_location_bg 2>/dev/null || echo "")
    motif_location_font_scale=$(get_toml meme motif_location_font_scale 2>/dev/null || echo "1.4")

    # MEME is multi-threaded: run floor(CPU / optimal_threads) genomes in parallel.
    # TOMTOM/FIMO are single-threaded: pass CPU as max-parallel inside the module.
    local parallel_meme_jobs=$(( CPU / optimal_threads ))
    (( parallel_meme_jobs < 1 )) && parallel_meme_jobs=1

    wait_for_meme_slot() { local limit="$1"; while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done; }

    local -a gtf_labels
    mapfile -t gtf_labels < <(get_toml gene_sequence gtf_extraction labels 2>/dev/null || true)

    if [[ ${#gtf_labels[@]} -eq 0 ]]; then
        log_warn "No gtf_extraction.labels found in config — skipping MEME"
        return 0
    fi

    local overwrite_arg=""
    [[ "$OVERWRITE" == true ]] && overwrite_arg="--overwrite"

    log_info "MEME threads per genome: $optimal_threads  |  parallel genome jobs: $parallel_meme_jobs"
    log_info "TOMTOM/FIMO max-parallel: $MAX_PARALLEL (single-threaded tools)"

    # ---------------------------------------------------------------------------
    # Priority 1: explicit [meme.inputs] — loop through pre-merged FASTA files
    # ---------------------------------------------------------------------------
    local -a meme_input_labels meme_input_fastas meme_input_phylo_orders
    mapfile -t meme_input_labels       < <(get_toml meme inputs labels 2>/dev/null || true)
    mapfile -t meme_input_fastas       < <(get_toml meme inputs fasta_files 2>/dev/null || true)
    mapfile -t meme_input_phylo_orders < <(get_toml meme inputs phylo_order_files 2>/dev/null || true)

    # Optional alphabet filter — when [meme].alphabets is set, only labels whose
    # detected alphabet appears in the list are processed. Default: both AA + NT.
    # Accepted tokens (case-insensitive): amino_acid, protein, nucleotide, dna.
    local -a _meme_alphabets _alph_filter=()
    mapfile -t _meme_alphabets < <(get_toml meme alphabets 2>/dev/null || true)
    if [[ ${#_meme_alphabets[@]} -eq 0 ]]; then
        _alph_filter=("amino_acid" "nucleotide")
    else
        local _a
        for _a in "${_meme_alphabets[@]}"; do
            case "${_a,,}" in
                amino_acid|protein|aa) _alph_filter+=("amino_acid") ;;
                nucleotide|dna|nt)     _alph_filter+=("nucleotide") ;;
                *) log_warn "Unknown [meme].alphabets token: '$_a' (skipping)" ;;
            esac
        done
    fi
    log_info "Alphabet filter: ${_alph_filter[*]}"

    if [[ ${#meme_input_labels[@]} -gt 0 ]]; then
        # Validate parallel arrays are same length
        if [[ ${#meme_input_labels[@]} -ne ${#meme_input_fastas[@]} ]]; then
            log_error "[meme.inputs] labels (${#meme_input_labels[@]}) and fasta_files (${#meme_input_fastas[@]}) count mismatch"
            log_error "Each label must have a corresponding fasta_files entry (comment/uncomment in pairs)"
            return 1
        fi
        log_info "Using explicit MEME inputs from [meme.inputs] (${#meme_input_labels[@]} dataset(s))"
        local -a meme_pids=()

        for i in "${!meme_input_labels[@]}"; do
            local dataset_label="${meme_input_labels[$i]}"
            local fasta_file="$PIPELINE_DIR/${meme_input_fastas[$i]}"

            # Apply alphabet filter — same regex as run_meme_pipeline.sh's
            # auto-detection so the orchestrator and module agree.
            local _label_alph="nucleotide"
            [[ "$dataset_label" =~ (amino_acid|protein) ]] && _label_alph="amino_acid"
            local _ok=false _fa
            for _fa in "${_alph_filter[@]}"; do
                [[ "$_fa" == "$_label_alph" ]] && { _ok=true; break; }
            done
            if ! $_ok; then
                log_info "Skip $dataset_label ($_label_alph not in alphabet filter)"
                continue
            fi

            # Output folder = label with alphabet suffix stripped, so AA + NT
            # variants share one folder (e.g. selected_GPE001970_{amino_acid,nucleotide}
            # both write into selected_GPE001970/04_MEME_Analysis/). Per-label file
            # prefixes inside the module keep AA/NT outputs from colliding.
            local output_label="$dataset_label"
            output_label="${output_label%_amino_acid}"
            output_label="${output_label%_nucleotide_fixed}"
            output_label="${output_label%_nucleotide}"
            output_label="${output_label%_protein}"
            output_label="${output_label%_aa}"
            output_label="${output_label%_nt}"
            output_label="${output_label%_dna}"
            local out_dir="$gene_seq_dir/$output_label/$MEME_STAGE_DIR"

            # Phylo-order resolution: auto-discover by convention, then let an
            # explicit non-empty [meme.inputs].phylo_order_files entry override.
            local phylo_order_file
            phylo_order_file=$(auto_find_phylo_order "$group" "$fasta_file")
            if [[ $i -lt ${#meme_input_phylo_orders[@]} ]]; then
                local _raw_phylo="${meme_input_phylo_orders[$i]}"
                [[ -n "$_raw_phylo" ]] && phylo_order_file="$PIPELINE_DIR/$_raw_phylo"
            fi
            [[ -n "$phylo_order_file" && ! -f "$phylo_order_file" ]] && {
                log_warn "Phylo-order file resolved but not on disk: $phylo_order_file"
                phylo_order_file=""
            }

            if [[ ! -f "$fasta_file" ]]; then
                log_warn "FASTA not found for $dataset_label at $fasta_file — skipping"
                continue
            fi

            log_info "Queuing MEME for $dataset_label"
            log_info "  FASTA  : $fasta_file"
            log_info "  Output : $out_dir"
            [[ -n "$phylo_order_file" ]] && log_info "  Phylo order: $phylo_order_file"

            local cmd=(
                bash "$MODULES/06_motif_analysis/meme_suite/run_meme_pipeline.sh"
                --fasta-file    "$fasta_file"
                --outdir        "$out_dir"
                --databases-dir "$databases_dir"
                --label         "$dataset_label"
                --steps         "$steps"
                --threads       "$optimal_threads"
                --max-parallel  "$MAX_PARALLEL"
                --nmotifs       "$nmotifs"
                --minw          "$minw"
                --maxw          "$maxw"
                --mod           "$mod"
                --time-limit    "$time_limit"
                --markov-order  "$markov_order"
                --tomtom-dbs    "$tomtom_dbs"
                --fimo-dbs      "$fimo_dbs"
                --color-scheme      "$color_scheme"
                --jpeg-dpi          "$jpeg_dpi"
                --jpeg-columns      "$jpeg_columns"
                --jpeg-quality      "$jpeg_quality"
                --jpeg-logo-padding "$jpeg_logo_padding"
                --motif-palette     "$motif_location_palette"
                --motif-bg          "$motif_location_bg"
                --motif-font-scale  "$motif_location_font_scale"
            )
            [[ -n "$phylo_order_file" ]] && cmd+=(--phylo-order-file "$phylo_order_file")
            [[ -n "$overwrite_arg" ]] && cmd+=("$overwrite_arg")

            wait_for_meme_slot "$parallel_meme_jobs"
            "${cmd[@]}" &
            meme_pids+=("$!:$dataset_label")
        done

        local meme_failed=0
        for entry in "${meme_pids[@]}"; do
            local pid="${entry%%:*}" label="${entry#*:}"
            if ! wait "$pid"; then
                log_error "MEME job failed for dataset: $label (PID $pid)"
                ((meme_failed++)) || true
            fi
        done
        if (( meme_failed > 0 )); then
            log_warn "$meme_failed of ${#meme_pids[@]} MEME dataset(s) failed — see errors above"
        fi
        log_step "MEME analysis complete: $group"
        return 0
    fi

    # ---------------------------------------------------------------------------
    # Fallback: derive FASTA dirs from upstream extraction stage (gtf_labels)
    # ---------------------------------------------------------------------------
    local -a meme_pids=()
    for genome_label in "${gtf_labels[@]}"; do
        # Input: FASTA sequences from the extraction stage
        local fasta_dir="$gene_seq_dir/$genome_label/$PLANTCARE_STAGE_DIR/$FASTA_STAGE_DIR"
        # Output: genome at top level, MEME as sibling of PlantCARE stage
        # Pattern: [output_stage]/[genome_name]/[specific_data_type]/[results]
        local out_dir="$gene_seq_dir/$genome_label/$MEME_STAGE_DIR"

        if [[ ! -d "$fasta_dir" ]]; then
            log_warn "No FASTA dir for $genome_label at $fasta_dir — skipping"
            continue
        fi

        log_info "Queuing MEME for $genome_label"
        log_info "  FASTA dir : $fasta_dir"
        log_info "  Output    : $out_dir"

        local cmd=(
            bash "$MODULES/06_motif_analysis/meme_suite/run_meme_pipeline.sh"
            --fasta-dir    "$fasta_dir"
            --outdir       "$out_dir"
            --databases-dir "$databases_dir"
            --label        "$genome_label"
            --steps        "$steps"
            --threads      "$optimal_threads"
            --max-parallel "$MAX_PARALLEL"
            --nmotifs      "$nmotifs"
            --minw         "$minw"
            --maxw         "$maxw"
            --mod          "$mod"
            --time-limit   "$time_limit"
            --markov-order "$markov_order"
            --tomtom-dbs   "$tomtom_dbs"
            --fimo-dbs     "$fimo_dbs"
            --color-scheme      "$color_scheme"
            --jpeg-dpi          "$jpeg_dpi"
            --jpeg-columns      "$jpeg_columns"
            --jpeg-quality      "$jpeg_quality"
            --jpeg-logo-padding "$jpeg_logo_padding"
            --motif-palette     "$motif_location_palette"
            --motif-bg          "$motif_location_bg"
            --motif-font-scale  "$motif_location_font_scale"
        )
        [[ -n "$overwrite_arg" ]] && cmd+=("$overwrite_arg")

        wait_for_meme_slot "$parallel_meme_jobs"
        "${cmd[@]}" &
        meme_pids+=("$!:$genome_label")
    done

    local meme_failed=0
    for entry in "${meme_pids[@]}"; do
        local pid="${entry%%:*}" label="${entry#*:}"
        if ! wait "$pid"; then
            log_error "MEME job failed for genome: $label (PID $pid)"
            ((meme_failed++)) || true
        fi
    done
    if (( meme_failed > 0 )); then
        log_warn "$meme_failed of ${#meme_pids[@]} MEME genome(s) failed — see errors above"
    fi

    log_step "MEME analysis complete: $group"
}

# ===========================================================================
# Combined JPEG: assemble per-dataset JPEG logos into one overview image
# showing motifs across ALL inputs in a single publication-ready figure.
# Gated by "RUN_JPEG_Export" operation.
# ===========================================================================
run_combined_jpeg() {
    local group="$1"
    local base_dir="$2"

    should_run "RUN_JPEG_Export" || return 0

    local gene_seq_dir
    gene_seq_dir=$(resolve_gene_sequence_dir "$base_dir")

    # Honor the same [meme].alphabets filter used in run_meme_analysis so that
    # commented-out alphabets don't get bundled in via stale on-disk JPEGs.
    local -a _meme_alphabets _alph_filter=()
    mapfile -t _meme_alphabets < <(get_toml meme alphabets 2>/dev/null || true)
    if [[ ${#_meme_alphabets[@]} -eq 0 ]]; then
        _alph_filter=("amino_acid" "nucleotide")
    else
        local _a
        for _a in "${_meme_alphabets[@]}"; do
            case "${_a,,}" in
                amino_acid|protein|aa) _alph_filter+=("amino_acid") ;;
                nucleotide|dna|nt)     _alph_filter+=("nucleotide") ;;
            esac
        done
    fi

    # Collect all per-dataset JPEG logo grids
    local -a jpeg_files=()
    local -a jpeg_labels=()
    while IFS= read -r jpg; do
        # Detect alphabet from the parent dir: .../05_JPEG/{amino_acid,nucleotide}/file.jpg
        local _parent _jpg_alph=""
        _parent=$(basename "$(dirname "$jpg")")
        case "$_parent" in
            amino_acid|nucleotide) _jpg_alph="$_parent" ;;
        esac

        # Filename fallback for the legacy flat layout (.../05_JPEG/<label>_motifs.jpg)
        if [[ -z "$_jpg_alph" ]]; then
            if [[ "$jpg" =~ (amino_acid|protein|_aa_) ]]; then
                _jpg_alph="amino_acid"
            elif [[ "$jpg" =~ (nucleotide|_nt_|_dna_) ]]; then
                _jpg_alph="nucleotide"
            fi
        fi

        if [[ -n "$_jpg_alph" ]]; then
            local _ok=false _fa
            for _fa in "${_alph_filter[@]}"; do
                [[ "$_fa" == "$_jpg_alph" ]] && { _ok=true; break; }
            done
            if ! $_ok; then
                log_info "Combined JPEG: skip $(basename "$jpg") ($_jpg_alph not in alphabet filter)"
                continue
            fi
        fi

        jpeg_files+=("$jpg")
        # Extract label from path: .../<label>/04_MEME_Analysis/05_JPEG/<label>_motifs.jpg
        local fname
        fname=$(basename "$jpg" _motifs.jpg)
        jpeg_labels+=("$fname")
    # Match both layouts:
    #   .../05_JPEG/<label>_motifs.jpg                         (legacy flat)
    #   .../05_JPEG/{amino_acid,nucleotide}/<label>_motifs.jpg (new alphabet-split)
    done < <(find "$gene_seq_dir" -path '*/05_JPEG/*' -name '*_motifs.jpg' 2>/dev/null | sort)

    if [[ ${#jpeg_files[@]} -lt 2 ]]; then
        log_info "Combined JPEG: fewer than 2 per-dataset JPEGs found — skipping combined view"
        return 0
    fi

    local combined_dir="$gene_seq_dir/05_Combined_JPEG"
    local combined_out="$combined_dir/${group}_all_motifs.jpg"

    if [[ "$OVERWRITE" != true && -f "$combined_out" ]]; then
        log_info "Combined JPEG exists, skipping (use OVERWRITE=true to redo): $combined_out"
        return 0
    fi

    # Prefer IMv7 'magick' over deprecated 'convert'
    local _im_cmd=""
    if command -v magick &>/dev/null; then
        _im_cmd="magick"
    elif command -v convert &>/dev/null; then
        _im_cmd="convert"
    fi

    if [[ -z "$_im_cmd" ]] && ! command -v montage &>/dev/null; then
        log_warn "ImageMagick not found — skipping combined JPEG"
        return 0
    fi

    mkdir -p "$combined_dir"
    log_step "Combined JPEG assembly: ${#jpeg_files[@]} datasets → $combined_out"

    local jpeg_dpi
    jpeg_dpi=$(get_toml meme jpeg_dpi 2>/dev/null || echo "600")
    local jpeg_quality
    jpeg_quality=$(get_toml meme jpeg_quality 2>/dev/null || echo "92")

    # Find a TTF/OTF font file on disk to pass as `-font /abs/path`. The
    # `magick -list font` route can list font *names* whose backing file is
    # not actually installed (then FreeType fails with `unable to read font`).
    # Searching the filesystem for an actual font file avoids that.
    local _font_path=""
    local _candidate
    for _candidate in \
        /usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf \
        /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf \
        /usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf \
        /usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf \
        /Library/Fonts/Helvetica.ttc \
        /System/Library/Fonts/Helvetica.ttc \
        /c/Windows/Fonts/arial.ttf \
        /mnt/c/Windows/Fonts/arial.ttf
    do
        [[ -f "$_candidate" ]] && { _font_path="$_candidate"; break; }
    done
    if [[ -z "$_font_path" ]]; then
        # Last resort: any .ttf or .otf under /usr/share/fonts.
        _font_path=$(find /usr/share/fonts -type f \( -name '*.ttf' -o -name '*.otf' \) 2>/dev/null | head -n 1 || true)
    fi

    if [[ -z "$_font_path" ]]; then
        log_warn "No TTF/OTF font found on host — combined JPEG will be assembled without per-dataset labels."
        log_warn "Install one with: sudo apt install fonts-dejavu-core"
    else
        log_info "Annotation font: $_font_path"
    fi

    # Label each per-dataset image (when a font is available), then stack
    # vertically. Without a font, just stack the originals.
    local -a labelled_parts=()
    local _annotate_log="$combined_dir/.annotate.log"
    : > "$_annotate_log"
    for i in "${!jpeg_files[@]}"; do
        local tmp_labelled="$combined_dir/.tmp_${jpeg_labels[$i]}.jpg"
        if [[ ! -f "${jpeg_files[$i]}" ]]; then
            log_warn "Per-dataset JPEG missing for ${jpeg_labels[$i]} — skipping"
            continue
        fi
        if [[ -z "$_font_path" ]]; then
            tmp_labelled="${jpeg_files[$i]}"
            labelled_parts+=("$tmp_labelled")
            continue
        fi
        if $_im_cmd "${jpeg_files[$i]}" \
            -gravity North \
            -background white \
            -splice 0x80 \
            -font "$_font_path" -fill black \
            -pointsize 48 \
            -annotate +0+15 "${jpeg_labels[$i]}" \
            -density "$jpeg_dpi" \
            "$tmp_labelled" 2>>"$_annotate_log"; then
            :
        else
            log_warn "Failed to annotate ${jpeg_labels[$i]} — using original (see $_annotate_log)"
            tmp_labelled="${jpeg_files[$i]}"
        fi
        labelled_parts+=("$tmp_labelled")
    done

    # Stack all labelled parts vertically
    $_im_cmd "${labelled_parts[@]}" \
        -background white \
        -gravity Center \
        -append \
        -density "$jpeg_dpi" \
        -quality "$jpeg_quality" \
        "$combined_out" 2>&1 | tee "$combined_dir/combined_montage.log" || true

    # Clean up temp files
    for f in "${labelled_parts[@]}"; do
        [[ "$f" == "$combined_dir/.tmp_"* ]] && rm -f "$f"
    done

    if [[ -f "$combined_out" ]]; then
        local sz
        sz=$(du -k "$combined_out" 2>/dev/null | cut -f1)
        log_info "Combined JPEG: $combined_out (${sz} KB) — ${#jpeg_files[@]} datasets"
    else
        log_warn "Combined JPEG not created — check $combined_dir/combined_montage.log"
    fi
}

run_secondary_structure_notice() {
    local group="$1"
    should_run "Secondary_structure_analysis" || return 0
    log_step "Secondary Structure Analysis: $group (not yet implemented)"
}

run_gtf_extraction() {
    local group="$1"
    local base_dir="$2"

    should_run "GTF_extraction" || return 0

    log_step "GTF Extraction from identification results: $group"

    local ident_dir="$base_dir/$(get_toml output_dirs identification)"
    local gene_seq_dir
    gene_seq_dir=$(resolve_gene_sequence_dir "$base_dir")
    local -a gtf_labels gtf_annotations gtf_genomes

    mapfile -t gtf_labels < <(get_toml gene_sequence gtf_extraction labels)
    mapfile -t gtf_annotations < <(get_toml gene_sequence gtf_extraction annotations)
    mapfile -t gtf_genomes < <(get_toml gene_sequence gtf_extraction genomes)

    local overwrite_arg=""
    [[ "$OVERWRITE" == true ]] && overwrite_arg="--overwrite"

    for i in "${!gtf_labels[@]}"; do
        local genome_label="${gtf_labels[$i]}"
        local annotation="$PIPELINE_DIR/${gtf_annotations[$i]}"
        local genome_dir="$ident_dir/$genome_label"

        if [[ ! -d "$genome_dir" ]]; then
            log_warn "No identification results for $genome_label — skipping"
            continue
        fi

        # Collect all HMMER hit ID files for this genome
        local -a hit_files=()
        while IFS= read -r f; do
            hit_files+=("$f")
        done < <(find "$genome_dir" -name "*_hit_ids.txt" 2>/dev/null)

        if [[ ${#hit_files[@]} -eq 0 ]]; then
            log_warn "No HMMER hit IDs found for $genome_label — skipping"
            continue
        fi

        # Combine all hit IDs into a single temp file
        local tmp_ids
        tmp_ids=$(mktemp "${TMPDIR:-/tmp}/${group}_${genome_label}_ids_XXXXXX.txt")
        TEMP_FILES+=("$tmp_ids")
        cat "${hit_files[@]}" | sort -u > "$tmp_ids"

        local num_ids
        num_ids=$(wc -l < "$tmp_ids")
        log_info "$genome_label: $num_ids unique hit IDs from ${#hit_files[@]} profile(s)"

        # Output directory per genome
        local out_dir
        out_dir="$gene_seq_dir/$genome_label/$PLANTCARE_STAGE_DIR/$GTF_STAGE_DIR"

        bash "$MODULES/06_motif_analysis/extract_gtf.sh" \
            --annotation "$annotation" \
            --ids "$tmp_ids" \
            --outdir "$out_dir" \
            --label "$genome_label" \
            $overwrite_arg

        log_info "$genome_label: GTF files -> $out_dir"
    done

    log_step "GTF extraction complete: $group"
}

process_gene_group() {
    local group="$1"
    local base_dir

    resolve_group_config "$group"
    load_pipeline_config
    base_dir="$PIPELINE_DIR/$(get_toml general base_dir)"
    migrate_gene_sequence_dir_if_needed "$base_dir"

    setup_logging

    run_gtf_extraction "$group" "$base_dir"
    run_sequence_extraction "$group" "$base_dir"
    run_ready_database "$group"
    run_meme_analysis "$group" "$base_dir"
    run_combined_jpeg "$group" "$base_dir"
    run_plantcare_analysis "$group" "$base_dir"
    run_secondary_structure_notice "$group"
    teardown_logging
}

# Load gene_groups from the parent shared TOML before any per-group merge,
# since we need the list to know which per-group configs to merge in turn.
SHARED_CONFIG="$PIPELINE_DIR/06_motif_analysisCONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null || true)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups not set in $SHARED_CONFIG" >&2
    exit 1
fi

for GENE_GROUP in "${GENE_GROUPS[@]}"; do
    process_gene_group "$GENE_GROUP"
done
