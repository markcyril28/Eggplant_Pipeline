#!/bin/bash
# ============================================================================
# Program 7: PlantCARE Promoter Analysis
# ============================================================================
# Extracts GTF entries for identified genes, extracts upstream/downstream
# promoter sequences, then runs PlantCARE cis-element analysis and
# heatmap visualisation.
#
# Edit gene_groups in 07_plant_care_analysisCONFIG.toml, then run:
#   bash g_plant_care_analysis.sh
# ============================================================================

set -euo pipefail

# ===================== IMPORTANT VARIABLES =====================
# All variables (GENE_GROUPS, CPU, OVERWRITE, OPERATIONS) are loaded from
# 07_plant_care_analysisCONFIG.toml — edit them there.
# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
TOML_MERGER="$MODULES/utils/merge_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

# Load GENE_GROUPS from shared config (read before the per-group loop)
SHARED_CONFIG="$PIPELINE_DIR/07_plant_care_analysisCONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups is empty in 07_plant_care_analysisCONFIG.toml" >&2
    exit 1
fi

should_run() { [[ " ${OPERATIONS[@]} " =~ " $1 " ]]; }

load_operations_from_config() {
    local ops_str
    ops_str=$(get_toml plantcare operations 2>/dev/null || true)
    if [[ -n "$ops_str" ]]; then
        mapfile -t OPERATIONS <<< "$ops_str"
    fi

    # CPU threads from TOML (machine/compute profile)
    local machine_val cpu_val
    machine_val=$(get_toml plantcare machine 2>/dev/null || true)
    machine_val="${machine_val:-Local}"
    cpu_val=$(get_toml plantcare compute "$machine_val" threads 2>/dev/null || true)
    [[ -n "$cpu_val" ]] && CPU="$cpu_val"

    # Overwrite flag from TOML (overrides the top-level default)
    local ow
    ow=$(get_toml plantcare overwrite 2>/dev/null || true)
    [[ -n "$ow" ]] && OVERWRITE="$ow"
}

TEMP_FILES=()
cleanup_temp() {
    for f in "${TEMP_FILES[@]+"${TEMP_FILES[@]}"}"; do rm -f "$f"; done
    safe_teardown_logging
}
trap cleanup_temp EXIT

PLANTCARE_OUTPUT_DIR="07_PlantCARE_Analysis"
GTF_STAGE_DIR="01_GTF_Extraction"
FASTA_STAGE_DIR="02_FASTA_with_upstream_and_downstream"
PLANTCARE_RAW_STAGE_DIR="03_PlantCARE_Results"
GENE_SEQUENCE_OUTPUT_DIR="06_Motif_Analysis"

# ---------- helpers ----------------------------------------------------------

resolve_top_folder() {
    local genome_path="$1"
    case "$genome_path" in
        *Solanum_melongena_v4.1*|*smel_v4_1*|*Eggplant*)
            echo "Solanum_melongena_v4.1" ;;
        *GPE001970*|*unito*)
            echo "GPE001970_SMEL5" ;;
        *)
            echo "shared" ;;
    esac
}

resolve_group_config() {
    local group="$1"
    local config_dir="$PIPELINE_DIR/config/${group}"
    local shared_dir="$PIPELINE_DIR/config/shared"

    if [[ -d "$config_dir" ]]; then
        CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${group}_plantcare_cfg_XXXXXX.toml")
        TEMP_FILES+=("$CONFIG_FILE")

        # Build list of TOML files to deep-merge (shared defaults first,
        # then group-specific overrides).  merge_toml.py silently skips
        # missing files.
        local -a toml_sources=(
            "$shared_dir/00_common.toml"
            "$PIPELINE_DIR/07_plant_care_analysisCONFIG.toml"
            "$config_dir/00_common_${group}.toml"
            "$config_dir/07_plantcare_analysis_${group}.toml"
        )
        [[ -f "$PIPELINE_DIR/01_hmmer_identifyCONFIG.toml" ]] && \
            toml_sources+=("$PIPELINE_DIR/01_hmmer_identifyCONFIG.toml")
        [[ -f "$config_dir/01_hmmer_gene_identification_${group}.toml" ]] && \
            toml_sources+=("$config_dir/01_hmmer_gene_identification_${group}.toml")

        python3 "$TOML_MERGER" "${toml_sources[@]}" > "$CONFIG_FILE"
    else
        CONFIG_FILE="$PIPELINE_DIR/config/${group}.toml"
    fi
}

resolve_plantcare_dir() {
    local base_dir="$1"
    echo "$base_dir/$PLANTCARE_OUTPUT_DIR"
}

resolve_plantcare_raw_dir() {
    local base_dir="$1"
    local pc_dir="$2"
    local configured_rel configured_abs fallback

    configured_rel=$(get_toml plantcare raw_results_dir 2>/dev/null || true)
    if [[ -n "$configured_rel" ]]; then
        configured_abs="$PIPELINE_DIR/$configured_rel"
        if [[ -d "$configured_abs" ]]; then
            echo "$configured_abs"
            return 0
        fi
    fi

    # Search inside the PlantCARE output dir first, then legacy locations
    fallback=$(find "$pc_dir" "$base_dir/$GENE_SEQUENCE_OUTPUT_DIR" \
        "$base_dir/05_Gene_Sequence_Analysis" "$base_dir/07_PlantCARE_Analysis" \
        -mindepth 1 -maxdepth 5 \
        -type d \( -name "$PLANTCARE_RAW_STAGE_DIR" -o -name "PlantCARE_Results" \) 2>/dev/null | head -n 1)

    if [[ -n "$fallback" ]]; then
        echo "$fallback"
        return 0
    fi

    echo ""
}

# ---------- stage functions --------------------------------------------------

run_gtf_extraction() {
    local group="$1"
    local base_dir="$2"
    local pc_dir="$3"

    should_run "GTF_extraction" || return 0

    log_step "GTF Extraction from identification results: $group"

    local ident_dir="$base_dir/$(get_toml output_dirs identification)"
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
            genome_dir="$ident_dir/shared/$genome_label"
        fi
        if [[ ! -d "$genome_dir" ]]; then
            log_warn "No identification results for $genome_label — skipping"
            continue
        fi

        local out_dir="$pc_dir/$genome_label/$GTF_STAGE_DIR"

        # If d_GENES/ exists, produce one GTF per gene family; otherwise combine all.
        local d_genes_dir
        d_genes_dir=$(find "$genome_dir" -type d -name "d_GENES" 2>/dev/null | head -n 1)

        if [[ -n "$d_genes_dir" ]]; then
            local -a family_dirs=()
            while IFS= read -r d; do
                family_dirs+=("$d")
            done < <(find "$d_genes_dir" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)

            for family_dir in "${family_dirs[@]}"; do
                local family
                family=$(basename "$family_dir")
                local -a hit_files=()
                while IFS= read -r f; do
                    hit_files+=("$f")
                done < <(find "$family_dir" -name "*_hit_ids.txt" 2>/dev/null)
                if [[ ${#hit_files[@]} -eq 0 ]]; then
                    log_warn "$genome_label/$family: no hit IDs — skipping"
                    continue
                fi
                local tmp_ids
                tmp_ids=$(mktemp "${TMPDIR:-/tmp}/${group}_${genome_label}_${family}_ids_XXXXXX.txt")
                TEMP_FILES+=("$tmp_ids")
                cat "${hit_files[@]}" | sort -u > "$tmp_ids"
                log_info "$genome_label/$family: $(wc -l < "$tmp_ids") unique hit IDs"
                bash "$MODULES/06_motif_analysis/extract_gtf.sh" \
                    --annotation "$annotation" \
                    --ids "$tmp_ids" \
                    --outdir "$out_dir/$family" \
                    --label "$genome_label" \
                    $overwrite_arg
                log_info "$genome_label/$family: GTF -> $out_dir/$family"
            done
        else
            local -a hit_files=()
            while IFS= read -r f; do
                hit_files+=("$f")
            done < <(find "$genome_dir" -name "*_hit_ids.txt" 2>/dev/null)
            if [[ ${#hit_files[@]} -eq 0 ]]; then
                log_warn "No HMMER hit IDs found for $genome_label — skipping"
                continue
            fi
            local tmp_ids
            tmp_ids=$(mktemp "${TMPDIR:-/tmp}/${group}_${genome_label}_ids_XXXXXX.txt")
            TEMP_FILES+=("$tmp_ids")
            cat "${hit_files[@]}" | sort -u > "$tmp_ids"
            log_info "$genome_label: $(wc -l < "$tmp_ids") unique hit IDs from ${#hit_files[@]} profile(s)"
            bash "$MODULES/06_motif_analysis/extract_gtf.sh" \
                --annotation "$annotation" \
                --ids "$tmp_ids" \
                --outdir "$out_dir" \
                --label "$genome_label" \
                $overwrite_arg
            log_info "$genome_label: GTF -> $out_dir"
        fi
    done

    log_step "GTF extraction complete: $group"
}

run_sequence_extraction() {
    local group="$1"
    local base_dir="$2"
    local pc_dir="$3"
    local upstream downstream

    should_run "Extraction_of_sequence_with_upstream_and_downstream" || return 0

    log_step "Gene Sequence Extraction: $group"

    upstream=$(get_toml gene_sequence upstream)
    downstream=$(get_toml gene_sequence downstream)

    # Try legacy single gtf_file key (monolithic configs: GRF_GIF, PLA)
    local gtf_file_rel
    gtf_file_rel=$(get_toml gene_sequence gtf_file 2>/dev/null || true)

    if [[ -n "$gtf_file_rel" ]]; then
        local genome gtf_file top_folder seq_dir
        genome="$PIPELINE_DIR/$(get_toml reference eggplant_v4_1_genome)"
        gtf_file="$PIPELINE_DIR/$gtf_file_rel"
        top_folder=$(resolve_top_folder "$genome")
        seq_dir="$pc_dir/$top_folder/$FASTA_STAGE_DIR"
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
            local gtf_base_dir="$pc_dir/$genome_label/$GTF_STAGE_DIR"
            local seq_base_dir="$pc_dir/$genome_label/$FASTA_STAGE_DIR"

            # Check for per-family GTF subdirs produced by run_gtf_extraction
            local -a family_gtf_dirs=()
            while IFS= read -r d; do
                [[ -f "$d/${genome_label}.gtf" ]] && family_gtf_dirs+=("$d")
            done < <(find "$gtf_base_dir" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)

            if [[ ${#family_gtf_dirs[@]} -gt 0 ]]; then
                for family_dir in "${family_gtf_dirs[@]}"; do
                    local family
                    family=$(basename "$family_dir")
                    local gtf_file="$family_dir/${genome_label}.gtf"
                    local seq_dir="$seq_base_dir/$family"
                    mkdir -p "$seq_dir"
                    log_info "Extracting sequences for $genome_label/$family"
                    bash "$MODULES/06_motif_analysis/extract_sequences.sh" \
                        --genome "$genome" \
                        --gtf "$gtf_file" \
                        --outdir "$seq_dir" \
                        --upstream "$upstream" \
                        --downstream "$downstream" \
                        --threads "$CPU"
                done
            else
                local gtf_file="$gtf_base_dir/${genome_label}.gtf"
                if [[ ! -f "$gtf_file" ]]; then
                    log_warn "No GTF found for $genome_label at $gtf_file"
                    log_warn "Skipping extraction"
                    continue
                fi
                mkdir -p "$seq_base_dir"
                log_info "Extracting sequences for $genome_label"
                bash "$MODULES/06_motif_analysis/extract_sequences.sh" \
                    --genome "$genome" \
                    --gtf "$gtf_file" \
                    --outdir "$seq_base_dir" \
                    --upstream "$upstream" \
                    --downstream "$downstream" \
                    --threads "$CPU"
            fi
        done
    fi

    log_step "Gene sequence extraction complete: $group"
}

run_plantcare_analysis() {
    local group="$1"
    local base_dir="$2"
    local pc_dir="$3"

    # Check if any PlantCARE analysis sub-step is enabled
    local any_plantcare=false
    for _op in "PlantCARE_post_processing" "PlantCARE_matrix_generation" \
               "PlantCARE_heatmap_generation" "PlantCARE_combined_heatmap" \
               "PlantCARE_groups_heatmap"; do
        should_run "$_op" && { any_plantcare=true; break; }
    done
    $any_plantcare || return 0

    log_step "PlantCARE Promoter Analysis: $group"

    # Build --steps string from enabled sub-operations
    local -a step_list=()
    should_run "PlantCARE_post_processing"    && step_list+=("post-process")
    should_run "PlantCARE_matrix_generation"  && step_list+=("matrix")
    should_run "PlantCARE_heatmap_generation" && step_list+=("heatmap")
    should_run "PlantCARE_combined_heatmap"   && step_list+=("combined-heatmap")
    should_run "PlantCARE_groups_heatmap"     && step_list+=("groups-heatmap")
    local steps_str
    steps_str=$(IFS=','; echo "${step_list[*]:-post-process,matrix,heatmap}")

    local mod_dir="$MODULES/07_plantcare_analysis"
    local dpi
    dpi=$(get_toml plantcare dpi 2>/dev/null || echo "900")

    local genes funcs
    genes=$(get_toml plantcare selected_genes 2>/dev/null || true)
    funcs=$(get_toml plantcare selected_functions 2>/dev/null || true)

    # Color configuration
    local color_palette cell_border_color
    color_palette=$(get_toml plantcare colors palette 2>/dev/null || true)
    cell_border_color=$(get_toml plantcare colors cell_border_color 2>/dev/null || echo "white")

    # Figure configuration
    local row_font col_font label_font cell_font cell_size min_freq
    local column_rotation cell_border_width legend_height gene_label
    row_font=$(get_toml plantcare figures row_font 2>/dev/null || echo "12")
    col_font=$(get_toml plantcare figures col_font 2>/dev/null || echo "10")
    label_font=$(get_toml plantcare figures label_font 2>/dev/null || echo "14")
    cell_font=$(get_toml plantcare figures cell_font 2>/dev/null || echo "10")
    cell_size=$(get_toml plantcare figures cell_size 2>/dev/null || echo "1.0")
    min_freq=$(get_toml plantcare figures min_freq 2>/dev/null || echo "0")
    column_rotation=$(get_toml plantcare figures column_rotation 2>/dev/null || echo "22.5")
    cell_border_width=$(get_toml plantcare figures cell_border_width 2>/dev/null || echo "0.5")
    legend_height=$(get_toml plantcare figures legend_height 2>/dev/null || echo "4")
    gene_label=$(get_toml plantcare figures gene_label 2>/dev/null || echo "Genes of Interest")

    # Helper: run the PlantCARE pipeline for one raw-dir → genome outdir
    _run_plantcare_pipeline() {
        local raw="$1" out="$2"
        local -a cmd=(
            bash "$MODULES/07_plantcare_analysis/run_pipeline.sh"
            --raw-dir "$raw"
            --outdir "$out"
            --module-dir "$mod_dir"
            --steps "$steps_str"
            --dpi "$dpi"
            --threads "$CPU"
            --row-font "$row_font"
            --col-font "$col_font"
            --label-font "$label_font"
            --cell-font "$cell_font"
            --cell-size "$cell_size"
            --min-freq "$min_freq"
            --column-rotation "$column_rotation"
            --cell-border-color "$cell_border_color"
            --cell-border-width "$cell_border_width"
            --legend-height "$legend_height"
            --gene-label "$gene_label"
        )
        [[ "$OVERWRITE" == true ]] && cmd+=(--overwrite)
        [[ -n "$genes" ]] && cmd+=(--genes "$(echo "$genes" | tr '\n' ',' | sed 's/,$//')")
        [[ -n "$funcs" ]] && cmd+=(--functions "$(echo "$funcs" | tr '\n' ',' | sed 's/,$//')")
        [[ -n "$color_palette" ]] && cmd+=(--color-palette "$(echo "$color_palette" | tr '\n' ',' | sed 's/,$//')")
        "${cmd[@]}"
    }

    # Multi-genome mode: iterate genomes so outputs land under
    # III_RESULT/[gene_group]/07_PlantCARE_Analysis/[genome_name]/
    local -a gtf_labels
    mapfile -t gtf_labels < <(get_toml gene_sequence gtf_extraction labels 2>/dev/null || true)

    if [[ ${#gtf_labels[@]} -gt 0 ]]; then
        for genome_label in "${gtf_labels[@]}"; do
            local genome_dir="$pc_dir/$genome_label"

            # Collect all raw PlantCARE result dirs (per-family or single).
            # Using dirname of each found dir as its output dir means per-family
            # dirs (e.g. genome_dir/GRF/03_PlantCARE_Results) output to
            # genome_dir/GRF/, while a top-level dir outputs to genome_dir/.
            local -a all_raw_dirs=()
            while IFS= read -r raw; do
                find "$raw" -type f \( -name "*.tab" -o -name "*.tar.gz" \) \
                    2>/dev/null | grep -q . && all_raw_dirs+=("$raw")
            done < <(find "$genome_dir" -mindepth 1 -maxdepth 3 \
                -type d \( -name "$PLANTCARE_RAW_STAGE_DIR" -o -name "PlantCARE_Results" \) \
                2>/dev/null | sort)

            if [[ ${#all_raw_dirs[@]} -eq 0 ]]; then
                log_info "No PlantCARE raw results for $genome_label — skipping analysis"
                continue
            fi

            for raw_dir in "${all_raw_dirs[@]}"; do
                local out_dir
                out_dir=$(dirname "$raw_dir")
                log_info "PlantCARE raw directory: $raw_dir"
                _run_plantcare_pipeline "$raw_dir" "$out_dir"
            done
        done
    else
        # Monolithic fallback (GRF_GIF / PLA): single raw directory
        local raw_dir
        raw_dir=$(resolve_plantcare_raw_dir "$base_dir" "$pc_dir")
        if [[ -z "$raw_dir" ]]; then
            log_error "PlantCARE raw results directory not found for $group"
            return 1
        fi
        if ! find "$raw_dir" -type f \( -name "*.tab" -o -name "*.tar.gz" \) | grep -q .; then
            log_error "No PlantCARE .tab or .tar.gz files found in: $raw_dir"
            return 1
        fi
        log_info "PlantCARE raw directory: $raw_dir"
        _run_plantcare_pipeline "$raw_dir" "$pc_dir"
    fi

    log_step "PlantCARE analysis complete: $group"
}

# ---------- main loop --------------------------------------------------------

process_gene_group() {
    local group="$1"
    local base_dir pc_dir

    resolve_group_config "$group"
    base_dir="$PIPELINE_DIR/$(get_toml general base_dir)"
    pc_dir=$(resolve_plantcare_dir "$base_dir")
    local pc_subdir
    pc_subdir=$(get_toml plantcare subdir 2>/dev/null || true)
    [[ -n "$pc_subdir" ]] && pc_dir="$pc_dir/$pc_subdir"
    mkdir -p "$pc_dir"

    setup_logging

    load_operations_from_config
    log_info "Operations: ${OPERATIONS[*]}"

    run_gtf_extraction "$group" "$base_dir" "$pc_dir"
    run_sequence_extraction "$group" "$base_dir" "$pc_dir"
    run_plantcare_analysis "$group" "$base_dir" "$pc_dir"

    teardown_logging
}

for GENE_GROUP in "${GENE_GROUPS[@]}"; do
    process_gene_group "$GENE_GROUP"
done
