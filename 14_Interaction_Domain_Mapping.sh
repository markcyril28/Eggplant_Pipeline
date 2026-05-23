#! /bin/bash

# Rationale, paper excerpts, and the HAP2 prefusion-postfusion conformational
# pathway (Feng et al. 2022, Nat Commun) that motivates the three AF3
# stoichiometries below have been moved out of this script to keep the
# orchestrator readable. See:
#     14_Interaction_Domain_Mapping_NOTES.md  (project root, next to this script)
#
# AF3 inputs already on disk (one HAP2 model per stoichiometry, paired with
# SmelDMPv5_10.610):
#     1 HAP2 + 1 DMP  - prefusion monomer analogue (most consistent with the
#                       Wang et al. 2022 MbY2H / co-IP setup)
#     2 HAP2 + 1 DMP  - on-pathway intermediate (not directly supported by
#                       the paper; included for completeness)
#     3 HAP2 + 1 DMP  - postfusion class II trimer analogue (tests whether
#                       DMP can still engage the activated assembly)

# ============================================================================
# Program 14: In Silico Domain Mapping of the HAP2-DMP Interaction
# ============================================================================
# Computational analogue of the wet-lab MbY2H deletion-binding assay reported
# by Wang et al. (2022, PMC9659367). The paper introduced truncations into
# AtHAP2/GCS1 and asked which deletions abolished AtDMP9 binding. Here we
# replace yeast two-hybrid with a structural / energetic readout against
# SmelDMPv5_10.610, the eggplant DMP candidate identified in this thesis.
#
# Strategy:
#   1. Build truncation variants of SmelHAP2 (or AtHAP2) matching the paper,
#      and (in silico-only) parallel truncations of DMP probing the regions
#      the paper could not test in MbY2H (notably the DMP N-terminal
#      cytoplasmic specificity region established by Wang et al. 2022 Fig. 5
#      chimera rescue). HAP2 and DMP ladders combine into AF3 jobs via the
#      [pairing] section of the TOML; default = orthogonal ladders.
#   2. Predict each (HAP2_variant x DMP_variant) complex at 1:1, 2:1, 3:1
#      stoichiometry via AlphaFold3.
#   3. Score the predicted interface (BSA, contacts, PRODIGY DG_pred, ipTM).
#   4. Run short GROMACS MD (NVT + NPT + 50 ns production) to relax each
#      complex and confirm the interface is stable, not a docking artefact.
#   5. Decompose binding free energy via gmx_MMPBSA on the MD ensemble and
#      via FoldX AnalyseComplex on the predicted snapshot.
#   6. Computational alanine scan of the WT/WT complex interface (HAP2 WT x
#      DMP WT) via MutaTeX to flag hot-spot residues; cross-reference with
#      the regions deleted in the HAP2 and DMP ladders.
#   7. Rank pairs by predicted DG and interface preservation; emit a
#      heatmap and a summary TSV that mirrors the paper's Fig. 3A panel.
#
# Operations (gated via [domain_mapping].operations in the TOML):
#   prepare_variants    - build truncated HAP2 FASTAs from variant table
#   prepare_complexes   - AlphaFold3 multimer prediction slot prep + drop-zone
#                         verification (USER ACTION required for backend=manual;
#                         set [run].show_manual = true in the TOML to print
#                         the submission checklist)
#   interface_analysis  - residue contacts, BSA, ipTM, PRODIGY DG_pred
#   md_equilibration    - short GROMACS MD via the existing PPI stage-10 modules
#   binding_energy      - gmx_MMPBSA on the trajectory + FoldX AnalyseComplex
#   alanine_scan        - MutaTeX alanine scan of the WT-complex interface
#   comparative_report  - cross-variant heatmap, ranking, summary TSV / PDF
#
# Output layout under III_RESULT/{GROUP}/14_Domain_Mapping/:
#   stoichiometry_comparison/   (Experiment 1: WT x WT at 1:1, 2:1, 3:1)
#     01_Variants/HAP2/                       (WT HAP2 FASTA)
#     01_Variants/DMP/                        (WT DMP FASTA)
#     02_Complexes/{stoich}/{pair}/           (AF3 .cif / ranking_debug.json)
#     03_Interfaces/{stoich}/{pair}/          (interface .tsv, BSA, contacts)
#     04_MD/{stoich}/{pair}/                  (GROMACS .gro/.xtc/.edr)
#     05_BindingEnergy/{stoich}/{pair}/       (MM-PBSA, FoldX dG)
#     07_Summary/                             (ranked_stoichiometry.tsv + figure)
#   deletion_ladder/            (Experiment 2: HAP2 ladder x DMP ladder at basis stoich)
#     01_Variants/HAP2/                       (truncated HAP2 FASTAs)
#     01_Variants/DMP/                        (truncated DMP FASTAs)
#     02_Complexes/{stoich}/{pair}/           (AF3 .cif / ranking_debug.json)
#     03_Interfaces/{stoich}/{pair}/          (interface .tsv, BSA, contacts)
#     04_MD/{stoich}/{pair}/                  (GROMACS .gro/.xtc/.edr)
#     05_BindingEnergy/{stoich}/{pair}/       (MM-PBSA, FoldX dG)
#     06_AlanineScan/                         (per-residue DDG; WT/WT pair only)
#     07_Summary/                             (cross-pair ranking, figures)
#
# {pair} = "{hap2_variant_name}__{dmp_variant_name}", e.g. "WT__WT",
# "delC_596_705__WT", "WT__delN_1_64". Built from [hap2_variants] x [dmp_variants]
# according to [pairing].mode (see TOML).
# Active experiments are gated by `active_comparison = [...]` in the TOML.
#
# References:
#   Tooling DOIs are listed at the bottom of 14_Interaction_Domain_MappingCONFIG.toml.
#   All biological background and paper-anchored design decisions (Wang 2022
#   deletion-binding paper; Feng 2022 HAP2 prefusion-monomer structure;
#   Cyprys 2019 DMP8/9 discovery; Shiba 2023 GCS1/HAP2 review; EC1 trigger
#   via Sprunck 2012; Zhang 2021 species-specific gamete recognition; DMP
#   haploid-induction lineage from Zhong 2019 onward) are documented with
#   full, citation-checked references in:
#       14_Interaction_Domain_Mapping_NOTES.md   (project root, Section 10)
#   Keep biological citations in sync there, not here.
# ============================================================================

set -euo pipefail

# ── Path setup ──────────────────────────────────────────────────────────────
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"
PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

SHARED_CONFIG="$PIPELINE_DIR/14_Interaction_Domain_MappingCONFIG.toml"
get_toml_shared() { python3 "$TOML_PARSER" "$SHARED_CONFIG" "$@"; }

# ── Run-mode toggles (all sourced from [run] in 14_Interaction_Domain_MappingCONFIG.toml) ──
# Edit the TOML, not this script. No CLI flags are accepted.
DRY_RUN=$(get_toml_shared run dry_run 2>/dev/null || echo "false")
SHOW_MANUAL=$(get_toml_shared run show_manual 2>/dev/null || echo "false")
mapfile -t ONLY_OPS          < <(get_toml_shared run only_ops 2>/dev/null)
mapfile -t ONLY_VARIANTS     < <(get_toml_shared run only_variants 2>/dev/null)
mapfile -t ONLY_DMP_VARIANTS < <(get_toml_shared run only_dmp_variants 2>/dev/null)

# ── Manual steps panel ──────────────────────────────────────────────────────
# Steps the orchestrator cannot perform on its own. AlphaFold3 in particular
# has no free programmatic API (as of 2026-05); each prediction must go
# through alphafoldserver.com by hand. This panel prints what to do and
# where to put the downloads so the rest of the pipeline can pick them up.
print_manual_steps() {
    cat <<'MANUAL'
============================================================================
USER ACTION REQUIRED - steps only you can do
============================================================================

  [M0] Scope reminder: this pipeline tests DMP-HAP2 *binding*, not the
       upstream EC1 signalling that triggers it.
       In vivo, gamete fusion is a three-stage signalling cascade
       (Sprunck 2012 Science 338:1093; Wang 2022 PMC9659367; Shiba 2023
       PMC9953686 - see Section 6 of 14_Interaction_Domain_Mapping_NOTES.md):
         (i)  EGG SIDE: EC1.1-EC1.5 peptides stored in egg cytoplasmic
              vesicles; released by regulated exocytosis upon sperm
              arrival.
         (ii) SPERM SIDE: an as-yet-unidentified receptor on the sperm
              PM senses EC1; this somehow activates DMP8/9. Wang 2022 SI
              Fig. S4 showed DMP9 does NOT bind EC1 directly, so the
              sensor is not DMP8/9 itself.
         (iii)DOWNSTREAM: DMP8/9 and HAP2 co-translocate from internal
              vesicles to the sperm PM. This is the binding step Stage
              14 maps.
       The AF3 complexes submitted from this pipeline are
       HAP2 + DMP (no EC1), corresponding to stage (iii) only. EC1 is
       neither modelled nor required for the binding measurement -
       Wang 2022 captured DMP-HAP2 binding in MbY2H and pull-down in
       the absence of EC1 too. Document this scope limit in your
       manuscript Discussion: the in silico map characterises
       DMP-HAP2 *interaction interface*, not the *trigger* that
       launches it.

  [M1] Decide HAP2 source sequence.
       The paper uses AtHAP2 (At4g11720, 705 aa). For an eggplant-native
       readout, use the SmelHAP2 ortholog you identified in Chapter IV.
       Whichever you choose, place its FASTA at the path declared in TOML
       [inputs].hap2_fasta. Verify residue numbering against UniProt before
       editing the [variants] deletion ranges - the paper's coordinates
       (596-705, 25-530, etc.) are AtHAP2 numbering and may shift in
       SmelHAP2.

  [M2] Map AtHAP2 deletion boundaries onto your HAP2 sequence.
       Run MAFFT --add or a pairwise Needleman-Wunsch alignment of AtHAP2
       to your HAP2, then transcribe the AtHAP2 25/530/596/609/705 cut
       points onto the aligned columns. Update [variants].deletions
       accordingly. The orchestrator will refuse to truncate past the
       sequence length, but it cannot tell you whether the cut falls inside
       a fold-critical helix - that is your call.

  [M3] Submit each variant + DMP complex to AlphaFold3.
       Free non-commercial AF3 inference is web-only at
           https://alphafoldserver.com/
       (20 jobs/day; multimer jobs supported.) For each row in
       [variants].names x [stoichiometry].chain_counts:
         - Open a new job.
         - Add the HAP2 variant sequence (1 copy if monomeric, 2 if
           dimeric, 3 if trimeric).
         - Add the SmelDMPv5_10.610 sequence (1 copy).
         - Submit, wait, download the ZIP.
         - Unzip into
             III_RESULT/{GROUP}/14_Domain_Mapping/02_Complexes/{stoich}/{variant}/
           preserving the AF3 default filenames (fold_*_model_0.cif,
           ranking_scores.json, summary_confidences_*.json).
       If you have a commercial AF3 license you can swap this for a local
       run by pointing [prepare_complexes].backend = "local" and providing
       the AF3 binary path - but the default backend is "manual".

  [M4] Install FoldX (free academic license).
       Request a license at https://foldxsuite.crg.eu/, download the
       binary, and place it at the path in [tools].foldx_binary. The
       MutaTeX dependency is already covered by setup_gromacs_and_mutatex.sh
       at the repo root; FoldX is not.

  [M5] Confirm GROMACS GPU build.
       The MD stage assumes CUDA-compiled GROMACS (see CLAUDE.md GROMACS
       notes). Run `gmx --version | grep -i cuda` and confirm "GPU support:
       CUDA". OpenCL builds on WSL2 will silently hang during NPT.

  [M6] Decide whether to keep the membrane.
       HAP2 is type-1 membrane; the paper's MbY2H captured the full
       trans-membrane context. If you want to honour that, run MD in a
       POPC/POPE bilayer via CHARMM-GUI Membrane Builder
       (https://charmm-gui.org/?doc=input/membrane.bilayer), download the
       prepared GROMACS box, and set [md_equilibration].membrane = true.
       The default in this pipeline is soluble (no membrane) for
       tractability; document this caveat in your manuscript.

  [M7] Sanity-check interface hotspots against the paper.
       After alanine_scan runs, manually compare the top DDG residues to
       the regions Wang et al. deleted. The pipeline cannot tell you
       whether your prediction "agrees with the paper" - that requires
       interpreting both the residue-level data and the published figures
       together.

============================================================================
MANUAL
}

if [[ "$SHOW_MANUAL" == "true" ]]; then
    print_manual_steps
    exit 0
fi

# ── Gene-group iteration ────────────────────────────────────────────────────
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups is empty in 14_Interaction_Domain_MappingCONFIG.toml" >&2
    exit 1
fi

wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done
}

op_enabled() {
    local op="$1"
    # [run].only_ops in the TOML overrides [domain_mapping].operations when set
    if (( ${#ONLY_OPS[@]} > 0 )); then
        for o in "${ONLY_OPS[@]}"; do [[ "$o" == "$op" ]] && return 0; done
        return 1
    fi
    for o in "${OPERATIONS[@]}"; do [[ "$o" == "$op" ]] && return 0; done
    return 1
}

variant_enabled() {
    local v="$1"
    (( ${#ONLY_VARIANTS[@]} == 0 )) && return 0
    for x in "${ONLY_VARIANTS[@]}"; do [[ "$x" == "$v" ]] && return 0; done
    return 1
}

dmp_variant_enabled() {
    local v="$1"
    (( ${#ONLY_DMP_VARIANTS[@]} == 0 )) && return 0
    for x in "${ONLY_DMP_VARIANTS[@]}"; do [[ "$x" == "$v" ]] && return 0; done
    return 1
}

run_or_echo() {
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] $*"
    else
        log_info "Running: $*"
        "$@"
    fi
}

# Dynamic per-tool conda-env switching. Each Stage 14 operation lives in
# its own conda env (Biopython / AlphaFold3 / PRODIGY / GROMACS / gmx_MMPBSA
# / MutaTeX / plotting); the env name is taken from [conda_envs] in the
# TOML and resolved by `env_for <op>` below. `conda_run_in <env> <cmd...>`
# wraps the command with `conda run -n <env> --no-capture-output` and
# honours $DRY_RUN. An empty env name skips the wrapping (runs in the
# orchestrator's inherited env), which is the right behaviour for tools
# that are static binaries (FoldX) or whose wrapper script handles its
# own activation.
conda_run_in() {
    local env="$1"; shift
    if [[ "$DRY_RUN" == "true" ]]; then
        if [[ -z "$env" ]]; then
            log_info "[DRY-RUN] $*"
        else
            log_info "[DRY-RUN] conda run -n $env --no-capture-output $*"
        fi
        return 0
    fi
    if [[ -z "$env" ]]; then
        log_info "Running (no conda wrap): $*"
        "$@"
    else
        log_info "Running (env=$env): $*"
        conda run -n "$env" --no-capture-output "$@"
    fi
}

# Trim leading/trailing whitespace from a string. Uses bash parameter
# expansion so it's quick and dependency-free.
trim_ws() {
    local s="$1"
    s="${s#"${s%%[![:space:]]*}"}"
    s="${s%"${s##*[![:space:]]}"}"
    printf '%s' "$s"
}

# Split a "name | description" row from [variants].rows / [dmp_variants].rows
# into the global _ROW_NAME / _ROW_DESC. Empty / whitespace-only rows
# return empty _ROW_NAME so the caller can skip.
parse_variant_row() {
    local row="$1"
    _ROW_NAME=""
    _ROW_DESC=""
    [[ -z "$(trim_ws "$row")" ]] && return
    local rest
    _ROW_NAME="${row%%|*}"
    rest="${row#*|}"
    if [[ "$rest" == "$row" ]]; then
        # No '|' in the row; treat the whole row as the name, blank desc.
        rest=""
    fi
    _ROW_NAME="$(trim_ws "$_ROW_NAME")"
    _ROW_DESC="$(trim_ws "$rest")"
}

# Resolve the conda env for a given operation key. Looks up
# [conda_envs].<op> first; falls back to legacy [tools].prodigy_conda_env
# / gmx_mmpbsa_env where appropriate so older configs still work. Returns
# empty string when neither is set.
env_for() {
    local op="$1" val=""
    val=$(get_toml conda_envs "$op" 2>/dev/null || echo "")
    if [[ -z "$val" ]]; then
        case "$op" in
            interface_analysis|prodigy_dg|comparative_report)
                val=$(get_toml tools prodigy_conda_env 2>/dev/null || echo "")
                ;;
            mmpbsa)
                val=$(get_toml tools gmx_mmpbsa_env 2>/dev/null || echo "")
                ;;
        esac
    fi
    printf '%s' "$val"
}

TMP_CONFIG_FILES=()
cleanup_tmp_configs() {
    local cfg
    for cfg in "${TMP_CONFIG_FILES[@]:-}"; do
        [[ -n "$cfg" && -f "$cfg" ]] && rm -f "$cfg"
    done
}
trap 'cleanup_tmp_configs; safe_teardown_logging; true' EXIT

# ============================================================================
# Per-group execution
# ============================================================================
for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# Config resolution: shared root + optional per-group override (deep-merged
# via merge_toml.py - same pattern as Stages 01/05/06/07). The earlier `cat`
# concat broke as soon as the override redeclared a section already present
# in the shared file because tomllib refuses duplicate tables.
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
PER_GROUP_CFG="$CONFIG_DIR/14_Interaction_Domain_Mapping_${GENE_GROUP}.toml"
if [[ -f "$PER_GROUP_CFG" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_domain_mapping_cfg_XXXXXX.toml")
    python3 "$MODULES/utils/merge_toml.py" "$SHARED_CONFIG" "$PER_GROUP_CFG" > "$CONFIG_FILE"
    TMP_CONFIG_FILES+=("$CONFIG_FILE")
else
    CONFIG_FILE="$SHARED_CONFIG"
fi

# Compute profile (Local vs HPC); honours host nproc as an upper bound
MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
HOST_CPU=$(nproc)
CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || echo "$HOST_CPU")
(( CPU > HOST_CPU )) && CPU=$HOST_CPU
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || echo "$CPU")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
(( MAX_PARALLEL > CPU )) && MAX_PARALLEL=$CPU
OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "false")

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir 2>/dev/null || echo "III_RESULT/${GENE_GROUP}")"
OUT_DIR="$BASE_DIR/14_Domain_Mapping"
mkdir -p "$OUT_DIR"

setup_logging

ENABLED=$(get_toml domain_mapping enabled 2>/dev/null || echo "true")
if [[ "$ENABLED" != "true" && "$ENABLED" != "True" ]]; then
    log_info "Domain mapping disabled for $GENE_GROUP. Skipping."
    teardown_logging
    continue
fi

# active_comparison gates which experiments run (subset of: stoichiometry_comparison | deletion_ladder).
# Each experiment gets its own output subtree under $OUT_DIR/$EXPERIMENT/ and its own operations
# list under [domain_mapping.operations].$EXPERIMENT.
mapfile -t ACTIVE_COMPARISON < <(get_toml_shared active_comparison 2>/dev/null)
if [[ ${#ACTIVE_COMPARISON[@]} -eq 0 ]]; then
    log_warn "active_comparison is empty in 14_Interaction_Domain_MappingCONFIG.toml; defaulting to both experiments."
    ACTIVE_COMPARISON=(stoichiometry_comparison deletion_ladder)
fi

# Inputs
HAP2_FASTA=$(get_toml inputs hap2_fasta 2>/dev/null)
DMP_FASTA=$(get_toml inputs dmp_fasta 2>/dev/null)
HAP2_REF_PDB=$(get_toml inputs hap2_reference_pdb 2>/dev/null || echo "")
[[ "$HAP2_FASTA" != /* && -n "$HAP2_FASTA" ]] && HAP2_FASTA="$PIPELINE_DIR/$HAP2_FASTA"
[[ "$DMP_FASTA"  != /* && -n "$DMP_FASTA"  ]] && DMP_FASTA="$PIPELINE_DIR/$DMP_FASTA"
if [[ ! -f "$HAP2_FASTA" || ! -f "$DMP_FASTA" ]]; then
    log_error "Required input FASTA missing. HAP2='$HAP2_FASTA' DMP='$DMP_FASTA'"
    log_error "Set [run].show_manual = true in $SHARED_CONFIG and re-run to print step [M1] (source-sequence selection)."
    teardown_logging
    continue
fi

# ── HAP2 variant table ──────────────────────────────────────────────────────
# Single-list layout: [hap2_variants].rows holds "name | description" rows;
# residue ranges live in [hap2_variants.coords.<at|sm>].<name>.
# coords_column selects which per-species coord table the orchestrator
# reads. Anything in coords with value "TBD" (or empty) is treated as
# "not yet remapped" and that variant is skipped with a warning.
HAP2_COORDS_COL=$(get_toml hap2_variants coords_column 2>/dev/null || echo "at")
case "$HAP2_COORDS_COL" in
    at|sm) ;;
    *)
        log_error "[hap2_variants].coords_column = '$HAP2_COORDS_COL'; must be 'at' or 'sm'"
        teardown_logging
        continue
        ;;
esac

mapfile -t _VARIANT_ROWS < <(get_toml hap2_variants rows 2>/dev/null)
if [[ ${#_VARIANT_ROWS[@]} -eq 0 ]]; then
    log_error "[hap2_variants].rows is empty - nothing to map."
    teardown_logging
    continue
fi

VARIANT_NAMES=()
VARIANT_DESCRIPTIONS=()
VARIANT_DELETIONS=()
for _row in "${_VARIANT_ROWS[@]}"; do
    parse_variant_row "$_row"
    [[ -z "$_ROW_NAME" ]] && continue
    if [[ "$_ROW_NAME" == "WT" ]]; then
        _coord=""
    else
        _coord=$(get_toml hap2_variants coords "$HAP2_COORDS_COL" "$_ROW_NAME" 2>/dev/null || echo "")
        if [[ -z "$_coord" || "$_coord" == "TBD" ]]; then
            log_warn "[hap2_variants] '$_ROW_NAME': no $HAP2_COORDS_COL coords (value '$_coord'); skipping. Add a range under [hap2_variants.coords.$HAP2_COORDS_COL].$_ROW_NAME to enable."
            continue
        fi
    fi
    VARIANT_NAMES+=("$_ROW_NAME")
    VARIANT_DESCRIPTIONS+=("$_ROW_DESC")
    VARIANT_DELETIONS+=("$_coord")
done
unset _row _coord

if [[ ${#VARIANT_NAMES[@]} -eq 0 ]]; then
    log_error "[hap2_variants]: all rows skipped (no usable $HAP2_COORDS_COL coords). Fix [hap2_variants.coords.$HAP2_COORDS_COL] and re-run."
    teardown_logging
    continue
fi

# ── DMP variant table ───────────────────────────────────────────────────────
# Symmetric peer of the HAP2 block above: same layout, same coords_column
# selector, same at/sm sub-tables. Default = "sm" since the shared
# [inputs].dmp_fasta is SmelDMP; the DMP_x_AtHAP2 per-run override flips
# it to "at" alongside its dmp_fasta override.
DMP_COORDS_COL=$(get_toml dmp_variants coords_column 2>/dev/null || echo "sm")
case "$DMP_COORDS_COL" in
    at|sm) ;;
    *)
        log_error "[dmp_variants].coords_column = '$DMP_COORDS_COL'; must be 'at' or 'sm'"
        teardown_logging
        continue
        ;;
esac

mapfile -t _DMP_VARIANT_ROWS < <(get_toml dmp_variants rows 2>/dev/null)
DMP_VARIANT_NAMES=()
DMP_VARIANT_DESCRIPTIONS=()
DMP_VARIANT_DELETIONS=()
if [[ ${#_DMP_VARIANT_ROWS[@]} -eq 0 ]]; then
    log_info "[dmp_variants].rows not declared - using single WT DMP partner from [inputs].dmp_fasta"
    DMP_VARIANT_NAMES=("WT")
    DMP_VARIANT_DESCRIPTIONS=("Full-length DMP (inferred default; [dmp_variants].rows not declared)")
    DMP_VARIANT_DELETIONS=("")
else
    for _row in "${_DMP_VARIANT_ROWS[@]}"; do
        parse_variant_row "$_row"
        [[ -z "$_ROW_NAME" ]] && continue
        if [[ "$_ROW_NAME" == "WT" ]]; then
            _coord=""
        else
            _coord=$(get_toml dmp_variants coords "$DMP_COORDS_COL" "$_ROW_NAME" 2>/dev/null || echo "")
            if [[ -z "$_coord" || "$_coord" == "TBD" ]]; then
                log_warn "[dmp_variants] '$_ROW_NAME': no $DMP_COORDS_COL coords (value '$_coord'); skipping. Add a range under [dmp_variants.coords.$DMP_COORDS_COL].$_ROW_NAME to enable."
                continue
            fi
        fi
        DMP_VARIANT_NAMES+=("$_ROW_NAME")
        DMP_VARIANT_DESCRIPTIONS+=("$_ROW_DESC")
        DMP_VARIANT_DELETIONS+=("$_coord")
    done
    unset _row _coord
    if [[ ${#DMP_VARIANT_NAMES[@]} -eq 0 ]]; then
        log_warn "[dmp_variants]: all rows skipped (no usable $DMP_COORDS_COL coords) - falling back to WT partner."
        DMP_VARIANT_NAMES=("WT")
        DMP_VARIANT_DESCRIPTIONS=("Full-length DMP (fallback; all rows skipped)")
        DMP_VARIANT_DELETIONS=("")
    fi
fi

# build_pairs builds PAIR_HAP2/PAIR_DMP/PAIR_LABEL from the active HAP2 +
# DMP variant lists (EXP_HAP2_NAMES / EXP_DMP_NAMES) according to mode.
# Called per-experiment so the pair list reflects the active ladder.
#   PAIR_HAP2[k]  = HAP2 variant name at pair index k
#   PAIR_DMP[k]   = DMP variant name at pair index k
#   PAIR_LABEL[k] = "{hap2}__{dmp}" composite label used in all output paths
PAIRING_MODE=$(get_toml pairing mode 2>/dev/null || echo "orthogonal")
build_pairs() {
    local mode="$1"
    local i j h d
    PAIR_HAP2=()
    PAIR_DMP=()
    PAIR_LABEL=()
    case "$mode" in
        orthogonal)
            # HAP2 ladder x WT DMP (DMP index 0 must be the WT row)
            for i in "${!EXP_HAP2_NAMES[@]}"; do
                PAIR_HAP2+=("${EXP_HAP2_NAMES[$i]}")
                PAIR_DMP+=("${EXP_DMP_NAMES[0]}")
                PAIR_LABEL+=("${EXP_HAP2_NAMES[$i]}__${EXP_DMP_NAMES[0]}")
            done
            # WT HAP2 x DMP ladder (skip j=0 to avoid duplicating WT__WT)
            for j in "${!EXP_DMP_NAMES[@]}"; do
                (( j == 0 )) && continue
                PAIR_HAP2+=("${EXP_HAP2_NAMES[0]}")
                PAIR_DMP+=("${EXP_DMP_NAMES[$j]}")
                PAIR_LABEL+=("${EXP_HAP2_NAMES[0]}__${EXP_DMP_NAMES[$j]}")
            done
            ;;
        matrix)
            for i in "${!EXP_HAP2_NAMES[@]}"; do
                for j in "${!EXP_DMP_NAMES[@]}"; do
                    PAIR_HAP2+=("${EXP_HAP2_NAMES[$i]}")
                    PAIR_DMP+=("${EXP_DMP_NAMES[$j]}")
                    PAIR_LABEL+=("${EXP_HAP2_NAMES[$i]}__${EXP_DMP_NAMES[$j]}")
                done
            done
            ;;
        pairwise)
            local nh=${#EXP_HAP2_NAMES[@]} nd=${#EXP_DMP_NAMES[@]} n
            (( nh > nd )) && n=$nh || n=$nd
            for ((i=0; i<n; i++)); do
                (( i < nh )) && h="${EXP_HAP2_NAMES[$i]}"      || h="${EXP_HAP2_NAMES[0]}"
                (( i < nd )) && d="${EXP_DMP_NAMES[$i]}" || d="${EXP_DMP_NAMES[0]}"
                PAIR_HAP2+=("$h")
                PAIR_DMP+=("$d")
                PAIR_LABEL+=("${h}__${d}")
            done
            ;;
        *)
            log_error "Unknown [pairing].mode = '$mode'. Valid: orthogonal | matrix | pairwise"
            return 1
            ;;
    esac
}

# A pair passes the user filter iff its HAP2 side passes only_variants AND
# its DMP side passes only_dmp_variants.
pair_enabled() {
    local hap2="$1" dmp="$2"
    variant_enabled "$hap2" && dmp_variant_enabled "$dmp"
}

# Tooling
PREPARE_BACKEND=$(get_toml prepare_complexes backend 2>/dev/null || echo "manual")
MD_BACKEND=$(get_toml md_equilibration backend 2>/dev/null || echo "gromacs")
DG_BACKENDS=$(get_toml binding_energy backends 2>/dev/null || echo "mmpbsa foldx prodigy")
FOLDX_BIN=$(get_toml tools foldx_binary 2>/dev/null || echo "")
GMX_BIN=$(get_toml tools gromacs_binary 2>/dev/null || echo "gmx")
PRODIGY_ENV=$(get_toml tools prodigy_conda_env 2>/dev/null || echo "egg")

# Per-operation conda envs resolved once per gene group. Each one is
# applied as `conda run -n <env>` around the matching invocation below.
# Empty values fall through to running in the orchestrator's inherited
# env (correct for static-binary tools like FoldX). See [conda_envs] in
# 14_Interaction_Domain_MappingCONFIG.toml for the per-op mapping.
ENV_PREPARE=$(env_for prepare_variants)
ENV_PREPARE_COMPLEXES=$(env_for prepare_complexes)
ENV_IFACE=$(env_for interface_analysis)
ENV_MD=$(env_for md_equilibration)
ENV_MMPBSA=$(env_for mmpbsa)
ENV_FOLDX=$(env_for foldx_dg)
ENV_PRODIGY_DG=$(env_for prodigy_dg)
ENV_ASCAN=$(env_for alanine_scan)
ENV_REPORT=$(env_for comparative_report)

log_step "Stage 14 (Domain Mapping): $GENE_GROUP  (active: ${ACTIVE_COMPARISON[*]})"
log_info "Complex prep backend: $PREPARE_BACKEND   MD backend: $MD_BACKEND"
log_info "Compute: host=${HOST_CPU}c  CPU=$CPU  MAX_PARALLEL=$MAX_PARALLEL"
log_info "Conda envs (per op):"
log_info "  prepare_variants   = '${ENV_PREPARE}'"
log_info "  prepare_complexes  = '${ENV_PREPARE_COMPLEXES}'   (used only when backend=local)"
log_info "  interface_analysis = '${ENV_IFACE}'"
log_info "  md_equilibration   = '${ENV_MD}'"
log_info "  mmpbsa             = '${ENV_MMPBSA}'"
log_info "  foldx_dg           = '${ENV_FOLDX}'   (empty = static binary, no env)"
log_info "  prodigy_dg         = '${ENV_PRODIGY_DG}'"
log_info "  alanine_scan       = '${ENV_ASCAN}'"
log_info "  comparative_report = '${ENV_REPORT}'"
log_info "Output root:   $OUT_DIR"

# ============================================================================
# Per-experiment execution loop
# ============================================================================
# Two experiments share the same per-group setup but write into separate
# output subtrees with their own pair/stoichiometry config and their own
# operations list. See [domain_mapping.operations] in the TOML.
for EXPERIMENT in "${ACTIVE_COMPARISON[@]}"; do
    case "$EXPERIMENT" in
        stoichiometry_comparison|deletion_ladder) ;;
        *)
            log_warn "Unknown experiment '$EXPERIMENT' in active_comparison; skipping."
            continue
            ;;
    esac

    # Per-experiment output dirs (the two top-level folders the user expects)
    EXP_OUT_DIR="$OUT_DIR/$EXPERIMENT"
    HAP2_VARIANTS_DIR="$EXP_OUT_DIR/01_Variants/HAP2"
    DMP_VARIANTS_DIR="$EXP_OUT_DIR/01_Variants/DMP"
    COMPLEX_DIR="$EXP_OUT_DIR/02_Complexes"
    IFACE_DIR="$EXP_OUT_DIR/03_Interfaces"
    MD_DIR="$EXP_OUT_DIR/04_MD"
    DG_DIR="$EXP_OUT_DIR/05_BindingEnergy"
    ASCAN_DIR="$EXP_OUT_DIR/06_AlanineScan"
    REPORT_DIR="$EXP_OUT_DIR/07_Summary"
    mkdir -p "$HAP2_VARIANTS_DIR" "$DMP_VARIANTS_DIR" \
             "$COMPLEX_DIR" "$IFACE_DIR" "$MD_DIR" "$DG_DIR" \
             "$ASCAN_DIR" "$REPORT_DIR"

    # Per-experiment operations list. The TOML structure is
    # [domain_mapping.operations].<experiment> = [ ... ], so we query that
    # nested array directly. Empty / missing list => no ops run.
    mapfile -t OPERATIONS < <(get_toml domain_mapping operations "$EXPERIMENT" 2>/dev/null)
    if [[ ${#OPERATIONS[@]} -eq 0 ]]; then
        log_warn "  [$EXPERIMENT] no operations declared under [domain_mapping.operations].$EXPERIMENT; skipping."
        continue
    fi

    # Per-experiment variant + pair + stoichiometry tables.
    # - stoichiometry_comparison: a single HAP2/DMP pair (default WT/WT) at
    #   3 stoichiometries (1, 2, 3 HAP2 copies per [stoichiometry_comparison]).
    # - deletion_ladder: full HAP2 x DMP variant matrix combined per
    #   [pairing].mode, evaluated at a SINGLE stoichiometry = [stoichiometry].basis.
    EXP_HAP2_NAMES=()
    EXP_HAP2_DESCS=()
    EXP_HAP2_DELETIONS=()
    EXP_DMP_NAMES=()
    EXP_DMP_DESCS=()
    EXP_DMP_DELETIONS=()
    PAIR_HAP2=()
    PAIR_DMP=()
    PAIR_LABEL=()
    STOICH_LABELS=()
    STOICH_COUNTS=()

    if [[ "$EXPERIMENT" == "stoichiometry_comparison" ]]; then
        SC_HAP2=$(get_toml stoichiometry_comparison hap2_variant 2>/dev/null || echo "WT")
        SC_DMP=$(get_toml stoichiometry_comparison dmp_variant 2>/dev/null || echo "WT")
        # Resolve the named HAP2 + DMP variants against the full ladders so
        # we inherit description + deletion range.
        for i in "${!VARIANT_NAMES[@]}"; do
            if [[ "${VARIANT_NAMES[$i]}" == "$SC_HAP2" ]]; then
                EXP_HAP2_NAMES+=("${VARIANT_NAMES[$i]}")
                EXP_HAP2_DESCS+=("${VARIANT_DESCRIPTIONS[$i]}")
                EXP_HAP2_DELETIONS+=("${VARIANT_DELETIONS[$i]}")
                break
            fi
        done
        for j in "${!DMP_VARIANT_NAMES[@]}"; do
            if [[ "${DMP_VARIANT_NAMES[$j]}" == "$SC_DMP" ]]; then
                EXP_DMP_NAMES+=("${DMP_VARIANT_NAMES[$j]}")
                EXP_DMP_DESCS+=("${DMP_VARIANT_DESCRIPTIONS[$j]}")
                EXP_DMP_DELETIONS+=("${DMP_VARIANT_DELETIONS[$j]}")
                break
            fi
        done
        if [[ ${#EXP_HAP2_NAMES[@]} -eq 0 ]]; then
            log_error "  [stoichiometry_comparison].hap2_variant='$SC_HAP2' not found in [hap2_variants].rows; skipping experiment."
            continue
        fi
        if [[ ${#EXP_DMP_NAMES[@]} -eq 0 ]]; then
            log_error "  [stoichiometry_comparison].dmp_variant='$SC_DMP' not found in [dmp_variants].rows; skipping experiment."
            continue
        fi
        PAIR_HAP2=("$SC_HAP2")
        PAIR_DMP=("$SC_DMP")
        PAIR_LABEL=("${SC_HAP2}__${SC_DMP}")
        mapfile -t STOICH_LABELS < <(get_toml stoichiometry_comparison labels 2>/dev/null \
            || printf '%s\n' monomeric dimeric postfusion_like)
        mapfile -t STOICH_COUNTS < <(get_toml stoichiometry_comparison chain_counts 2>/dev/null \
            || printf '%s\n' 1 2 3)
    else
        # deletion_ladder uses the full HAP2 + DMP ladders
        EXP_HAP2_NAMES=("${VARIANT_NAMES[@]}")
        EXP_HAP2_DESCS=("${VARIANT_DESCRIPTIONS[@]}")
        EXP_HAP2_DELETIONS=("${VARIANT_DELETIONS[@]}")
        EXP_DMP_NAMES=("${DMP_VARIANT_NAMES[@]}")
        EXP_DMP_DESCS=("${DMP_VARIANT_DESCRIPTIONS[@]}")
        EXP_DMP_DELETIONS=("${DMP_VARIANT_DELETIONS[@]}")
        if ! build_pairs "$PAIRING_MODE"; then
            log_error "  [deletion_ladder] build_pairs failed for mode='$PAIRING_MODE'; skipping experiment."
            continue
        fi
        BASIS=$(get_toml stoichiometry basis 2>/dev/null || echo "monomeric")
        BASIS_COPIES=$(get_toml stoichiometry hap2_copies 2>/dev/null || echo "1")
        STOICH_LABELS=("$BASIS")
        STOICH_COUNTS=("$BASIS_COPIES")
    fi

    log_step "  Experiment '$EXPERIMENT' -> $EXP_OUT_DIR"
    log_info "    Operations:    ${OPERATIONS[*]}"
    log_info "    HAP2 variants: ${EXP_HAP2_NAMES[*]}"
    log_info "    DMP variants:  ${EXP_DMP_NAMES[*]}"
    log_info "    Pairs (${#PAIR_LABEL[@]}): ${PAIR_LABEL[*]}"
    log_info "    Stoichiometry: ${STOICH_LABELS[*]} (${STOICH_COUNTS[*]} HAP2 copies)"

# ============================================================================
# Operation 1: Build truncated HAP2 and DMP FASTAs
# ============================================================================
if op_enabled "prepare_variants"; then
    log_step "  Op 1/7: prepare_variants  (HAP2 + DMP ladders)"
    GEN_SCRIPT="$MODULES/14_special_pipeline/generate_variants.py"

    log_info "    HAP2 ladder -> $HAP2_VARIANTS_DIR"
    for i in "${!EXP_HAP2_NAMES[@]}"; do
        VNAME="${EXP_HAP2_NAMES[$i]}"
        variant_enabled "$VNAME" || continue
        VDESC="${EXP_HAP2_DESCS[$i]:-no description}"
        VDEL="${EXP_HAP2_DELETIONS[$i]:-}"   # empty = WT
        OUT_FA="$HAP2_VARIANTS_DIR/${VNAME}.fasta"
        if [[ -f "$OUT_FA" && "$OVERWRITE" != "true" ]]; then
            log_info "      [$VNAME] exists, skip (OVERWRITE=false)"
            continue
        fi
        log_info "      [$VNAME] $VDESC  (delete='$VDEL')"
        conda_run_in "$ENV_PREPARE" python3 "$GEN_SCRIPT" \
            --input "$HAP2_FASTA" \
            --name "$VNAME" \
            --description "$VDESC" \
            --deletions "$VDEL" \
            --output "$OUT_FA"
    done

    log_info "    DMP ladder  -> $DMP_VARIANTS_DIR"
    for j in "${!EXP_DMP_NAMES[@]}"; do
        DNAME="${EXP_DMP_NAMES[$j]}"
        dmp_variant_enabled "$DNAME" || continue
        DDESC="${EXP_DMP_DESCS[$j]:-no description}"
        DDEL="${EXP_DMP_DELETIONS[$j]:-}"
        OUT_FA="$DMP_VARIANTS_DIR/${DNAME}.fasta"
        if [[ -f "$OUT_FA" && "$OVERWRITE" != "true" ]]; then
            log_info "      [$DNAME] exists, skip (OVERWRITE=false)"
            continue
        fi
        log_info "      [$DNAME] $DDESC  (delete='$DDEL')"
        conda_run_in "$ENV_PREPARE" python3 "$GEN_SCRIPT" \
            --input "$DMP_FASTA" \
            --name "$DNAME" \
            --description "$DDESC" \
            --deletions "$DDEL" \
            --output "$OUT_FA"
    done
fi

# ============================================================================
# Operation 2: AlphaFold3 multimer prediction
# ============================================================================
# `backend = "manual"` is the realistic default - AF3 has no free programmatic
# inference API as of 2026-05. The orchestrator stages the work, prints the
# submission checklist, and confirms each download exists before moving on.
#
# Iterates over the HAP2 x DMP pair matrix built earlier (PAIR_HAP2/PAIR_DMP/
# PAIR_LABEL). Each (pair, stoich) cell becomes one AF3 job and one output
# directory tagged with the composite "{hap2}__{dmp}" pair label.
if op_enabled "prepare_complexes"; then
    log_step "  Op 2/7: prepare_complexes  (backend=$PREPARE_BACKEND, ${#PAIR_LABEL[@]} pairs)"
    MANIFEST="$COMPLEX_DIR/_submission_manifest.tsv"
    DROP_README="$COMPLEX_DIR/_HOW_TO_DROP_AF3_DOWNLOADS.txt"
    : > "$MANIFEST"
    printf 'pair\thap2_variant\tdmp_variant\tstoich\tn_hap2\tn_dmp\thap2_fasta\tdmp_fasta\texpected_output\n' >> "$MANIFEST"

    MISSING=0
    for k in "${!PAIR_LABEL[@]}"; do
        H="${PAIR_HAP2[$k]}"
        D="${PAIR_DMP[$k]}"
        PLAB="${PAIR_LABEL[$k]}"
        pair_enabled "$H" "$D" || continue
        H_FASTA="$HAP2_VARIANTS_DIR/${H}.fasta"
        D_FASTA="$DMP_VARIANTS_DIR/${D}.fasta"
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            NHAP="${STOICH_COUNTS[$j]}"
            JOB_DIR="$COMPLEX_DIR/$SLAB/$PLAB"
            mkdir -p "$JOB_DIR"
            EXPECT_CIF="$JOB_DIR/fold_${PLAB}_${SLAB}_model_0.cif"
            printf '%s\t%s\t%s\t%s\t%s\t1\t%s\t%s\t%s\n' \
                "$PLAB" "$H" "$D" "$SLAB" "$NHAP" "$H_FASTA" "$D_FASTA" "$EXPECT_CIF" >> "$MANIFEST"

            # Per-slot drop-zone sentinel. Always rewritten (cheap); edits to
            # the variant table propagate on the next prepare_complexes run.
            cat > "$JOB_DIR/_SUBMIT.txt" <<EOF
AlphaFold3 submission slot
  pair          : $PLAB
  stoichiometry : $SLAB

Submit at: https://alphafoldserver.com/
  HAP2 sequence : $H_FASTA   (paste $NHAP copies)
  DMP sequence  : $D_FASTA   (paste 1 copy)

After the job finishes:
  1. Download the AF3 result ZIP.
  2. Unzip into this directory, preserving the AF3 filenames
     (fold_*_model_0.cif, ranking_scores.json, summary_confidences_*.json).
  3. The orchestrator expects the model-0 file at:
       $EXPECT_CIF
     If AF3 names it differently, rename or symlink to that exact path.

Re-run "bash 14_Interaction_Domain_Mapping.sh" once the drop is in place.
EOF

            # deletion_ladder WT__WT: lift the CIF from stoichiometry_comparison so we
            # don't need a separate AF3 submission for the anchor pair. The monomeric
            # WT__WT complex is identical across both experiments for the same run label.
            if [[ "$EXPERIMENT" == "deletion_ladder" && "$PLAB" == "WT__WT" \
                    && ( ! -s "$EXPECT_CIF" || "$OVERWRITE" == "true" ) ]]; then
                _SIBLING_CIF="$OUT_DIR/stoichiometry_comparison/02_Complexes/$SLAB/$PLAB/fold_${PLAB}_${SLAB}_model_0.cif"
                if [[ -s "$_SIBLING_CIF" ]]; then
                    log_info "  [LIFT] WT__WT / $SLAB: copying CIF from stoichiometry_comparison (no new AF3 job needed)"
                    cp "$_SIBLING_CIF" "$EXPECT_CIF"
                fi
            fi

            # Auto-extract any downloaded AF3 zip in this slot directory.
            # Zip files are always kept as backup - never deleted after extraction.
            # If the extracted CIF has a different name than expected, it is
            # renamed to the canonical EXPECT_CIF path so downstream ops find it.
            _ZIPS=()
            while IFS= read -r _z; do _ZIPS+=("$_z"); done \
                < <(find "$JOB_DIR" -maxdepth 1 -name '*.zip' 2>/dev/null | sort)
            if (( ${#_ZIPS[@]} > 0 )); then
                _ZIP="${_ZIPS[-1]}"   # use the most recent zip if multiple exist
                if [[ ! -s "$EXPECT_CIF" || "$OVERWRITE" == "true" ]]; then
                    log_info "  [ZIP] Extracting $(basename "$_ZIP") -> $JOB_DIR  (zip kept as backup)"
                    unzip -o -q "$_ZIP" -d "$JOB_DIR" \
                        || log_warn "  [ZIP] unzip failed for $(basename "$_ZIP")"
                    # Normalise CIF path: AF3 server may name the file differently.
                    if [[ ! -s "$EXPECT_CIF" ]]; then
                        _FOUND_CIF=$(find "$JOB_DIR" -maxdepth 2 -name '*.cif' \
                            ! -path "$EXPECT_CIF" 2>/dev/null | head -1)
                        if [[ -n "$_FOUND_CIF" ]]; then
                            mv "$_FOUND_CIF" "$EXPECT_CIF"
                            log_info "  [ZIP] Renamed $(basename "$_FOUND_CIF") -> $(basename "$EXPECT_CIF")"
                        fi
                    fi
                else
                    log_info "  [ZIP] $(basename "$_ZIP") present; CIF already extracted (OVERWRITE=false, zip kept)"
                fi
            fi

            if [[ "$PREPARE_BACKEND" == "manual" ]]; then
                if [[ ! -s "$EXPECT_CIF" ]]; then
                    log_warn "  [MISSING] $PLAB / $SLAB ($NHAP HAP2[$H] + 1 DMP[$D]) - upload to https://alphafoldserver.com/ and save to $JOB_DIR/"
                    MISSING=$((MISSING+1))
                else
                    log_info "  [OK]      $PLAB / $SLAB"
                fi
            else
                # Local AF3 wrapper (commercial licence required)
                PRED_WRAP="$MODULES/14_special_pipeline/predict_af3_local.sh"
                conda_run_in "$ENV_PREPARE_COMPLEXES" bash "$PRED_WRAP" \
                    --hap2-fasta "$H_FASTA" \
                    --hap2-copies "$NHAP" \
                    --dmp-fasta "$D_FASTA" \
                    --dmp-copies 1 \
                    --out-dir "$JOB_DIR" \
                    --threads "$CPU"
            fi
        done
    done

    # Top-level drop-zone README explaining the whole tree
    cat > "$DROP_README" <<EOF
AlphaFold3 drop-zone for Stage 14 ($GENE_GROUP)
================================================

Layout:
  $COMPLEX_DIR/
    {stoichiometry}/                  # monomeric, dimeric, trimeric
      {pair_label}/                   # e.g. WT__WT, delC_596_705__WT,
                                      #      WT__delN_1_64
        _SUBMIT.txt                   # this slot's submission instructions
        fold_<pair>_<stoich>_model_0.cif  # drop the AF3 download here

Manifest (machine-readable queue): $MANIFEST

Pair label = "{hap2_variant}__{dmp_variant}". Everything before the double
underscore is the HAP2 truncation; everything after is the DMP truncation.
"WT" on either side means full-length wild-type.

Portal: https://alphafoldserver.com/  (free non-commercial tier; ~20
jobs/day; multimer supported). For each pair-and-stoichiometry slot:
  1. Open a new AF3 job.
  2. Paste the HAP2 sequence (copies per the slot's _SUBMIT.txt) and the
     DMP sequence (1 copy).
  3. Submit, wait, download the ZIP.
  4. Extract into the matching slot directory; the orchestrator picks it
     up on the next run.
EOF

    if [[ "$PREPARE_BACKEND" == "manual" && "$MISSING" -gt 0 ]]; then
        log_warn "$MISSING AF3 prediction(s) still need to be uploaded by hand."
        log_warn "Manifest:      $MANIFEST"
        log_warn "Drop-zone doc: $DROP_README"
        log_warn "Per-slot doc:  $COMPLEX_DIR/{stoich}/{pair}/_SUBMIT.txt"
        log_warn "Re-run this stage after the AF3 downloads are in place."
        if [[ "$DRY_RUN" != "true" ]]; then
            log_warn "Skipping downstream operations for $GENE_GROUP until predictions exist."
            teardown_logging
            continue
        fi
    fi
fi

# ============================================================================
# Operation 3: Interface analysis
# ============================================================================
if op_enabled "interface_analysis"; then
    log_step "Op 3/7: interface_analysis (${#PAIR_LABEL[@]} pairs)"
    ANALYZER="$MODULES/14_special_pipeline/analyze_interface.py"
    IFACE_CUTOFF=$(get_toml interface contact_cutoff_angstrom 2>/dev/null || echo "5.0")
    SASA_PROBE=$(get_toml interface sasa_probe_radius 2>/dev/null || echo "1.4")

    IFACE_PIDS=()
    for k in "${!PAIR_LABEL[@]}"; do
        H="${PAIR_HAP2[$k]}"
        D="${PAIR_DMP[$k]}"
        PLAB="${PAIR_LABEL[$k]}"
        pair_enabled "$H" "$D" || continue
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            CIF="$COMPLEX_DIR/$SLAB/$PLAB/fold_${PLAB}_${SLAB}_model_0.cif"
            OUT_SUB="$IFACE_DIR/$SLAB/$PLAB"
            mkdir -p "$OUT_SUB"
            if [[ ! -s "$CIF" ]]; then
                log_warn "  [skip] missing $CIF"
                continue
            fi
            wait_for_slot "$MAX_PARALLEL"
            (
                conda_run_in "$ENV_IFACE" python3 "$ANALYZER" \
                    --complex "$CIF" \
                    --hap2-chain-prefix A \
                    --dmp-chain-prefix B \
                    --contact-cutoff "$IFACE_CUTOFF" \
                    --sasa-probe "$SASA_PROBE" \
                    --output "$OUT_SUB"
            ) &
            IFACE_PIDS+=($!)
        done
    done
    # Wait only on explicit job PIDs - bare `wait` deadlocks against the
    # tee process-substitutions set up by setup_logging (those tees only
    # exit when the script's stdout closes, which never happens while
    # `wait` is blocking on them).
    (( ${#IFACE_PIDS[@]} > 0 )) && wait "${IFACE_PIDS[@]}"
fi

# ============================================================================
# Operation 4: Short GROMACS MD equilibration (reuses PPI stage-10 modules)
# ============================================================================
if op_enabled "md_equilibration"; then
    log_step "Op 4/7: md_equilibration  (backend=$MD_BACKEND, ${#PAIR_LABEL[@]} pairs)"
    PPI_GROMACS="$PIPELINE_DIR/10_run_gromacs_pipeline.sh"
    MD_WRAP="$MODULES/14_special_pipeline/run_short_md.sh"
    NS=$(get_toml md_equilibration production_ns 2>/dev/null || echo "50")
    USE_MEMBRANE=$(get_toml md_equilibration membrane 2>/dev/null || echo "false")

    for k in "${!PAIR_LABEL[@]}"; do
        H="${PAIR_HAP2[$k]}"
        D="${PAIR_DMP[$k]}"
        PLAB="${PAIR_LABEL[$k]}"
        pair_enabled "$H" "$D" || continue
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            CIF="$COMPLEX_DIR/$SLAB/$PLAB/fold_${PLAB}_${SLAB}_model_0.cif"
            MD_SUB="$MD_DIR/$SLAB/$PLAB"
            mkdir -p "$MD_SUB"
            [[ ! -s "$CIF" ]] && { log_warn "  [skip] $CIF missing"; continue; }
            log_info "  [$PLAB/$SLAB] MD = ${NS} ns, membrane=$USE_MEMBRANE"
            conda_run_in "$ENV_MD" bash "$MD_WRAP" \
                --complex "$CIF" \
                --output "$MD_SUB" \
                --gmx "$GMX_BIN" \
                --threads "$CPU" \
                --production-ns "$NS" \
                --membrane "$USE_MEMBRANE" \
                --ppi-pipeline "$PPI_GROMACS" \
                --config "$CONFIG_FILE"
        done
    done
fi

# ============================================================================
# Operation 5: Binding free energy
# ============================================================================
if op_enabled "binding_energy"; then
    log_step "Op 5/7: binding_energy  (backends=$DG_BACKENDS, ${#PAIR_LABEL[@]} pairs)"
    MMPBSA_WRAP="$MODULES/14_special_pipeline/compute_mmpbsa.sh"
    FOLDX_WRAP="$MODULES/14_special_pipeline/compute_foldx_dg.sh"
    PRODIGY_WRAP="$MODULES/14_special_pipeline/compute_prodigy_dg.sh"

    for k in "${!PAIR_LABEL[@]}"; do
        H="${PAIR_HAP2[$k]}"
        D="${PAIR_DMP[$k]}"
        PLAB="${PAIR_LABEL[$k]}"
        pair_enabled "$H" "$D" || continue
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            CIF="$COMPLEX_DIR/$SLAB/$PLAB/fold_${PLAB}_${SLAB}_model_0.cif"
            MD_SUB="$MD_DIR/$SLAB/$PLAB"
            DG_SUB="$DG_DIR/$SLAB/$PLAB"
            mkdir -p "$DG_SUB"
            [[ ! -s "$CIF" ]] && { log_warn "  [skip] $CIF missing"; continue; }

            if [[ "$DG_BACKENDS" == *mmpbsa* ]]; then
                conda_run_in "$ENV_MMPBSA" bash "$MMPBSA_WRAP" \
                    --md-dir "$MD_SUB" --out "$DG_SUB/mmpbsa" --threads "$CPU"
            fi
            if [[ "$DG_BACKENDS" == *foldx* ]]; then
                if [[ -z "$FOLDX_BIN" || ! -x "$FOLDX_BIN" ]]; then
                    log_warn "  FoldX binary not set or not executable (set [run].show_manual = true in the TOML to print step [M4]); skipping FoldX backend."
                else
                    conda_run_in "$ENV_FOLDX" bash "$FOLDX_WRAP" \
                        --complex "$CIF" --foldx "$FOLDX_BIN" --out "$DG_SUB/foldx"
                fi
            fi
            if [[ "$DG_BACKENDS" == *prodigy* ]]; then
                # --env kept for back-compat with wrappers that activate
                # PRODIGY themselves; the outer conda_run_in is authoritative.
                conda_run_in "$ENV_PRODIGY_DG" bash "$PRODIGY_WRAP" \
                    --complex "$CIF" --env "$ENV_PRODIGY_DG" --out "$DG_SUB/prodigy"
            fi
        done
    done
fi

# ============================================================================
# Operation 6: Computational alanine scan (WT only - this is the hot-spot map)
# ============================================================================
if op_enabled "alanine_scan"; then
    log_step "Op 6/7: alanine_scan  (WT HAP2 / WT DMP pair only)"
    SCAN_WRAP="$MODULES/14_special_pipeline/alanine_scan.sh"
    WT_HAP2=$(get_toml alanine_scan reference_variant 2>/dev/null || echo "WT")
    WT_DMP=$(get_toml alanine_scan reference_dmp_variant 2>/dev/null || echo "WT")
    WT_STOICH=$(get_toml alanine_scan reference_stoichiometry 2>/dev/null || echo "monomeric")
    WT_PLAB="${WT_HAP2}__${WT_DMP}"
    WT_CIF="$COMPLEX_DIR/$WT_STOICH/$WT_PLAB/fold_${WT_PLAB}_${WT_STOICH}_model_0.cif"
    if [[ ! -s "$WT_CIF" ]]; then
        log_warn "  Reference complex $WT_CIF missing; skip alanine scan."
    else
        if [[ -z "$FOLDX_BIN" || ! -x "$FOLDX_BIN" ]]; then
            log_warn "  FoldX binary missing (set [run].show_manual = true in the TOML to print step [M4]); skip alanine scan."
        else
            conda_run_in "$ENV_ASCAN" bash "$SCAN_WRAP" \
                --complex "$WT_CIF" \
                --foldx "$FOLDX_BIN" \
                --interface-tsv "$IFACE_DIR/$WT_STOICH/$WT_PLAB/interface_residues.tsv" \
                --threads "$CPU" \
                --output "$ASCAN_DIR"
        fi
    fi
fi

# ============================================================================
# Operation 7: Comparative report
# ============================================================================
# Reports the headline ranking against [stoichiometry].primary_stoichiometry
# (default "monomeric") - the 1:1 case is the biologically meaningful
# DMP-HAP2 binding measurement; the 2:1 and 3:1 stoichiometries are reported
# alongside but never override the primary case. See NOTES.md Section 7.
if op_enabled "comparative_report"; then
    log_step "  Op 7/7: comparative_report"
    REPORT_WRAP="$MODULES/14_special_pipeline/compare_variants.py"
    PRIMARY_STOICH=$(get_toml stoichiometry basis 2>/dev/null || echo "monomeric")
    conda_run_in "$ENV_REPORT" python3 "$REPORT_WRAP" \
        --hap2-variants "${EXP_HAP2_NAMES[*]}" \
        --dmp-variants "${EXP_DMP_NAMES[*]}" \
        --pair-labels "${PAIR_LABEL[*]}" \
        --pair-hap2 "${PAIR_HAP2[*]}" \
        --pair-dmp "${PAIR_DMP[*]}" \
        --pairing-mode "$PAIRING_MODE" \
        --stoich "${STOICH_LABELS[*]}" \
        --primary-stoich "$PRIMARY_STOICH" \
        --iface-dir "$IFACE_DIR" \
        --dg-dir "$DG_DIR" \
        --ascan-dir "$ASCAN_DIR" \
        --output "$REPORT_DIR" \
        --experiment "$EXPERIMENT"
fi

log_info "  Experiment '$EXPERIMENT' finished. Summary at $REPORT_DIR/"

done  # EXPERIMENT loop (active_comparison)

log_info "Stage 14 finished for $GENE_GROUP. Outputs under $OUT_DIR/"
teardown_logging

done  # GENE_GROUP loop
