#! /bin/bash

:
<<'EOF'
Result to test:

Interactions between HAP2/GCS1 and DMP9.

Arabidopsis HAP2/GCS1 is a type-1 membrane protein of 705 amino acids, with a long extracellular domain, a single TMD, and a rather short and histidine-rich cytoplasmic tail (7, 15, 42). To identify DMP9-interacting domains of HAP2/GCS1, we introduced deletions in its ectodomain and the cytoplasmic carboxy(C)-terminus for testing binding ability. Notably, their interaction with DMP9 was not compromised when HAP2 lacked the C terminus (596-705, leaving 13 amino acids adjacent to TM domain), or when defined regions of the HAP2 ectodomain were deleted according the HAP2 crystal structure (43) (Fig. 3A and SI Appendix, Fig. S7). Even the HAP2 variant without the ectodomain (25-530) was able to weakly interact with DMP9 in yeast (Fig. 3A and SI Appendix, Fig. S8). This suggests large interaction interfaces involving both extra- and intracellular parts of HAP2, and potentially also involving the membrane-embedded domains of HAP2. Thus, our interaction studies support a previous genetic study in Arabidopsis showing that the extra- and intracellular parts of HAP2 are equally important to rescue the fertility in hap2/gcs1 mutant plants (42). Although this study did not examine the subcellular location of truncated HAP2 proteins, the molecular dissection of C. reinhardtii HAP2 showed that its ectodomain is required for trafficking to the cell surface, whereas the cytoplasmic region is crucial to target CrHAP2 to the minus gamete mating structure and to regulate its fusion activity (44). Currently, it remains unclear whether the HAP2 forms expressed in heterologous systems could fold into a native conformation or the unproper folding forms could be generated since large segments of HAP2 were deleted in the MbY2H assays.

Reference Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9659367/

Then, I already have a AlphaFold3 PPI interaction model of;
    monomeric HAP2 with SmelDMPv5_10.610
    dimeric HAP2 with SmelDMPv5_10.610
    trimeric HAP2 with SmelDMPv5_10.610
EOF

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
#   1. Build truncation variants of SmelHAP2 (or AtHAP2) matching the paper.
#   2. Predict each variant + SmelDMPv5_10.610 complex at 1:1, 2:1, 3:1
#      stoichiometry via AlphaFold3.
#   3. Score the predicted interface (BSA, contacts, PRODIGY DG_pred, ipTM).
#   4. Run short GROMACS MD (NVT + NPT + 50 ns production) to relax each
#      complex and confirm the interface is stable, not a docking artefact.
#   5. Decompose binding free energy via gmx_MMPBSA on the MD ensemble and
#      via FoldX AnalyseComplex on the predicted snapshot.
#   6. Computational alanine scan of the WT-complex interface (MutaTeX) to
#      flag hot-spot residues; cross-reference with the regions deleted in
#      the paper.
#   7. Rank variants by predicted DG and interface preservation; emit a
#      heatmap and a summary TSV that mirrors the paper's Fig. 3A panel.
#
# Operations (gated via [domain_mapping].operations in the TOML):
#   prepare_variants    - build truncated HAP2 FASTAs from variant table
#   predict_complexes   - AlphaFold3 multimer prediction (USER ACTION required;
#                         see --show-manual-steps)
#   interface_analysis  - residue contacts, BSA, ipTM, PRODIGY DG_pred
#   md_equilibration    - short GROMACS MD via the existing PPI stage-10 modules
#   binding_energy      - gmx_MMPBSA on the trajectory + FoldX AnalyseComplex
#   alanine_scan        - MutaTeX alanine scan of the WT-complex interface
#   comparative_report  - cross-variant heatmap, ranking, summary TSV / PDF
#
# Output layout under III_RESULT/{GROUP}/14_Domain_Mapping/:
#   01_Variants/                            (truncated HAP2 FASTAs)
#   02_Complexes/{stoich}/{variant}/        (AF3 .cif / ranking_debug.json)
#   03_Interfaces/{stoich}/{variant}/       (interface .tsv, BSA, contacts)
#   04_MD/{stoich}/{variant}/               (GROMACS .gro/.xtc/.edr)
#   05_BindingEnergy/{stoich}/{variant}/    (MM-PBSA, FoldX dG)
#   06_AlanineScan/                         (per-residue DDG; WT only)
#   07_Summary/                             (cross-variant ranking, figures)
#
# References (cite when reusing this pipeline in the manuscript):
#   AlphaFold3       Abramson et al. 2024     10.1038/s41586-024-07487-w
#   GROMACS          Abraham et al. 2015      10.1016/j.softx.2015.06.001
#   gmx_MMPBSA       Valdes-Tresanco 2021     10.1021/acs.jctc.1c00645
#   FoldX            Schymkowitz et al. 2005  10.1093/nar/gki387
#   MutaTeX          Tiberti et al. 2022      10.5281/zenodo.6346203
#   freesasa         Mitternacht 2016         10.12688/f1000research.7931.1
#   PRODIGY          Xue et al. 2016          10.1093/bioinformatics/btw514
#   PyMOL            Schroedinger LLC (no DOI; cite version)
#   Paper basis      Wang et al. 2022         PMC9659367
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

SHARED_CONFIG="$PIPELINE_DIR/14_HAP2_Domain_MappingCONFIG.toml"

# ── CLI flags ───────────────────────────────────────────────────────────────
DRY_RUN=false
SHOW_MANUAL=false
ONLY_OP=""
ONLY_VARIANT=""

print_usage() {
    cat <<USAGE
Usage: bash 14_Special_Pipeline.sh [options]

Options:
  --dry-run                Print the actions for each stage without running them.
  --show-manual-steps      Print the steps you (the operator) must do by hand,
                           then exit. Use this to scope wet-lab-equivalent work
                           (AlphaFold3 submissions, FoldX licensing, etc.).
  --op <name>              Run only one operation (e.g. interface_analysis).
                           Repeatable.
  --variant <name>         Restrict variant iteration to one name from
                           [variants].names. Repeatable.
  -h, --help               Show this message.
USAGE
}

ONLY_OPS=()
ONLY_VARIANTS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)            DRY_RUN=true; shift ;;
        --show-manual-steps)  SHOW_MANUAL=true; shift ;;
        --op)                 ONLY_OPS+=("$2"); shift 2 ;;
        --variant)            ONLY_VARIANTS+=("$2"); shift 2 ;;
        -h|--help)            print_usage; exit 0 ;;
        *)  echo "ERROR: unknown flag '$1'" >&2; print_usage; exit 1 ;;
    esac
done

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
       run by pointing [predict_complexes].backend = "local" and providing
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
    echo "ERROR: pipeline.gene_groups is empty in 14_special_pipelineCONFIG.toml" >&2
    exit 1
fi

wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done
}

op_enabled() {
    local op="$1"
    # --op filter overrides TOML when present
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

run_or_echo() {
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] $*"
    else
        log_info "Running: $*"
        "$@"
    fi
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

# Config resolution: shared root + optional per-group override (cat-merged)
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
PER_GROUP_CFG="$CONFIG_DIR/14_special_pipeline_${GENE_GROUP}.toml"
if [[ -f "$PER_GROUP_CFG" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_special_cfg_XXXXXX.toml")
    cat "$SHARED_CONFIG" "$PER_GROUP_CFG" > "$CONFIG_FILE"
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
VARIANTS_DIR="$OUT_DIR/01_Variants"
COMPLEX_DIR="$OUT_DIR/02_Complexes"
IFACE_DIR="$OUT_DIR/03_Interfaces"
MD_DIR="$OUT_DIR/04_MD"
DG_DIR="$OUT_DIR/05_BindingEnergy"
ASCAN_DIR="$OUT_DIR/06_AlanineScan"
REPORT_DIR="$OUT_DIR/07_Summary"
mkdir -p "$VARIANTS_DIR" "$COMPLEX_DIR" "$IFACE_DIR" "$MD_DIR" "$DG_DIR" "$ASCAN_DIR" "$REPORT_DIR"

setup_logging

ENABLED=$(get_toml domain_mapping enabled 2>/dev/null || echo "true")
if [[ "$ENABLED" != "true" && "$ENABLED" != "True" ]]; then
    log_info "Domain mapping disabled for $GENE_GROUP. Skipping."
    teardown_logging
    continue
fi

mapfile -t OPERATIONS < <(get_toml domain_mapping operations 2>/dev/null \
    || printf '%s\n' prepare_variants predict_complexes interface_analysis \
                     md_equilibration binding_energy alanine_scan comparative_report)

# Inputs
HAP2_FASTA=$(get_toml inputs hap2_fasta 2>/dev/null)
DMP_FASTA=$(get_toml inputs dmp_fasta 2>/dev/null)
HAP2_REF_PDB=$(get_toml inputs hap2_reference_pdb 2>/dev/null || echo "")
[[ "$HAP2_FASTA" != /* && -n "$HAP2_FASTA" ]] && HAP2_FASTA="$PIPELINE_DIR/$HAP2_FASTA"
[[ "$DMP_FASTA"  != /* && -n "$DMP_FASTA"  ]] && DMP_FASTA="$PIPELINE_DIR/$DMP_FASTA"
if [[ ! -f "$HAP2_FASTA" || ! -f "$DMP_FASTA" ]]; then
    log_error "Required input FASTA missing. HAP2='$HAP2_FASTA' DMP='$DMP_FASTA'"
    log_error "See --show-manual-steps [M1] for guidance on selecting source sequences."
    teardown_logging
    continue
fi

# Variant table (parallel arrays so parse_toml.py can read them)
mapfile -t VARIANT_NAMES        < <(get_toml variants names 2>/dev/null)
mapfile -t VARIANT_DESCRIPTIONS < <(get_toml variants descriptions 2>/dev/null)
mapfile -t VARIANT_DELETIONS    < <(get_toml variants deletions 2>/dev/null)
if [[ ${#VARIANT_NAMES[@]} -eq 0 ]]; then
    log_error "[variants].names is empty - nothing to map."
    teardown_logging
    continue
fi

# Stoichiometry (HAP2 copy counts to predict)
mapfile -t STOICH_LABELS < <(get_toml stoichiometry labels 2>/dev/null \
    || printf '%s\n' monomeric dimeric trimeric)
mapfile -t STOICH_COUNTS < <(get_toml stoichiometry chain_counts 2>/dev/null \
    || printf '%s\n' 1 2 3)

# Tooling
PREDICT_BACKEND=$(get_toml predict_complexes backend 2>/dev/null || echo "manual")
MD_BACKEND=$(get_toml md_equilibration backend 2>/dev/null || echo "gromacs")
DG_BACKENDS=$(get_toml binding_energy backends 2>/dev/null || echo "mmpbsa foldx prodigy")
FOLDX_BIN=$(get_toml tools foldx_binary 2>/dev/null || echo "")
GMX_BIN=$(get_toml tools gromacs_binary 2>/dev/null || echo "gmx")
PRODIGY_ENV=$(get_toml tools prodigy_conda_env 2>/dev/null || echo "egg")

log_step "Stage 14 (Domain Mapping): $GENE_GROUP"
log_info "Operations:    ${OPERATIONS[*]}"
log_info "Variants:      ${VARIANT_NAMES[*]}"
log_info "Stoichiometry: ${STOICH_LABELS[*]} (${STOICH_COUNTS[*]} HAP2 copies)"
log_info "Predict backend: $PREDICT_BACKEND   MD backend: $MD_BACKEND"
log_info "Compute: host=${HOST_CPU}c  CPU=$CPU  MAX_PARALLEL=$MAX_PARALLEL"
log_info "Output:        $OUT_DIR"

# ============================================================================
# Operation 1: Build truncated HAP2 FASTAs
# ============================================================================
if op_enabled "prepare_variants"; then
    log_step "Op 1/7: prepare_variants"
    GEN_SCRIPT="$MODULES/14_special_pipeline/generate_variants.py"
    for i in "${!VARIANT_NAMES[@]}"; do
        VNAME="${VARIANT_NAMES[$i]}"
        variant_enabled "$VNAME" || continue
        VDESC="${VARIANT_DESCRIPTIONS[$i]:-no description}"
        VDEL="${VARIANT_DELETIONS[$i]:-}"   # empty = WT
        OUT_FA="$VARIANTS_DIR/${VNAME}.fasta"
        if [[ -f "$OUT_FA" && "$OVERWRITE" != "true" ]]; then
            log_info "  [$VNAME] exists, skip (OVERWRITE=false)"
            continue
        fi
        log_info "  [$VNAME] $VDESC  (delete='$VDEL')"
        run_or_echo python3 "$GEN_SCRIPT" \
            --input "$HAP2_FASTA" \
            --name "$VNAME" \
            --description "$VDESC" \
            --deletions "$VDEL" \
            --output "$OUT_FA"
    done
fi

# ============================================================================
# Operation 2: AlphaFold3 multimer prediction
# ============================================================================
# `backend = "manual"` is the realistic default - AF3 has no free programmatic
# inference API as of 2026-05. The orchestrator stages the work, prints the
# submission checklist, and confirms each download exists before moving on.
if op_enabled "predict_complexes"; then
    log_step "Op 2/7: predict_complexes  (backend=$PREDICT_BACKEND)"
    MANIFEST="$COMPLEX_DIR/_submission_manifest.tsv"
    : > "$MANIFEST"
    printf 'variant\tstoich\tn_hap2\tn_dmp\tvariant_fasta\tdmp_fasta\texpected_output\n' >> "$MANIFEST"

    MISSING=0
    for i in "${!VARIANT_NAMES[@]}"; do
        VNAME="${VARIANT_NAMES[$i]}"
        variant_enabled "$VNAME" || continue
        V_FASTA="$VARIANTS_DIR/${VNAME}.fasta"
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            NHAP="${STOICH_COUNTS[$j]}"
            JOB_DIR="$COMPLEX_DIR/$SLAB/$VNAME"
            mkdir -p "$JOB_DIR"
            EXPECT_CIF="$JOB_DIR/fold_${VNAME}_${SLAB}_model_0.cif"
            printf '%s\t%s\t%s\t1\t%s\t%s\t%s\n' \
                "$VNAME" "$SLAB" "$NHAP" "$V_FASTA" "$DMP_FASTA" "$EXPECT_CIF" >> "$MANIFEST"

            if [[ "$PREDICT_BACKEND" == "manual" ]]; then
                if [[ ! -s "$EXPECT_CIF" ]]; then
                    log_warn "  [MISSING] $VNAME / $SLAB ($NHAP HAP2 + 1 DMP) - upload to https://alphafoldserver.com/ and save to $JOB_DIR/"
                    MISSING=$((MISSING+1))
                else
                    log_info "  [OK]      $VNAME / $SLAB"
                fi
            else
                # Local AF3 wrapper (commercial licence required)
                PRED_WRAP="$MODULES/14_special_pipeline/predict_af3_local.sh"
                run_or_echo bash "$PRED_WRAP" \
                    --hap2-fasta "$V_FASTA" \
                    --hap2-copies "$NHAP" \
                    --dmp-fasta "$DMP_FASTA" \
                    --dmp-copies 1 \
                    --out-dir "$JOB_DIR" \
                    --threads "$CPU"
            fi
        done
    done

    if [[ "$PREDICT_BACKEND" == "manual" && "$MISSING" -gt 0 ]]; then
        log_warn "$MISSING AF3 prediction(s) still need to be uploaded by hand."
        log_warn "Manifest written to $MANIFEST"
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
    log_step "Op 3/7: interface_analysis"
    ANALYZER="$MODULES/14_special_pipeline/analyze_interface.py"
    IFACE_CUTOFF=$(get_toml interface contact_cutoff_angstrom 2>/dev/null || echo "5.0")
    SASA_PROBE=$(get_toml interface sasa_probe_radius 2>/dev/null || echo "1.4")

    for i in "${!VARIANT_NAMES[@]}"; do
        VNAME="${VARIANT_NAMES[$i]}"
        variant_enabled "$VNAME" || continue
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            CIF="$COMPLEX_DIR/$SLAB/$VNAME/fold_${VNAME}_${SLAB}_model_0.cif"
            OUT_SUB="$IFACE_DIR/$SLAB/$VNAME"
            mkdir -p "$OUT_SUB"
            if [[ ! -s "$CIF" ]]; then
                log_warn "  [skip] missing $CIF"
                continue
            fi
            wait_for_slot "$MAX_PARALLEL"
            (
                run_or_echo conda run -n "$PRODIGY_ENV" python3 "$ANALYZER" \
                    --complex "$CIF" \
                    --hap2-chain-prefix A \
                    --dmp-chain-prefix B \
                    --contact-cutoff "$IFACE_CUTOFF" \
                    --sasa-probe "$SASA_PROBE" \
                    --output "$OUT_SUB"
            ) &
        done
    done
    wait
fi

# ============================================================================
# Operation 4: Short GROMACS MD equilibration (reuses PPI stage-10 modules)
# ============================================================================
if op_enabled "md_equilibration"; then
    log_step "Op 4/7: md_equilibration  (backend=$MD_BACKEND)"
    PPI_GROMACS="$PIPELINE_DIR/10_run_gromacs_pipeline.sh"
    MD_WRAP="$MODULES/14_special_pipeline/run_short_md.sh"
    NS=$(get_toml md_equilibration production_ns 2>/dev/null || echo "50")
    USE_MEMBRANE=$(get_toml md_equilibration membrane 2>/dev/null || echo "false")

    for i in "${!VARIANT_NAMES[@]}"; do
        VNAME="${VARIANT_NAMES[$i]}"
        variant_enabled "$VNAME" || continue
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            CIF="$COMPLEX_DIR/$SLAB/$VNAME/fold_${VNAME}_${SLAB}_model_0.cif"
            MD_SUB="$MD_DIR/$SLAB/$VNAME"
            mkdir -p "$MD_SUB"
            [[ ! -s "$CIF" ]] && { log_warn "  [skip] $CIF missing"; continue; }
            log_info "  [$VNAME/$SLAB] MD = ${NS} ns, membrane=$USE_MEMBRANE"
            run_or_echo bash "$MD_WRAP" \
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
    log_step "Op 5/7: binding_energy  (backends=$DG_BACKENDS)"
    MMPBSA_WRAP="$MODULES/14_special_pipeline/compute_mmpbsa.sh"
    FOLDX_WRAP="$MODULES/14_special_pipeline/compute_foldx_dg.sh"
    PRODIGY_WRAP="$MODULES/14_special_pipeline/compute_prodigy_dg.sh"

    for i in "${!VARIANT_NAMES[@]}"; do
        VNAME="${VARIANT_NAMES[$i]}"
        variant_enabled "$VNAME" || continue
        for j in "${!STOICH_LABELS[@]}"; do
            SLAB="${STOICH_LABELS[$j]}"
            CIF="$COMPLEX_DIR/$SLAB/$VNAME/fold_${VNAME}_${SLAB}_model_0.cif"
            MD_SUB="$MD_DIR/$SLAB/$VNAME"
            DG_SUB="$DG_DIR/$SLAB/$VNAME"
            mkdir -p "$DG_SUB"
            [[ ! -s "$CIF" ]] && { log_warn "  [skip] $CIF missing"; continue; }

            if [[ "$DG_BACKENDS" == *mmpbsa* ]]; then
                run_or_echo bash "$MMPBSA_WRAP" \
                    --md-dir "$MD_SUB" --out "$DG_SUB/mmpbsa" --threads "$CPU"
            fi
            if [[ "$DG_BACKENDS" == *foldx* ]]; then
                if [[ -z "$FOLDX_BIN" || ! -x "$FOLDX_BIN" ]]; then
                    log_warn "  FoldX binary not set or not executable (see --show-manual-steps [M4]); skipping FoldX backend."
                else
                    run_or_echo bash "$FOLDX_WRAP" \
                        --complex "$CIF" --foldx "$FOLDX_BIN" --out "$DG_SUB/foldx"
                fi
            fi
            if [[ "$DG_BACKENDS" == *prodigy* ]]; then
                run_or_echo bash "$PRODIGY_WRAP" \
                    --complex "$CIF" --env "$PRODIGY_ENV" --out "$DG_SUB/prodigy"
            fi
        done
    done
fi

# ============================================================================
# Operation 6: Computational alanine scan (WT only - this is the hot-spot map)
# ============================================================================
if op_enabled "alanine_scan"; then
    log_step "Op 6/7: alanine_scan  (WT only)"
    SCAN_WRAP="$MODULES/14_special_pipeline/alanine_scan.sh"
    WT_NAME=$(get_toml alanine_scan reference_variant 2>/dev/null || echo "WT")
    WT_STOICH=$(get_toml alanine_scan reference_stoichiometry 2>/dev/null || echo "monomeric")
    WT_CIF="$COMPLEX_DIR/$WT_STOICH/$WT_NAME/fold_${WT_NAME}_${WT_STOICH}_model_0.cif"
    if [[ ! -s "$WT_CIF" ]]; then
        log_warn "  Reference complex $WT_CIF missing; skip alanine scan."
    else
        if [[ -z "$FOLDX_BIN" || ! -x "$FOLDX_BIN" ]]; then
            log_warn "  FoldX binary missing (see --show-manual-steps [M4]); skip alanine scan."
        else
            run_or_echo bash "$SCAN_WRAP" \
                --complex "$WT_CIF" \
                --foldx "$FOLDX_BIN" \
                --interface-tsv "$IFACE_DIR/$WT_STOICH/$WT_NAME/interface_residues.tsv" \
                --threads "$CPU" \
                --output "$ASCAN_DIR"
        fi
    fi
fi

# ============================================================================
# Operation 7: Comparative report
# ============================================================================
if op_enabled "comparative_report"; then
    log_step "Op 7/7: comparative_report"
    REPORT_WRAP="$MODULES/14_special_pipeline/compare_variants.py"
    run_or_echo conda run -n "$PRODIGY_ENV" python3 "$REPORT_WRAP" \
        --variants "${VARIANT_NAMES[*]}" \
        --stoich "${STOICH_LABELS[*]}" \
        --iface-dir "$IFACE_DIR" \
        --dg-dir "$DG_DIR" \
        --ascan-dir "$ASCAN_DIR" \
        --output "$REPORT_DIR"
fi

log_info "Stage 14 finished for $GENE_GROUP. Summary at $REPORT_DIR/"
teardown_logging

done  # GENE_GROUP loop
