#!/usr/bin/env bash
# Build the HAP2 ectodomain-deletion grouped figures for the Stage-14
# deletion_ladder experiment. Each grouped figure stitches together
# (deletion-red highlight on the left, AF3 model overview on the right)
# for every HAP2 deletion variant, one row per variant, two columns.
#
# Produces three figures per conformation (monomeric, postfusion):
#   *_grouped_<conformation>.png          all 7 variants
#   *_grouped_<conformation>_single.png   single-domain deletions (4)
#   *_grouped_<conformation>_combined.png combined-domain deletions (3)
#
# Variant order within each figure is N-to-C by deletion start position.
# Single = one contiguous region removed. Combined = multiple regions /
# multiple domains removed in the same construct (names contain "And").
#
# Usage:
#   bash modules/14_special_pipeline/make_grouped_deletion_figures.sh
#   bash modules/14_special_pipeline/make_grouped_deletion_figures.sh DMP_x_AtHAP2

set -euo pipefail

GENE_GROUP="${1:-DMP_x_SmelHAP2}"

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMBINE_PY="$PIPELINE_DIR/modules/14_special_pipeline/combine_panels.py"

HAP2_VARIANTS_DIR="$PIPELINE_DIR/III_RESULT/$GENE_GROUP/14_Domain_Mapping/deletion_ladder/07_Summary/deletion_variants/HAP2_variants"
RED_MONO="$HAP2_VARIANTS_DIR/deletion_red_monomeric"
RED_POST="$HAP2_VARIANTS_DIR/deletion_red_postfusion"
MODEL_MONO="$HAP2_VARIANTS_DIR/model_monomeric"
MODEL_POST="$HAP2_VARIANTS_DIR/model_postfusion"

for d in "$RED_MONO" "$RED_POST" "$MODEL_MONO" "$MODEL_POST"; do
    [[ -d "$d" ]] || { echo "[ERROR] missing input directory: $d" >&2; exit 1; }
done
[[ -f "$COMBINE_PY" ]] || { echo "[ERROR] missing $COMBINE_PY" >&2; exit 1; }

# N-to-C order (AtHAP2 start residues):
#   delEcto         25-530
#   delEctoAndC     25-530,596-705
#   delEctoD2       189-330
#   delFL           234-250
#   delPreTMD       531-560
#   delPreTMDAndTMD 531-582
#   delPreTMDAndTMDAndJuxtaMem 531-595
SINGLE=(delEcto delEctoD2 delFL delPreTMD)
COMBINED=(delEctoAndC delPreTMDAndTMD delPreTMDAndTMDAndJuxtaMem)
ALL=(delEcto delEctoAndC delEctoD2 delFL delPreTMD delPreTMDAndTMD delPreTMDAndTMDAndJuxtaMem)

# Resolve the (red, model) panel pair for one variant in one conformation.
# Sets _RED_PNG and _MODEL_PNG; aborts if either is missing.
resolve_pair() {
    local variant="$1" conf="$2"
    local red_dir model_dir red_png model_png
    case "$conf" in
        monomeric)
            red_dir="$RED_MONO"
            model_dir="$MODEL_MONO"
            red_png="$red_dir/hap2_dmp_${variant}_red_monomeric_labeled.png"
            model_png="$model_dir/hap2_dmp_${variant}_model_monomeric_structure_cropped.png"
            ;;
        postfusion)
            red_dir="$RED_POST"
            model_dir="$MODEL_POST"
            red_png="$red_dir/hap2_dmp_${variant}_red_labeled.png"
            model_png="$model_dir/hap2_dmp_${variant}_model_postfusion_structure_cropped.png"
            ;;
        *)
            echo "[ERROR] unknown conformation: $conf" >&2
            return 1
            ;;
    esac
    [[ -f "$red_png"   ]] || { echo "[ERROR] missing red panel:   $red_png"   >&2; return 1; }
    [[ -f "$model_png" ]] || { echo "[ERROR] missing model panel: $model_png" >&2; return 1; }
    _RED_PNG="$red_png"
    _MODEL_PNG="$model_png"
}

# Build one grouped figure. Args: conformation, suffix (""/_single/_combined), variant names...
build_figure() {
    local conf="$1" suffix="$2"; shift 2
    local variants=("$@")
    local out="$HAP2_VARIANTS_DIR/hap2_dmp_ectodomain_deletions_grouped_${conf}${suffix}.png"
    local -a args=()
    local v
    for v in "${variants[@]}"; do
        resolve_pair "$v" "$conf"
        args+=(--image "$_RED_PNG" --caption "$v (red)")
        args+=(--image "$_MODEL_PNG" --caption "$v (model)")
    done
    echo "[INFO] building $(basename "$out")  (${#variants[@]} variants, ${#args[@]} args)"
    python3 "$COMBINE_PY" "${args[@]}" --out "$out" --cols 2
}

for CONF in monomeric postfusion; do
    build_figure "$CONF" ""          "${ALL[@]}"
    build_figure "$CONF" "_single"   "${SINGLE[@]}"
    build_figure "$CONF" "_combined" "${COMBINED[@]}"
done

echo "[OK] grouped deletion figures written under $HAP2_VARIANTS_DIR"
