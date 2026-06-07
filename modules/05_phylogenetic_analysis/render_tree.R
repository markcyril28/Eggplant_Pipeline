#!/usr/bin/env Rscript
# ============================================================================
# Phylogenetic Tree Visualization using ggtree
# ============================================================================
# Suppress fontconfig warnings on headless Linux / WSL2
options(bitmapType = "cairo")
Sys.setenv(FONTCONFIG_PATH = "/etc/fonts")
# Renders publication-quality phylogenetic tree figures from Newick tree files.
#
# Usage:
#   Rscript render_tree.R \
#       --input <tree.treefile> \
#       --output <output.png> \
#       [--layout rectangular] \
#       [--show-bootstrap true] \
#       [--bootstrap-threshold 70] \
#       [--bootstrap-style text]        # text (default) | dots (colored circles)
#       [--label-style replace]         # replace (default) | append
#       [--node-color-high "#B71C1C"] \
#       [--node-color-medium "#E65100"] \
#       [--node-size-high 3.0] \
#       [--node-size-medium 2.0] \
#       [--tip-label-size 4.5] \
#       [--bootstrap-label-size 3.0] \
#       [--width 14] \
#       [--height 10] \
#       [--dpi 600] \
#       [--branch-width 0.8] \
#       [--title ""] \
#       [--highlight-eggplant true] \
#       [--outgroup-pattern "Outgroup"] \
#       [--root-outgroup false]
# ============================================================================

suppressPackageStartupMessages({
    library(ggtree)
    library(treeio)
    library(ggplot2)
    library(argparse)
})

# ======================== Argument Parsing ========================

parser <- ArgumentParser(description = "Render phylogenetic tree figures with ggtree")
parser$add_argument("--input",               required = TRUE,
    help = "Input Newick tree file (.treefile, .contree, .raxml.support, .nwk)")
parser$add_argument("--output",              required = TRUE,
    help = "Output image file path (.png)")
parser$add_argument("--layout",              default = "rectangular",
    help = "Tree layout: rectangular, circular, fan, daylight, equal_angle (default: rectangular)")
parser$add_argument("--show-bootstrap",      default = "true",
    help = "Show bootstrap support on nodes (default: true)")
parser$add_argument("--bootstrap-threshold", type = "integer", default = 70L,
    help = "Only show bootstrap values >= threshold (default: 70)")
parser$add_argument("--bootstrap-style",     default = "dots",
    help = "Bootstrap display: dots (colored circles, default) or text")
parser$add_argument("--label-style",         default = "replace",
    help = "Tip label style: replace (clean short names, default) or append (legacy brackets)")
parser$add_argument("--tip-label-size",      type = "double", default = 4.5,
    help = "Tip label font size (default: 4.5)")
parser$add_argument("--bootstrap-label-size", type = "double", default = 3.0,
    help = "Bootstrap text label font size — text mode only (default: 3.0)")
parser$add_argument("--node-color-high",     default = "#B71C1C",
    help = "Node dot color for BS >= 95 (default: #B71C1C deep red)")
parser$add_argument("--node-color-medium",   default = "#E65100",
    help = "Node dot color for BS >= threshold but < 95 (default: #E65100 deep orange)")
parser$add_argument("--node-size-high",      type = "double", default = 3.0,
    help = "Node dot size for high-confidence nodes (default: 3.0)")
parser$add_argument("--node-size-medium",    type = "double", default = 2.0,
    help = "Node dot size for medium-confidence nodes (default: 2.0)")
parser$add_argument("--width",               type = "double", default = 14.0,
    help = "Figure width in inches (default: 14)")
parser$add_argument("--height",              type = "double", default = 10.0,
    help = "Figure height in inches (default: 10)")
parser$add_argument("--dpi",                 type = "integer", default = 600L,
    help = "Figure resolution in DPI (default: 600)")
parser$add_argument("--branch-width",        type = "double", default = 0.8,
    help = "Branch line width (default: 0.8)")
parser$add_argument("--title",               default = "",
    help = "Optional figure title")
parser$add_argument("--highlight-eggplant",  default = "true",
    help = "Color-code eggplant (SMEL5) and key ortholog tips (default: true)")
parser$add_argument("--outgroup-pattern",    default = "Outgroup",
    help = "Regex pattern to identify outgroup taxa (default: Outgroup)")
parser$add_argument("--root-outgroup",       default = "false",
    help = "Root tree on outgroup taxa (default: false)")
parser$add_argument("--exclude-tips",        default = "",
    help = "Comma-separated regex patterns of tips to drop before rendering (default: none)")
parser$add_argument("--open-angle",          type = "double", default = 15.0,
    help = "Fan layout open angle in degrees (default: 15)")
parser$add_argument("--bootstrap-color",     default = "#D32F2F",
    help = "Color for bootstrap text labels — text mode only (default: #D32F2F)")
parser$add_argument("--branch-color",        default = "grey20",
    help = "Branch line color (default: grey20)")
parser$add_argument("--tip-label-offset",    type = "double", default = 0.005,
    help = "Horizontal offset of tip labels; negative = left, positive = right (default: 0.005)")
parser$add_argument("--tip-point-size",      type = "double", default = 2.5,
    help = "Tip point symbol size (default: 2.5)")
parser$add_argument("--color-smeldmp",       default = "#4527A0",
    help = "Tip color for SmelDMP eggplant genes (default: #4527A0 deep purple)")
parser$add_argument("--color-haploid",       default = "#B71C1C",
    help = "Tip color for validated haploid inducers (default: #B71C1C deep red)")
parser$add_argument("--color-ortholog",      default = "#37474F",
    help = "Tip color for other DMP orthologs (default: #37474F blue-gray)")
parser$add_argument("--color-outgroup",      default = "#78909C",
    help = "Tip color for outgroup taxa (default: #78909C blue-gray)")
parser$add_argument("--treescale-fontsize",  type = "double", default = 3.0,
    help = "Scale bar font size (default: 3.0)")
parser$add_argument("--treescale-offset",    type = "double", default = 0.3,
    help = "Scale bar vertical offset (default: 0.3)")
parser$add_argument("--treescale-color",     default = "grey35",
    help = "Scale bar color (default: grey35)")
parser$add_argument("--phylo-software",      default = "",
    help = "Phylogenetic software name (e.g., IQ-TREE2, RAxML-NG, MEGA_CC)")
parser$add_argument("--phylo-model",         default = "",
    help = "Substitution model used (e.g., LG+F+R4, GTR+F+R4)")
parser$add_argument("--phylo-bootstrap",     default = "",
    help = "Bootstrap method and replicates (e.g., UFBoot 5000 / SH-aLRT 5000)")
parser$add_argument("--sequence-type",       default = "",
    help = "Sequence type: AA or NT")

args <- parser$parse_args()

# ======================== Validate Inputs ========================

if (!file.exists(args$input)) {
    stop(sprintf("Input tree file not found: %s", args$input))
}

valid_layouts <- c("rectangular", "circular", "fan", "daylight", "equal_angle")
if (!(args$layout %in% valid_layouts)) {
    stop(sprintf("Invalid layout '%s'. Choose from: %s",
                 args$layout, paste(valid_layouts, collapse = ", ")))
}

output_dir <- dirname(args$output)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ======================== Read Tree ========================

cat(sprintf("Reading tree: %s\n", args$input))
tree <- read.tree(args$input)

if (is.null(tree)) {
    stop("Failed to parse tree file")
}

n_tips <- length(tree$tip.label)
cat(sprintf("  Tips: %d\n", n_tips))

# ======================== Exclude Tips (optional) ========================

tree_was_modified <- FALSE

if (nchar(args$exclude_tips) > 0) {
    patterns <- trimws(strsplit(args$exclude_tips, ",")[[1]])
    tips_to_drop <- c()
    for (pat in patterns) {
        if (nchar(pat) > 0) {
            matched <- grep(pat, tree$tip.label, value = TRUE)
            tips_to_drop <- c(tips_to_drop, matched)
        }
    }
    tips_to_drop <- unique(tips_to_drop)
    if (length(tips_to_drop) > 0) {
        cat(sprintf("  Excluding %d tip(s): %s\n",
                    length(tips_to_drop), paste(tips_to_drop, collapse = ", ")))
        tree <- ape::drop.tip(tree, tips_to_drop)
        n_tips <- length(tree$tip.label)
        cat(sprintf("  Tips after exclusion: %d\n", n_tips))
        tree_was_modified <- TRUE
    }
}

# ======================== Auto-size ========================
# Scale figure height with number of tips for readability

if (n_tips > 30 && args$height == 10.0) {
    args$height <- max(10, n_tips * 0.4)
    cat(sprintf("  Auto-scaled height to %.1f inches for %d tips\n", args$height, n_tips))
}
if (n_tips > 30 && args$tip_label_size == 4.5) {
    args$tip_label_size <- max(3.5, 5.5 - (n_tips / 100))
    cat(sprintf("  Auto-scaled tip label size to %.1f\n", args$tip_label_size))
}

# ======================== Root on Outgroup (optional) ========================

if (tolower(args$root_outgroup) == "true" && nchar(args$outgroup_pattern) > 0) {
    outgroup_tips <- grep(args$outgroup_pattern, tree$tip.label, value = TRUE)
    if (length(outgroup_tips) > 0) {
        cat(sprintf("  Rooting on %d outgroup taxa matching '%s'\n",
                    length(outgroup_tips), args$outgroup_pattern))
        if (length(outgroup_tips) == 1L) {
            # Single outgroup tip: root directly on it
            tree <- ape::root(tree, outgroup = outgroup_tips, resolve.root = TRUE)
        } else {
            # Multiple outgroup taxa: find their MRCA and root on that branch.
            # Guard: if the outgroups span the current root (MRCA == root node),
            # root(node=) throws "ambiguous resolution" — fall back to the first
            # outgroup tip so at least one anchor is placed correctly.
            mrca_node <- ape::getMRCA(tree, outgroup_tips)
            root_node <- ape::Ntip(tree) + 1L
            if (!is.null(mrca_node) && mrca_node != root_node) {
                tree <- ape::root(tree, node = mrca_node, resolve.root = TRUE)
            } else {
                cat(sprintf("  Note: outgroup MRCA is the root node (outgroups span root); rooting on first outgroup tip '%s'\n",
                            outgroup_tips[1]))
                tree <- ape::root(tree, outgroup = outgroup_tips[1], resolve.root = TRUE)
            }
        }
        tree_was_modified <- TRUE
        # The new root node has no evaluable bootstrap support; clear its label so
        # it does not appear as a spurious value in the bootstrap display layer.
        if (!is.null(tree$node.label) && length(tree$node.label) > 0L) {
            tree$node.label[1L] <- ""
        }
    } else {
        stop(sprintf(paste0(
            "Rooting requested (root_outgroup = true) but no tips match outgroup ",
            "pattern '%s'. Refusing to render an UNROOTED tree as if it were rooted. ",
            "Add the outgroup to the alignment and re-infer the tree, or set ",
            "visualization.root_outgroup = false to render unrooted on purpose."),
            args$outgroup_pattern))
    }
}

# ======================== Save Modified Tree ========================
# Write pruned/rooted tree as a versioned copy — original file untouched.

if (tree_was_modified) {
    # Strip any existing _rooted chains to prevent cascading suffixes
    clean_input <- sub("(_rooted)+\\.nwk$", ".nwk", args$input)
    mod_path <- sub("\\.[^.]+$", "_rooted.nwk", clean_input)
    # Handle double extensions like .raxml.support
    if (grepl("\\.raxml\\.", clean_input)) {
        mod_path <- sub("\\.raxml\\.[^.]+$", "_rooted.nwk", clean_input)
    }
    ape::write.tree(tree, file = mod_path)
    cat(sprintf("  Saved modified tree: %s\n", mod_path))
}

# ======================== Parse Bootstrap Values ========================
# IQ-TREE2 .treefile: labels like "97.8/100" (SH-aLRT/UFBoot) — take last value
# IQ-TREE2 .contree / RAxML .raxml.support: single integer

has_bootstrap <- !is.null(tree$node.label) && length(tree$node.label) > 0

extract_last_boot <- function(x) {
    if (is.na(x) || nchar(trimws(as.character(x))) == 0) return(NA_real_)
    parts <- strsplit(as.character(x), "/")[[1]]
    suppressWarnings(as.numeric(parts[length(parts)]))
}

if (has_bootstrap) {
    boot_numeric <- sapply(tree$node.label, extract_last_boot, USE.NAMES = FALSE)
    valid_boot   <- boot_numeric[!is.na(boot_numeric)]
    if (length(valid_boot) > 0) {
        cat(sprintf("  Bootstrap range: %.0f - %.0f\n", min(valid_boot), max(valid_boot)))
    }
}

# ======================== Gene Category Classification ========================

# SmelDMP short name map, manuscript v5 nomenclature.
smeldmp_name_map <- c(
    "SMEL5_01g008730.1" = "SmelDMPv5_01.730",
    "SMEL5_01g026030.1" = "SmelDMPv5_01.030",
    "SMEL5_02g013320.1" = "SmelDMPv5_02.320",
    "SMEL5_04g005390.1" = "SmelDMPv5_04.390",
    "SMEL5_10g003660.1" = "SmelDMPv5_10.660",
    "SMEL5_10g017610.1" = "SmelDMPv5_10.610",
    "SMEL5_12g005350.1" = "SmelDMPv5_12.350"
)

# Functionally validated haploid-inducing DMP genes — drawn from
# II_INPUTS/DMP_HI_registry.tsv. Patterns are matched as substrings (fixed=TRUE)
# against ORIGINAL tip labels, so both raw locus IDs and short gene symbols
# are listed to cover pre- and post-shortened forms.
# NOTE: "CsDMP9"/XP_006482605 is Citrus sinensis — NOT Cucumis sativus.
# Yin et al. 2024 cucumber HI is CsaV3_1G028660 (which IS listed below).
haploid_inducer_patterns <- c(
    "AtDMP8",   "AT1G09157",                                    # AtDMP8 Zhong et al. 2020
    "AtDMP9",   "AT5G39650",                                    # AtDMP9 Zhong et al. 2020
    "GmDMP",   "Glyma.18G097400", "Glyma.18G098300",            # Zhong et al. 2024
    "NtDMP",                                                    # X. Zhang et al. 2022
    "SlDMP3",  "SlDMP8",  "Solyc05g007920",                     # Zhong 2022b / Deng 2025
    "StDMP",   "Soltu.DM.05G005100",                            # J. Zhang et al. 2022
    "ClDMP3",  "Cla97C06G121370",                               # Chen et al. 2023
    "CsaV3_1G028660",                                           # Yin et al. 2024 (cucumber CsDMP)
    "OsDMP1",  "OsDMP3", "LOC_Os08g01530", "LOC_Os01g29240",    # Liang et al. 2025
    "GhDMPa",  "GhDMPd", "Gh_A11G3045", "Gh_D11G0735",          # Long et al. 2024
    "MtDMP",   "Medtr7g010890", "Medtr5g044580",                # N. Wang et al. 2022
    "ZmDMP",   "Zm00001d044822",                                # Zhong et al. 2019
    "BoDMP",   "LOC106333617", "LOC106333853",                  # Zhao et al. 2022
    "BjuDMP",  "BjuA04g10430S", "BjuA03g54090S",                # Chu et al. 2025
    "BjuB08g57390S", "BjuB01g27600S"
)

# Classify using ORIGINAL labels (before any renaming)
classify_tip_orig <- function(label) {
    if (grepl("^SMEL5_|^Smel_|^SMELG|^SmelDMP", label)) return("SmelDMP")
    for (pat in haploid_inducer_patterns) {
        if (grepl(pat, label, fixed = TRUE)) return("Haploid Inducer")
    }
    if (grepl(args$outgroup_pattern, label)) return("Outgroup")
    return("Other DMP Ortholog")
}

orig_categories <- sapply(tree$tip.label, classify_tip_orig, USE.NAMES = FALSE)

# ======================== Label Cleaning ========================

# Canonical HI-registry alias map — locus / accession token -> gene symbol.
# Source: II_INPUTS/DMP_HI_registry.tsv (validated DMP haploid-inducer genes).
# Substring match: any tip label containing one of these keys is replaced
# with the corresponding gene symbol so phylogeny figures show readable names.
hi_alias_map <- c(
    "BjuA04g10430S"      = "BjuDMP1",
    "BjuA03g54090S"      = "BjuDMP2",
    "BjuB08g57390S"      = "BjuDMP3",
    "BjuB01g27600S"      = "BjuDMP4",
    "LOC106333617"       = "BoDMP9",
    "LOC106333853"       = "BoDMP9",
    "Cla97C06G121370"    = "ClDMP3",
    "CsaV3_1G028660"     = "CsDMP",
    "Glyma.18G097400"    = "GmDMP1",
    "Glyma.18G098300"    = "GmDMP2",
    "Gh_A11G3045"        = "GhDMPa",
    "Gh_D11G0735"        = "GhDMPd",
    "LOC_Os08g01530"     = "OsDMP1",
    "LOC_Os01g29240"     = "OsDMP3",
    "Solyc05g007920"     = "SlDMP3",
    "Soltu.DM.05G005100" = "StDMP",
    "Medtr7g010890"      = "MtDMP8",
    "Medtr5g044580"      = "MtDMP9",
    "Zm00001d044822"     = "ZmDMP",
    # Arabidopsis thaliana DMP paralogs (AT locus -> gene symbol)
    "AT1G09157"          = "AtDMP8",
    "AT3G02430"          = "AtDMP5",
    "AT3G21520"          = "AtDMP1",
    "AT3G21550"          = "AtDMP2",
    "AT4G18425"          = "AtDMP4",
    "AT4G24310"          = "AtDMP3",
    "AT4G28485"          = "AtDMP7",
    "AT5G27370"          = "AtDMP10",
    "AT5G39650"          = "AtDMP9",
    "AT5G46090"          = "AtDMP6"
)

# Cosmetic shortening: removes purely-redundant accession suffixes from
# already-shortened tip labels (e.g., "BoDMP_LOC106333617" -> "BoDMP9",
# "AtDMP8_AtDMP9" -> "AtDMP8/9", "GmDMP1_2" -> "GmDMP1.2"). Applied to BOTH
# short and full forms so the rendered label collapses to one clean name when
# the only difference is cosmetic.
shorten_label <- function(lab) {
    # 1. Registry alias: substring match on canonical locus tokens
    for (key in names(hi_alias_map)) {
        if (grepl(key, lab, fixed = TRUE)) return(hi_alias_map[[key]])
    }
    # 2. Strip trailing NCBI LOC accession appended to a gene symbol
    lab <- sub("_LOC[0-9]+$", "", lab)
    # 3. Collapse dual-gene labels: "AtDMP8_AtDMP9" -> "AtDMP8/9"
    m <- regmatches(lab,
        regexec("^([A-Za-z]+)DMP([0-9]+)_\\1DMP([0-9]+)$", lab))[[1]]
    if (length(m) == 4L) {
        return(sprintf("%sDMP%s/%s", m[2], m[3], m[4]))
    }
    # 4. Paralog disambiguator: trailing "_N" (single digit) -> ".N".
    # Restricted to one digit so this never collapses NCBI accession suffixes
    # like "XM_008658777" (which would otherwise become "XM.008658777").
    lab <- sub("_([0-9])$", ".\\1", lab)
    lab
}

# Short name: gene symbol only — strips accessions, AT loci, Glyma IDs, and metadata.
short_gene_name <- function(lab) {
    lab <- shorten_label(lab)
    lab <- sub("_CDS_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_(XP|XM|XR|WP)_[0-9]+\\.?[0-9]*.*$", "", lab)
    lab <- sub("_AT[0-9]+G[0-9]+\\.?[0-9]*.*$", "", lab)
    lab <- sub("_Glyma\\.[0-9A-Z]+\\.?[0-9]*.*$", "", lab)
    lab <- sub("_G[0-9]+(_Outgroup)?$", "", lab)
    # Trailing citation strip. Includes '.' so years converted from "_2022"
    # to ".2022" by the earlier digit-to-dot rule in shorten_label() also peel.
    lab <- sub("_(Wang|Zhong|Deng|Liu|Yang|Lin)[._A-Za-z0-9]*$", "", lab)
    return(lab)
}

# Full name: gene symbol + accession/locus ID — strips only metadata notes.
full_gene_name <- function(lab) {
    lab <- shorten_label(lab)
    lab <- sub("_CDS_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_G[0-9]+(_Outgroup)?$", "", lab)
    # Trailing citation strip. Includes '.' so years converted from "_2022"
    # to ".2022" by the earlier digit-to-dot rule in shorten_label() also peel.
    lab <- sub("_(Wang|Zhong|Deng|Liu|Yang|Lin)[._A-Za-z0-9]*$", "", lab)
    return(lab)
}

# Build new tip labels in "Short_Name (Full Name)" format.
new_labels <- sapply(tree$tip.label, function(lab) {
    if (tolower(args$label_style) == "replace") {
        # SmelDMP: short = manuscript name; full = clean SMEL5 locus ID
        for (smel_id in names(smeldmp_name_map)) {
            if (grepl(smel_id, lab, fixed = TRUE)) {
                short <- smeldmp_name_map[[smel_id]]
                full  <- sub("_CDS_[0-9]+-[0-9]+$", "", lab)
                full  <- sub("_[0-9]+-[0-9]+$", "", full)
                return(paste0(short, " (", full, ")"))
            }
        }
        # Others: "GeneName (GeneName_Accession)"
        short <- short_gene_name(lab)
        full  <- full_gene_name(lab)
        if (short == full) return(short)
        return(paste0(short, " (", full, ")"))
    } else {
        # Legacy append mode: keep original label, add [SmelDMP...] bracket
        for (smel_id in names(smeldmp_name_map)) {
            if (grepl(smel_id, lab, fixed = TRUE)) {
                return(paste0(lab, "  [", smeldmp_name_map[[smel_id]], "]"))
            }
        }
        return(lab)
    }
}, USE.NAMES = FALSE)

tree$tip.label <- new_labels

# Map original category calls to new (cleaned) labels
tip_categories <- setNames(orig_categories, new_labels)

# ======================== Visual Encoding ========================

cat_levels <- c("SmelDMP", "Haploid Inducer", "Other DMP Ortholog", "Outgroup")

# Category colors — driven by TOML/args
category_colors <- c(
    "SmelDMP"            = args$color_smeldmp,
    "Haploid Inducer"    = args$color_haploid,
    "Other DMP Ortholog" = args$color_ortholog,
    "Outgroup"           = args$color_outgroup
)

# Distinct shapes per category — aids colorblind readers
# 18 = filled diamond, 16 = filled circle, 15 = filled square, 17 = filled triangle
category_shapes <- c(
    "SmelDMP"            = 18,
    "Haploid Inducer"    = 16,
    "Other DMP Ortholog" = 15,
    "Outgroup"           = 17
)

category_labels_map <- c(
    "SmelDMP"            = "SmelDMP (eggplant)",
    "Haploid Inducer"    = "Haploid inducer (validated)",
    "Other DMP Ortholog" = "Other DMP ortholog",
    "Outgroup"           = "Outgroup"
)

# ======================== Build Tree Plot ========================

show_bootstrap     <- tolower(args$show_bootstrap) == "true" && has_bootstrap
highlight_eggplant <- tolower(args$highlight_eggplant) == "true"
use_dot_bootstrap  <- tolower(args$bootstrap_style) == "dots"

# Base tree — branch color from config (default grey20: softer than black, crisp at print resolution)
p <- ggtree(tree, layout = args$layout, linewidth = args$branch_width, color = args$branch_color)

if (args$layout == "fan") {
    p <- ggtree(tree, layout = "fan", linewidth = args$branch_width,
                open.angle = args$open_angle, color = args$branch_color)
}

# --- Inject category + fontface into tree data (before adding geom layers) ---
p$data$category <- ifelse(
    p$data$isTip,
    tip_categories[p$data$label],
    NA_character_
)
p$data$category <- factor(p$data$category, levels = cat_levels)

# Fontface: outgroup = italic; eggplant / haploid = bold; others = plain
p$data$tip_face <- ifelse(
    p$data$isTip & p$data$category == "Outgroup",
    "italic",
    ifelse(
        p$data$isTip & p$data$category %in% c("SmelDMP", "Haploid Inducer"),
        "bold",
        "plain"
    )
)

# --- Bootstrap dot preprocessing (compute bs_val for all nodes) ---
if (show_bootstrap && use_dot_bootstrap) {
    p$data$bs_val <- sapply(p$data$label, extract_last_boot, USE.NAMES = FALSE)
}

# --- Bootstrap dots at internal nodes (dot mode) ---
# Two tiers: high (BS >= 95) = deep red; medium (>= threshold, < 95) = orange
if (show_bootstrap && use_dot_bootstrap) {
    threshold <- args$bootstrap_threshold
    p <- p +
        geom_nodepoint(
            aes(subset = (!isTip & !is.na(bs_val) & bs_val >= 95)),
            shape = 21, fill = args$node_color_high, color = "white",
            size = args$node_size_high, stroke = 0.35
        ) +
        geom_nodepoint(
            aes(subset = (!isTip & !is.na(bs_val) & bs_val >= threshold & bs_val < 95)),
            shape = 21, fill = args$node_color_medium, color = "white",
            size = args$node_size_medium, stroke = 0.35
        )
}

# --- Tip points (shape + color aes — drives a combined, merged legend) ---
if (highlight_eggplant) {
    p <- p +
        geom_tippoint(
            aes(color = category, shape = category),
            size = args$tip_point_size,
            show.legend = TRUE
        ) +
        geom_tiplab(
            aes(color = category, fontface = tip_face),
            size   = args$tip_label_size,
            offset = args$tip_label_offset,
            show.legend = FALSE
        ) +
        scale_color_manual(
            name   = "Gene Category",
            values = category_colors,
            labels = category_labels_map,
            na.value = "grey50",
            drop = TRUE
        ) +
        scale_shape_manual(
            name   = "Gene Category",
            values = category_shapes,
            labels = category_labels_map,
            drop = TRUE
        ) +
        guides(
            color = guide_legend(
                override.aes   = list(size = 4.5),
                title.position = "top",
                title.hjust    = 0.5
            ),
            shape = guide_legend(
                title.position = "top",
                title.hjust    = 0.5
            )
        )
} else {
    p <- p + geom_tiplab(size = args$tip_label_size, offset = args$tip_label_offset, color = "grey20")
}

# --- Bootstrap text labels (legacy text mode only) ---
if (show_bootstrap && !use_dot_bootstrap) {
    threshold <- args$bootstrap_threshold
    get_last <- function(x) {
        parts <- strsplit(as.character(x), "/")[[1]]
        parts[length(parts)]
    }
    p <- p + geom_nodelab(
        aes(label = ifelse(
            !is.na(suppressWarnings(as.numeric(sapply(label, get_last)))) &
            suppressWarnings(as.numeric(sapply(label, get_last))) >= threshold,
            sapply(label, get_last),
            ""
        )),
        size  = args$bootstrap_label_size,
        color = args$bootstrap_color,
        hjust = 1.5,
        vjust = -0.4
    )
}

# --- Scale bar ---
p <- p + geom_treescale(
    fontsize  = args$treescale_fontsize,
    offset    = args$treescale_offset,
    color     = args$treescale_color
)

# --- Title ---
if (nchar(args$title) > 0) {
    p <- p + ggtitle(args$title) +
        theme(plot.title = element_text(
            size = 13, face = "bold", hjust = 0.5, color = "grey15",
            margin = margin(b = 4)
        ))
}

# --- Center tree + labels (rectangular only) ---
if (args$layout == "rectangular") {
    max_label_len <- max(nchar(tree$tip.label))
    x_range   <- range(p$data$x, na.rm = TRUE)
    x_span    <- diff(x_range)
    label_ext <- max_label_len * args$tip_label_size * 0.006 * x_span
    root_pad  <- label_ext * 0.5   # left padding centers tree+labels composite
    p <- p + coord_cartesian(
        xlim = c(x_range[1] - root_pad, x_range[2] + label_ext * 1.1),
        clip = "off"
    )
}

# --- Phylogenetic parameters + inline bootstrap key in caption ---
param_parts <- c()
if (nchar(args$phylo_software)  > 0) param_parts <- c(param_parts, args$phylo_software)
if (nchar(args$sequence_type)   > 0) param_parts <- c(param_parts, args$sequence_type)
if (nchar(args$phylo_model)     > 0) param_parts <- c(param_parts, paste0("Model: ", args$phylo_model))
if (nchar(args$phylo_bootstrap) > 0) param_parts <- c(param_parts, paste0("Bootstrap: ", args$phylo_bootstrap))

# Inline node-dot key (dot mode only)
if (show_bootstrap && use_dot_bootstrap) {
    param_parts <- c(
        param_parts,
        sprintf("Node support:  \u25cf BS \u2265 95 (red)   \u25cf BS \u2265 %d (orange)",
                args$bootstrap_threshold)
    )
}

if (length(param_parts) > 0) {
    p <- p + labs(caption = paste(param_parts, collapse = "  |  "))
}

# --- Theme ---
p <- p + theme(
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(15, 15, 10, 15, unit = "pt"),
    plot.caption      = element_text(
        size   = 8,
        color  = "grey40",
        hjust  = 0,
        face   = "italic",
        margin = margin(t = 8)
    ),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 9),
    legend.key         = element_blank(),
    legend.background  = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.margin      = margin(6, 10, 6, 10),
    legend.box.margin  = margin(5, 0, 5, 0)
)

# ======================== Save Figure ========================

cat(sprintf("Saving: %s (%dx%d @ %d DPI)\n",
            args$output, as.integer(args$width), as.integer(args$height), args$dpi))

ggsave(
    filename  = args$output,
    plot      = p,
    width     = args$width,
    height    = args$height,
    dpi       = args$dpi,
    units     = "in",
    bg        = "white",
    limitsize = FALSE
)

cat(sprintf("Done: %s (%.1f MB)\n", args$output, file.info(args$output)$size / 1e6))
