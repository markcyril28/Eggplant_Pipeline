#!/usr/bin/env Rscript
# ============================================================================
# Combined Bootstrap Tree — IQ-TREE2 UFBoot + RAxML-NG standard BS
# ============================================================================
# Renders a publication-quality phylogenetic tree (IQ-TREE2 topology) with
# combined bootstrap support from both IQ-TREE2 (UFBoot) and RAxML-NG at
# every internal node, formatted as "UFBoot/RAX-BS" (e.g., "98/95").
#
# Clade matching: for each internal node in the IQ-TREE2 tree, the set of
# descendant tips is computed via BFS on the edge matrix.  The same tip-set
# key is looked up in the RAxML tree, giving the paired bootstrap value.
# This is robust to any difference in node indexing order between the two
# programs.
#
# Usage:
#   Rscript render_combined_bootstrap.R \
#       --tree-iqtree <stem_IQTREE2.treefile> \
#       --tree-raxml  <stem_RAXML.raxml.support> \
#       --output      <out.png> \
#       [--bootstrap-threshold 70] \
#       [--combined-style fraction]   # fraction="98/95"; dual="IQ:98 RAX:95"
#       [--combined-sep /]            # separator between IQ and RAX values
#       [... all shared render_tree.R visual args ...]
# ============================================================================
options(bitmapType = "cairo")
Sys.setenv(FONTCONFIG_PATH = "/etc/fonts")

suppressPackageStartupMessages({
    library(ggtree)
    library(treeio)
    library(ggplot2)
    library(ape)
    library(argparse)
})

# ======================== Argument Parsing ========================

parser <- ArgumentParser(
    description = "Render combined IQ-TREE2 UFBoot + RAxML-NG bootstrap tree figure")

# --- Required ---
parser$add_argument("--tree-iqtree",  required = TRUE,
    help = "IQ-TREE2 .treefile (node labels as SH-aLRT/UFBoot)")
parser$add_argument("--tree-raxml",   required = TRUE,
    help = "RAxML-NG .raxml.support file (standard bootstrap integers)")
parser$add_argument("--output",       required = TRUE,
    help = "Output PNG file path")

# --- Combined-bootstrap display ---
parser$add_argument("--combined-style", default = "fraction",
    help = "Node label format: fraction='98/95' (default) or dual='IQ:98 RAX:95'")
parser$add_argument("--combined-sep",   default = "/",
    help = "Separator between IQ-TREE2 and RAxML values (default: /)")

# --- Bootstrap ---
parser$add_argument("--bootstrap-threshold", type = "integer", default = 70L,
    help = "Only show labels where max(UFBoot, RAX-BS) >= threshold (default: 70)")
parser$add_argument("--bootstrap-label-size", type = "double", default = 3.0,
    help = "Node label font size (default: 3.0)")
parser$add_argument("--bootstrap-color", default = "#D32F2F",
    help = "Node label color (default: #D32F2F)")
parser$add_argument("--show-bootstrap", default = "true",
    help = "Show bootstrap values on nodes (default: true)")
parser$add_argument("--bootstrap-style", default = "text",
    help = "Bootstrap display: text (combined labels, default) or dots (colored circles)")
parser$add_argument("--node-color-high", default = "#B71C1C",
    help = "Node dot color for BS >= 95 — dots mode only (default: #B71C1C deep red)")
parser$add_argument("--node-color-medium", default = "#E65100",
    help = "Node dot color for moderate support — dots mode only (default: #E65100 deep orange)")
parser$add_argument("--node-size-high", type = "double", default = 3.0,
    help = "Node dot size for high-confidence nodes — dots mode only (default: 3.0)")
parser$add_argument("--node-size-medium", type = "double", default = 2.0,
    help = "Node dot size for moderate support — dots mode only (default: 2.0)")

# --- Tip labels ---
parser$add_argument("--label-style", default = "replace",
    help = "Tip label style: replace (clean short names, default) or append")
parser$add_argument("--tip-label-size",  type = "double", default = 4.5,
    help = "Tip label font size (default: 4.5)")
parser$add_argument("--tip-point-size",  type = "double", default = 2.5,
    help = "Tip point symbol size (default: 2.5)")
parser$add_argument("--tip-label-offset", type = "double", default = 0.005,
    help = "Horizontal offset of tip labels (default: 0.005)")

# --- Canvas ---
parser$add_argument("--layout",  default = "rectangular",
    help = "Tree layout: rectangular (default), circular, fan, daylight, equal_angle")
parser$add_argument("--open-angle", type = "double", default = 15.0,
    help = "Fan layout open angle in degrees (default: 15)")
parser$add_argument("--width",   type = "double", default = 14.0,
    help = "Figure width in inches (default: 14)")
parser$add_argument("--height",  type = "double", default = 10.0,
    help = "Figure height in inches (default: 10)")
parser$add_argument("--dpi",     type = "integer", default = 600L,
    help = "Figure resolution in DPI (default: 600)")
parser$add_argument("--xlim-expand", type = "double", default = 1.0,
    help = "Multiplier for x-axis expansion beyond tree tips; >1 = more room for labels / wider tree, <1 = tighter (default: 1.0)")

# --- Branches ---
parser$add_argument("--branch-width", type = "double", default = 0.8,
    help = "Branch line width (default: 0.8)")
parser$add_argument("--branch-color", default = "grey20",
    help = "Branch line color (default: grey20)")

# --- Category colors ---
parser$add_argument("--highlight-eggplant", default = "true",
    help = "Color-code tips by gene category (default: true)")
parser$add_argument("--outgroup-pattern",   default = "Outgroup",
    help = "Regex to identify outgroup tips (default: Outgroup)")
parser$add_argument("--root-outgroup",      default = "false",
    help = "Root tree on outgroup before rendering (default: false)")
parser$add_argument("--exclude-tips",       default = "",
    help = "Comma-separated regex patterns of tips to drop before rendering (default: none)")
parser$add_argument("--color-smeldmp",  default = "#4527A0",
    help = "Tip color for SmelDMP eggplant genes (default: #4527A0 deep purple)")
parser$add_argument("--color-haploid",  default = "#B71C1C",
    help = "Tip color for validated haploid inducers (default: #B71C1C deep red)")
parser$add_argument("--color-ortholog", default = "#37474F",
    help = "Tip color for other DMP orthologs (default: #37474F blue-gray)")
parser$add_argument("--color-outgroup", default = "#78909C",
    help = "Tip color for outgroup taxa (default: #78909C gray)")

# --- Scale bar ---
parser$add_argument("--treescale-fontsize", type = "double", default = 3.0)
parser$add_argument("--treescale-offset",   type = "double", default = 0.3)
parser$add_argument("--treescale-color",    default = "grey35")

# --- Caption metadata ---
parser$add_argument("--phylo-model",   default = "",
    help = "Substitution model string for figure caption")
parser$add_argument("--sequence-type", default = "",
    help = "Sequence type (AA or NT) for figure caption")
parser$add_argument("--title",         default = "",
    help = "Plot title (auto-generated from sequence-type and model if blank)")
parser$add_argument("--topology-source", default = "IQ-TREE2",
    help = "Topology tree source label for subtitles and captions")
parser$add_argument("--support-source",  default = "RAxML-NG",
    help = "Support tree source label for subtitles and captions")
parser$add_argument("--tree-third", default = "",
    help = "Optional third support tree Newick file")
parser$add_argument("--third-source", default = "",
    help = "Third support tree source label")

args <- parser$parse_args()

# ======================== Input Validation ========================

if (!file.exists(args$tree_iqtree)) {
    stop(sprintf("IQ-TREE2 file not found: %s", args$tree_iqtree))
}
if (!file.exists(args$tree_raxml)) {
    stop(sprintf("RAxML file not found: %s", args$tree_raxml))
}
if (nchar(args$tree_third) > 0 && !file.exists(args$tree_third)) {
    stop(sprintf("Third support tree file not found: %s", args$tree_third))
}

has_third_tree <- nchar(args$tree_third) > 0
if (has_third_tree && nchar(args$third_source) == 0) {
    args$third_source <- "Third"
}

output_dir <- dirname(args$output)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

valid_layouts <- c("rectangular", "circular", "fan", "daylight", "equal_angle")
if (!(args$layout %in% valid_layouts)) {
    stop(sprintf("Invalid layout '%s'. Choose from: %s",
                 args$layout, paste(valid_layouts, collapse = ", ")))
}

# ======================== Read Trees ========================

cat(sprintf("Reading topology tree (%s): %s\n", args$topology_source, args$tree_iqtree))
iq_tree  <- ape::read.tree(args$tree_iqtree)
cat(sprintf("Reading support tree (%s):  %s\n", args$support_source, args$tree_raxml))
rax_tree <- ape::read.tree(args$tree_raxml)
if (has_third_tree) {
    cat(sprintf("Reading third support tree (%s): %s\n", args$third_source, args$tree_third))
    third_tree <- ape::read.tree(args$tree_third)
}

cat(sprintf("  %s: %d tips, %d internal nodes\n",
            args$topology_source, length(iq_tree$tip.label), iq_tree$Nnode))
cat(sprintf("  %s: %d tips, %d internal nodes\n",
            args$support_source, length(rax_tree$tip.label), rax_tree$Nnode))
if (has_third_tree) {
    cat(sprintf("  %s: %d tips, %d internal nodes\n",
                args$third_source, length(third_tree$tip.label), third_tree$Nnode))
}

common_tips <- intersect(iq_tree$tip.label, rax_tree$tip.label)
if (length(common_tips) < 3) {
    stop(sprintf("Too few common tips (%d) between IQ-TREE2 and RAxML trees.",
                 length(common_tips)))
}
if (length(common_tips) < length(iq_tree$tip.label)) {
    cat(sprintf("  Warning: %d %s tips absent from %s tree\n",
                length(iq_tree$tip.label) - length(common_tips),
                args$topology_source, args$support_source))
}
if (has_third_tree) {
    common_tips_third <- intersect(iq_tree$tip.label, third_tree$tip.label)
    if (length(common_tips_third) < 3) {
        stop(sprintf("Too few common tips (%d) between topology and third support trees.",
                     length(common_tips_third)))
    }
    if (length(common_tips_third) < length(iq_tree$tip.label)) {
        cat(sprintf("  Warning: %d %s tips absent from %s tree\n",
                    length(iq_tree$tip.label) - length(common_tips_third),
                    args$topology_source, args$third_source))
    }
}

# ======================== Exclude Tips (optional) ========================

if (nchar(args$exclude_tips) > 0) {
    patterns <- trimws(strsplit(args$exclude_tips, ",")[[1]])
    drop_matched <- function(tree, tree_name) {
        tips_to_drop <- c()
        for (pat in patterns) {
            if (nchar(pat) > 0) {
                matched <- grep(pat, tree$tip.label, value = TRUE)
                tips_to_drop <- c(tips_to_drop, matched)
            }
        }
        tips_to_drop <- unique(tips_to_drop)
        if (length(tips_to_drop) > 0) {
            cat(sprintf("[INFO] Excluding from %s: %s\n",
                        tree_name, paste(tips_to_drop, collapse = ", ")))
            tree <- ape::drop.tip(tree, tips_to_drop)
        }
        tree
    }
    iq_tree  <- drop_matched(iq_tree,  args$topology_source)
    rax_tree <- drop_matched(rax_tree, args$support_source)
    if (has_third_tree) {
        third_tree <- drop_matched(third_tree, args$third_source)
    }

    # Recompute common tips after exclusion
    common_tips <- intersect(iq_tree$tip.label, rax_tree$tip.label)
    if (length(common_tips) < 3) {
        stop(sprintf("Too few common tips (%d) after exclusion.", length(common_tips)))
    }
    cat(sprintf("  Tips after exclusion: %s=%d, %s=%d\n",
                args$topology_source, length(iq_tree$tip.label),
                args$support_source, length(rax_tree$tip.label)))
    if (has_third_tree) {
        cat(sprintf("  Tips after exclusion: %s=%d\n",
                    args$third_source, length(third_tree$tip.label)))
    }
}

# ======================== Rooting (Optional) ========================

root_on_outgroup <- function(tree, outgroup_pattern) {
    outgroup_tips <- grep(outgroup_pattern, tree$tip.label, value = TRUE)
    if (length(outgroup_tips) == 0) {
        cat(sprintf("  Warning: no tips matching outgroup pattern '%s'\n",
                    outgroup_pattern))
        return(tree)
    }
    cat(sprintf("  Rooting on %d outgroup taxa matching '%s'\n",
                length(outgroup_tips), outgroup_pattern))
    if (length(outgroup_tips) == 1L) {
        tree <- ape::root(tree, outgroup = outgroup_tips, resolve.root = TRUE)
    } else {
        # Guard: if outgroups span the current root (MRCA == root node),
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
    # Clear the root node label so it is not rendered as a spurious support value
    if (!is.null(tree$node.label) && length(tree$node.label) > 0L) {
        tree$node.label[1L] <- ""
    }
    tree
}

if (tolower(args$root_outgroup) == "true" && nchar(args$outgroup_pattern) > 0) {
    cat(sprintf("Rooting %s tree...\n", args$topology_source))
    iq_tree  <- root_on_outgroup(iq_tree,  args$outgroup_pattern)
    cat(sprintf("Rooting %s tree...\n", args$support_source))
    rax_tree <- root_on_outgroup(rax_tree, args$outgroup_pattern)
    if (has_third_tree) {
        cat(sprintf("Rooting %s tree...\n", args$third_source))
        third_tree <- root_on_outgroup(third_tree, args$outgroup_pattern)
    }
}

n_tips <- length(iq_tree$tip.label)

# ======================== Auto-size ========================

if (n_tips > 30) {
    auto_height <- max(10, n_tips * 0.4)
    if (auto_height > args$height) {
        args$height <- auto_height
        cat(sprintf("  Auto-scaled height to %.1f inches for %d tips\n",
                    args$height, n_tips))
    }
}
if (n_tips > 30) {
    auto_tip_size <- max(3.5, 5.5 - (n_tips / 100))
    if (auto_tip_size > args$tip_label_size) {
        args$tip_label_size <- auto_tip_size
        cat(sprintf("  Auto-scaled tip label size to %.1f\n", args$tip_label_size))
    }
}

# ======================== Clade Matching ========================
# BFS from a given node to collect all descendant tip labels.
# Works correctly for both rooted and unrooted trees as stored by ape
# (edges always run from parent to child in tree$edge).

get_subtree_tips <- function(tree, node) {
    n_tips_local <- length(tree$tip.label)
    visited <- logical(n_tips_local + tree$Nnode)
    visited[node] <- TRUE
    queue <- node
    tips  <- integer(0)
    while (length(queue) > 0L) {
        curr  <- queue[1L]; queue <- queue[-1L]
        children <- tree$edge[tree$edge[, 1L] == curr, 2L]
        for (ch in children) {
            if (!visited[ch]) {
                visited[ch] <- TRUE
                if (ch <= n_tips_local) tips <- c(tips, ch)
                else                    queue <- c(queue, ch)
            }
        }
    }
    sort(tree$tip.label[tips])
}

# Canonical bipartition key: always use the smaller half of the split.
# This ensures correct matching regardless of which side of an arbitrary
# root a node falls on (IQ-TREE2 and RAxML-NG may root Newick differently
# when root_outgroup = false).
canonical_key <- function(tips, all_tips) {
    complement <- setdiff(all_tips, tips)
    n_t <- length(tips)
    n_c <- length(complement)
    if (n_t < n_c) {
        return(paste(tips, collapse = "|"))
    } else if (n_t > n_c) {
        return(paste(complement, collapse = "|"))
    } else {
        # Equal size: pick lexicographically smaller
        key_t <- paste(tips, collapse = "|")
        key_c <- paste(complement, collapse = "|")
        return(min(key_t, key_c))
    }
}

# Build a lookup table: canonical tip-set key -> RAxML bootstrap value.
build_support_map <- function(support_tree) {
    n_tips_support <- length(support_tree$tip.label)
    all_tips       <- sort(support_tree$tip.label)
    support_map    <- list()
    for (i in seq_len(support_tree$Nnode)) {
        bs_raw <- support_tree$node.label[i]
        bs_num <- extract_support(bs_raw)
        if (is.na(bs_num)) next
        tips <- get_subtree_tips(support_tree, n_tips_support + i)
        if (length(tips) < 2L) next       # skip trivial (single-tip) clades
        key <- canonical_key(tips, all_tips)
        support_map[[key]] <- bs_num
    }
    support_map
}

# ======================== Combined Label Construction ========================

# Extract support from node label ("SH-aLRT/UFBoot", MEGA decimal supports,
# or plain integer).
extract_support <- function(x) {
    if (is.na(x) || nchar(trimws(as.character(x))) == 0L) return(NA_real_)
    parts <- strsplit(as.character(x), "/")[[1L]]
    raw <- parts[length(parts)]
    raw <- sub("\\[.*$", "", raw)
    raw <- gsub("[^0-9.eE+-]", "", raw)
    val <- suppressWarnings(as.numeric(raw))
    if (is.na(val)) return(NA_real_)
    if (val > 0 && val <= 1) val <- val * 100
    val
}

build_combined_labels <- function(iq_tree, rax_map, third_map, threshold, sep, style) {
    n_tips_iq <- length(iq_tree$tip.label)
    all_tips  <- sort(iq_tree$tip.label)
    n_nodes   <- iq_tree$Nnode
    labels    <- character(n_nodes)

    for (i in seq_len(n_nodes)) {
        ufboot <- extract_support(iq_tree$node.label[i])

        tips <- get_subtree_tips(iq_tree, n_tips_iq + i)
        if (length(tips) < 2L) { labels[i] <- ""; next }

        key    <- canonical_key(tips, all_tips)
        rax_bs <- rax_map[[key]]
        if (is.null(rax_bs)) rax_bs <- NA_real_
        third_bs <- third_map[[key]]
        if (is.null(third_bs)) third_bs <- NA_real_

        # Display only if at least one source meets the threshold
        ufboot_ok <- !is.na(ufboot) && ufboot >= threshold
        rax_ok    <- !is.na(rax_bs) && rax_bs >= threshold
        third_ok  <- !is.na(third_bs) && third_bs >= threshold
        if (!ufboot_ok && !rax_ok && !third_ok) { labels[i] <- ""; next }

        iq_str    <- if (!is.na(ufboot)) as.character(round(ufboot)) else "\u2013"
        rax_str   <- if (!is.na(rax_bs)) as.character(round(rax_bs)) else "\u2013"
        third_str <- if (!is.na(third_bs)) as.character(round(third_bs)) else "\u2013"

        labels[i] <- if (style == "dual") {
            label_parts <- c(
                paste0(args$topology_source, ":", iq_str),
                paste0(args$support_source, ":", rax_str)
            )
            if (length(third_map) > 0) {
                label_parts <- c(label_parts, paste0(args$third_source, ":", third_str))
            }
            paste(label_parts, collapse = " ")
        } else {
            if (length(third_map) > 0) {
                paste0(iq_str, sep, rax_str, sep, third_str)
            } else {
                paste0(iq_str, sep, rax_str)
            }
        }
    }
    labels
}

# Max of UFBoot and RAX-BS per internal node (drives dots mode tier assignment)
build_combined_max_bs <- function(iq_tree, rax_map, third_map) {
    n_tips_iq <- length(iq_tree$tip.label)
    all_tips  <- sort(iq_tree$tip.label)
    n_nodes   <- iq_tree$Nnode
    max_bs    <- rep(NA_real_, n_nodes)

    for (i in seq_len(n_nodes)) {
        ufboot <- extract_support(iq_tree$node.label[i])
        tips   <- get_subtree_tips(iq_tree, n_tips_iq + i)
        if (length(tips) < 2L) next
        key    <- canonical_key(tips, all_tips)
        rax_bs <- rax_map[[key]]
        if (is.null(rax_bs)) rax_bs <- NA_real_
        third_bs <- third_map[[key]]
        if (is.null(third_bs)) third_bs <- NA_real_
        vals   <- c(ufboot, rax_bs, third_bs)
        vals   <- vals[!is.na(vals)]
        if (length(vals) > 0L) max_bs[i] <- max(vals)
    }
    max_bs
}

cat(sprintf("Matching %s bootstrap values onto %s topology...\n",
            args$support_source, args$topology_source))
rax_map <- build_support_map(rax_tree)
cat(sprintf("  %s clades indexed: %d\n", args$support_source, length(rax_map)))
third_map <- list()
if (has_third_tree) {
    cat(sprintf("Matching %s bootstrap values onto %s topology...\n",
                args$third_source, args$topology_source))
    third_map <- build_support_map(third_tree)
    cat(sprintf("  %s clades indexed: %d\n", args$third_source, length(third_map)))
}

combined_labels <- build_combined_labels(
    iq_tree   = iq_tree,
    rax_map   = rax_map,
    third_map = third_map,
    threshold = args$bootstrap_threshold,
    sep       = args$combined_sep,
    style     = args$combined_style
)

combined_max_bs <- build_combined_max_bs(iq_tree, rax_map, third_map)

n_shown <- sum(nchar(combined_labels) > 0L)
has_dash <- any(grepl("\u2013", combined_labels[nchar(combined_labels) > 0L], fixed = TRUE))
cat(sprintf("  Combined labels to display: %d / %d internal nodes\n",
            n_shown, iq_tree$Nnode))
if (has_dash) cat("  Note: some nodes have \u2013 (one method lacks a matching clade)\n")

iq_tree$node.label <- combined_labels

# ======================== Gene Category Classification ========================

smeldmp_name_map <- c(
    "SMEL5_01g008730.1" = "SmelDMPv5_01.730",
    "SMEL5_01g026030.1" = "SmelDMPv5_01.030",
    "SMEL5_02g013320.1" = "SmelDMPv5_02.320",
    "SMEL5_04g005390.1" = "SmelDMPv5_04.390",
    "SMEL5_10g003660.1" = "SmelDMPv5_10.660",
    "SMEL5_10g017610.1" = "SmelDMPv5_10.610",
    "SMEL5_12g005350.1" = "SmelDMPv5_12.350"
)

# Validated HI loci from II_INPUTS/DMP_HI_registry.tsv (kept in sync with render_tree.R).
# NOTE: "CsDMP9"/XP_006482605 is Citrus sinensis — NOT Cucumis sativus.
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

classify_tip_orig <- function(label) {
    if (grepl("^SMEL5_|^Smel_|^SMELG", label)) return("SmelDMP")
    for (pat in haploid_inducer_patterns) {
        if (grepl(pat, label, fixed = TRUE)) return("Haploid Inducer")
    }
    if (grepl(args$outgroup_pattern, label)) return("Outgroup")
    return("Other DMP Ortholog")
}

orig_categories <- sapply(iq_tree$tip.label, classify_tip_orig, USE.NAMES = FALSE)

# ======================== Label Cleaning ========================

# Canonical HI-registry alias map — locus / accession token -> gene symbol.
# Source: II_INPUTS/DMP_HI_registry.tsv. Kept in sync with render_tree.R.
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

shorten_label <- function(lab) {
    for (key in names(hi_alias_map)) {
        if (grepl(key, lab, fixed = TRUE)) return(hi_alias_map[[key]])
    }
    lab <- sub("_LOC[0-9]+$", "", lab)
    m <- regmatches(lab,
        regexec("^([A-Za-z]+)DMP([0-9]+)_\\1DMP([0-9]+)$", lab))[[1]]
    if (length(m) == 4L) {
        return(sprintf("%sDMP%s/%s", m[2], m[3], m[4]))
    }
    # Single-digit only so NCBI accessions like "XM_008658777" pass through.
    lab <- sub("_([0-9])$", ".\\1", lab)
    lab
}

short_gene_name <- function(lab) {
    lab <- shorten_label(lab)
    lab <- sub("_CDS_[0-9]+-[0-9]+$",              "", lab)
    lab <- sub("_[0-9]+-[0-9]+$",                  "", lab)
    lab <- sub("_(XP|XM|XR|WP)_[0-9]+\\.?[0-9]*.*$", "", lab)
    lab <- sub("_AT[0-9]+G[0-9]+\\.?[0-9]*.*$",    "", lab)
    lab <- sub("_Glyma\\.[0-9A-Z]+\\.?[0-9]*.*$",  "", lab)
    lab <- sub("_G[0-9]+(_Outgroup)?$",             "", lab)
    # Trailing citation strip. Includes '.' in the suffix class so years that
    # were converted from "_2022" to ".2022" by the earlier digit-to-dot rule
    # in shorten_label() also get peeled off.
    lab <- sub("_(Wang|Zhong|Deng|Liu|Yang|Lin)[._A-Za-z0-9]*$", "", lab)
    return(lab)
}

full_gene_name <- function(lab) {
    lab <- shorten_label(lab)
    lab <- sub("_CDS_[0-9]+-[0-9]+$",  "", lab)
    lab <- sub("_[0-9]+-[0-9]+$",      "", lab)
    lab <- sub("_G[0-9]+(_Outgroup)?$", "", lab)
    # Trailing citation strip. Includes '.' in the suffix class so years that
    # were converted from "_2022" to ".2022" by the earlier digit-to-dot rule
    # in shorten_label() also get peeled off.
    lab <- sub("_(Wang|Zhong|Deng|Liu|Yang|Lin)[._A-Za-z0-9]*$", "", lab)
    return(lab)
}

new_labels <- sapply(iq_tree$tip.label, function(lab) {
    if (tolower(args$label_style) == "replace") {
        for (smel_id in names(smeldmp_name_map)) {
            if (grepl(smel_id, lab, fixed = TRUE)) {
                short <- smeldmp_name_map[[smel_id]]
                full  <- sub("_CDS_[0-9]+-[0-9]+$", "", lab)
                full  <- sub("_[0-9]+-[0-9]+$", "", full)
                return(paste0(short, " (", full, ")"))
            }
        }
        short <- short_gene_name(lab)
        full  <- full_gene_name(lab)
        if (short == full) return(short)
        return(paste0(short, " (", full, ")"))
    } else {
        for (smel_id in names(smeldmp_name_map)) {
            if (grepl(smel_id, lab, fixed = TRUE)) {
                return(paste0(lab, "  [", smeldmp_name_map[[smel_id]], "]"))
            }
        }
        return(lab)
    }
}, USE.NAMES = FALSE)

iq_tree$tip.label <- new_labels
tip_categories    <- setNames(orig_categories, new_labels)

# ======================== Visual Encoding ========================

cat_levels <- c("SmelDMP", "Haploid Inducer", "Other DMP Ortholog", "Outgroup")

category_colors <- c(
    "SmelDMP"            = args$color_smeldmp,
    "Haploid Inducer"    = args$color_haploid,
    "Other DMP Ortholog" = args$color_ortholog,
    "Outgroup"           = args$color_outgroup
)

# 18=filled diamond, 16=filled circle, 15=filled square, 17=filled triangle
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

# ======================== Pre-compute Haploid Tip Indices ========================

haploid_tip_idx <- which(orig_categories == "Haploid Inducer")

# ======================== Build Tree Plot ========================

highlight_eggplant <- tolower(args$highlight_eggplant) == "true"
show_bootstrap     <- tolower(args$show_bootstrap) == "true"
use_dot_bootstrap  <- tolower(args$bootstrap_style) == "dots"

p <- ggtree(iq_tree, layout = args$layout,
            linewidth = args$branch_width, color = args$branch_color)

if (args$layout == "fan") {
    p <- ggtree(iq_tree, layout = "fan",
                linewidth = args$branch_width,
                open.angle = args$open_angle,
                color = args$branch_color)
}

# ======================== Rotate: Haploid-Inducer-Rich Clade to Top ========================
# At every internal node, put the child clade containing MORE haploid inducers
# on top. Ties broken by higher haploid fraction, then by existing y-order.
# Uses ggtree::rotate() on the plot object so y-coordinates update immediately.

if (length(haploid_tip_idx) >= 1L) {
    n_t <- length(iq_tree$tip.label)
    root_nd <- n_t + 1L

    tips_under <- function(nd) {
        if (nd <= n_t) return(nd)
        stack <- nd; out <- integer(0)
        while (length(stack) > 0L) {
            v <- stack[1L]; stack <- stack[-1L]
            ch <- iq_tree$edge[iq_tree$edge[, 1] == v, 2]
            for (cc in ch) {
                if (cc <= n_t) out <- c(out, cc) else stack <- c(stack, cc)
            }
        }
        out
    }

    # Walk internal nodes root-to-tip so parent rotations settle before children
    internal_nodes <- unique(iq_tree$edge[, 1])
    depth <- sapply(internal_nodes, function(nd) {
        d <- 0L; cur <- nd
        while (cur != root_nd) {
            par <- iq_tree$edge[iq_tree$edge[, 2] == cur, 1]
            if (length(par) == 0L) break
            cur <- par; d <- d + 1L
        }
        d
    })
    internal_nodes <- internal_nodes[order(depth)]

    n_rotated <- 0L
    for (anc in internal_nodes) {
        ch <- iq_tree$edge[iq_tree$edge[, 1] == anc, 2]
        if (length(ch) != 2L) next
        c1_tips <- tips_under(ch[1L])
        c2_tips <- tips_under(ch[2L])
        c1_cnt <- sum(c1_tips %in% haploid_tip_idx)
        c2_cnt <- sum(c2_tips %in% haploid_tip_idx)
        if (c1_cnt == 0L && c2_cnt == 0L) next

        # Decide which child should be on top: more haploids wins;
        # ties broken by higher haploid fraction.
        top_is_c1 <- if (c1_cnt != c2_cnt) {
            c1_cnt > c2_cnt
        } else {
            (c1_cnt / length(c1_tips)) >= (c2_cnt / length(c2_tips))
        }
        top_tips <- if (top_is_c1) c1_tips else c2_tips
        bot_tips <- if (top_is_c1) c2_tips else c1_tips

        top_mean_y <- mean(p$data$y[p$data$node %in% top_tips], na.rm = TRUE)
        bot_mean_y <- mean(p$data$y[p$data$node %in% bot_tips], na.rm = TRUE)
        if (is.finite(top_mean_y) && is.finite(bot_mean_y) && top_mean_y < bot_mean_y) {
            p <- ggtree::rotate(p, anc)
            n_rotated <- n_rotated + 1L
        }
    }
    cat(sprintf("  Rotated %d node(s) to place haploid-inducer-rich clades at top\n", n_rotated))
}

# Inject category and fontface into tree data
p$data$category <- ifelse(
    p$data$isTip,
    tip_categories[p$data$label],
    NA_character_
)
p$data$category <- factor(p$data$category, levels = cat_levels)

p$data$tip_face <- ifelse(
    p$data$isTip & p$data$category == "Outgroup",
    "italic",
    ifelse(
        p$data$isTip & p$data$category %in% c("SmelDMP", "Haploid Inducer"),
        "bold",
        "plain"
    )
)

# --- Combined bootstrap at internal nodes ---
if (show_bootstrap && use_dot_bootstrap) {
    # Dots mode: mapped fill + size so a proper legend renders automatically
    bs_max_lookup <- rep(NA_real_, n_tips + iq_tree$Nnode)
    bs_max_lookup[(n_tips + 1):(n_tips + iq_tree$Nnode)] <- combined_max_bs
    p$data$bs_max <- bs_max_lookup[p$data$node]

    threshold       <- args$bootstrap_threshold
    tier_high_label <- "\u2265 95 (strong)"
    tier_mod_label  <- sprintf("\u2265 %d (moderate)", threshold)

    p$data$bs_tier <- NA_character_
    p$data$bs_tier[!p$data$isTip & !is.na(p$data$bs_max) & p$data$bs_max >= 95] <-
        tier_high_label
    p$data$bs_tier[!p$data$isTip & !is.na(p$data$bs_max) &
                   p$data$bs_max >= threshold & p$data$bs_max < 95] <- tier_mod_label
    p$data$bs_tier <- factor(p$data$bs_tier, levels = c(tier_high_label, tier_mod_label))

    node_support_data <- p$data[!is.na(p$data$bs_tier) & !p$data$isTip, , drop = FALSE]

    p <- p +
        geom_point(
            data        = node_support_data,
            aes(x = x, y = y, fill = bs_tier, size = bs_tier),
            shape       = 21,
            color       = "white",
            stroke      = 0.35,
            inherit.aes = FALSE,
            show.legend = TRUE
        ) +
        scale_fill_manual(
            name   = "Node Support\n(max UFBoot / RAxML-BS)",
            values = setNames(
                c(args$node_color_high, args$node_color_medium),
                c(tier_high_label, tier_mod_label)
            ),
            na.value = NA,
            drop = TRUE
        ) +
        scale_size_manual(
            name   = "Node Support\n(max UFBoot / RAxML-BS)",
            values = setNames(
                c(args$node_size_high, args$node_size_medium),
                c(tier_high_label, tier_mod_label)
            ),
            drop  = TRUE,
            guide = "none"
        ) +
        guides(
            fill = guide_legend(
                title.position = "top",
                title.hjust    = 0.5,
                override.aes   = list(size = 4.5, shape = 21, color = "white")
            )
        )
} else if (show_bootstrap) {
    # Text mode: combined "UFBoot/RAX-BS" labels (default)
    p <- p + geom_nodelab(
        aes(label = ifelse(!isTip & nchar(label) > 0L, label, "")),
        size  = args$bootstrap_label_size,
        color = args$bootstrap_color,
        hjust = 1.3,
        vjust = -0.4
    )
}

# --- Tip points and labels ---
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
            name     = "Gene Category",
            values   = category_colors,
            labels   = category_labels_map,
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
            color = guide_legend(override.aes   = list(size = 4.5),
                                 title.position = "top", title.hjust = 0.5),
            shape = guide_legend(title.position = "top", title.hjust = 0.5)
        )
} else {
    p <- p + geom_tiplab(size = args$tip_label_size, offset = args$tip_label_offset, color = "grey20")
}

# --- Scale bar ---
p <- p + geom_treescale(
    fontsize = args$treescale_fontsize,
    offset   = args$treescale_offset,
    color    = args$treescale_color
)

# --- Center tree + labels (rectangular only) ---
if (args$layout == "rectangular") {
    max_label_len <- max(nchar(iq_tree$tip.label))
    x_range   <- range(p$data$x, na.rm = TRUE)
    x_span    <- diff(x_range)
    label_ext <- max_label_len * args$tip_label_size * 0.006 * x_span * args$xlim_expand
    root_pad  <- label_ext * 0.5   # left padding centers tree+labels composite
    p <- p + coord_cartesian(
        xlim = c(x_range[1] - root_pad, x_range[2] + label_ext * 1.1),
        clip = "off"
    )
}

# --- Title, subtitle, and caption ---
sep_display <- args$combined_sep
support_label <- if (has_third_tree) {
    paste(args$topology_source, args$support_source, args$third_source, sep = sep_display)
} else {
    paste(args$topology_source, args$support_source, sep = sep_display)
}
support_sources_text <- if (has_third_tree) {
    sprintf("%s support + %s support + %s standard BS",
            args$topology_source, args$support_source, args$third_source)
} else {
    sprintf("%s support + %s standard BS", args$topology_source, args$support_source)
}

# Build plot title (auto-derive from sequence type + model when --title not given)
if (nchar(args$title) > 0) {
    plot_title <- args$title
} else {
    seq_label <- switch(trimws(toupper(args$sequence_type)),
        "AA"         = "Amino Acid",
        "PROTEIN"    = "Amino Acid",
        "NT"         = "Nucleotide",
        "NUCLEOTIDE" = "Nucleotide",
        if (nchar(args$sequence_type) > 0) args$sequence_type else "Phylogenetic"
    )
    model_str  <- if (nchar(args$phylo_model) > 0) paste0(" - ", args$phylo_model) else ""
    plot_title <- paste0(seq_label, " Phylogeny", model_str)
}

# Subtitle: always explains the dual-bootstrap format so it reads as a legend key
if (show_bootstrap && use_dot_bootstrap) {
    plot_subtitle <- sprintf(
        "Node dots: max(%s) - threshold \u2265 %d",
        support_label, args$bootstrap_threshold
    )
} else if (show_bootstrap) {
    dash_note <- if (has_dash) paste0("  (\u2013 = no matching clade)") else ""
    plot_subtitle <- sprintf(
        "Node labels: %s%s - values \u2265 %d shown",
        support_label, dash_note, args$bootstrap_threshold
    )
} else {
    plot_subtitle <- sprintf("%s topology with combined bootstrap support",
                             args$topology_source)
}

# Caption: concise parameter/provenance line
caption_parts <- c()
if (nchar(args$phylo_model) > 0) {
    caption_parts <- c(caption_parts, paste0("Model: ", args$phylo_model))
}
if (nchar(args$sequence_type) > 0) {
    caption_parts <- c(caption_parts, paste0("Sequence: ", args$sequence_type))
}
caption_parts <- c(caption_parts,
    sprintf("Topology: %s  |  Bootstrap sources: %s",
            args$topology_source, support_sources_text))

p <- p + labs(
    title    = plot_title,
    subtitle = plot_subtitle,
    caption  = paste(caption_parts, collapse = "  |  ")
)

# --- Theme ---
p <- p + theme(
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(15, 15, 15, 15, unit = "pt"),
    plot.title        = element_text(size = 13, face = "bold",    hjust = 0.5,
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 9,  color = "grey25", hjust = 0.5,
                                     face = "italic", margin = margin(b = 10)),
    plot.caption      = element_text(size = 7.5, color = "grey45", hjust = 0,
                                     face = "italic", margin = margin(t = 8)),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 9),
    legend.key         = element_blank(),
    legend.background  = element_rect(fill = "white", color = "grey80",
                                      linewidth = 0.3),
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

cat(sprintf("Done: %s (%.1f MB)\n",
            args$output, file.info(args$output)$size / 1e6))
