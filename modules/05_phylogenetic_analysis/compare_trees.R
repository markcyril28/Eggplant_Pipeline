#!/usr/bin/env Rscript
# ============================================================================
# compare_trees.R — Tanglegram comparison of two phylogenetic trees
# ============================================================================
# Suppress fontconfig warnings on headless Linux / WSL2
options(bitmapType = "cairo")
Sys.setenv(FONTCONFIG_PATH = "/etc/fonts")
# Produces a side-by-side tanglegram comparing two Newick trees (e.g. IQ-TREE2
# vs RAxML) with connecting lines between shared tips and a normalized
# Robinson-Foulds (RF) distance annotation.
#
# Primary method:  phytools::cophyloplot() — proper tanglegram with connecting
#                  lines and ladderized tip ordering.
# Fallback method: ggtree + patchwork side-by-side panel when phytools is not
#                  available.
#
# Usage:
#   Rscript compare_trees.R \
#       --tree1 path/to/BASENAME_IQTREE2.treefile \
#       --tree2 path/to/BASENAME_RAXML.raxml.support \
#       --label1 "IQ-TREE2" \
#       --label2 "RAxML" \
#       --output path/to/comparison.png \
#       [--dpi 600] [--width 14] [--height 10]
# ============================================================================

suppressPackageStartupMessages({
    library(argparse)
    library(ape)
})

# ======================== Argument Parsing ========================

parser <- ArgumentParser(
    description = "Tanglegram comparison of two phylogenetic trees"
)
parser$add_argument("--tree1",
    required = TRUE,
    help = "First tree file (e.g. IQ-TREE2 .treefile)")
parser$add_argument("--tree2",
    required = TRUE,
    help = "Second tree file (e.g. RAxML .raxml.support or .raxml.bestTree)")
parser$add_argument("--label1",
    default = "IQ-TREE2",
    help = "Label for tree1 [IQ-TREE2]")
parser$add_argument("--label2",
    default = "RAxML",
    help = "Label for tree2 [RAxML]")
parser$add_argument("--output",
    required = TRUE,
    help = "Output PNG file path")
parser$add_argument("--dpi",
    type = "integer", default = 600,
    help = "Figure resolution in DPI [600]")
parser$add_argument("--width",
    type = "double", default = 14,
    help = "Figure width in inches [14]")
parser$add_argument("--height",
    type = "double", default = 10,
    help = "Figure height in inches [10]")
parser$add_argument("--root-outgroup",
    default = "false",
    help = "Root trees on outgroup? [false]")
parser$add_argument("--outgroup-pattern",
    default = NULL,
    help = "Regex pattern for outgroup tip label(s)")
parser$add_argument("--exclude-tips",
    default = "",
    help = "Comma-separated regex patterns of tips to drop before comparison [none]")
parser$add_argument("--highlight-eggplant",
    default = "true",
    help = "Colour eggplant tips red? [true]")
parser$add_argument("--eggplant-pattern",
    default = "^SMEL5_|^Smel_|^SMELG",
    help = "Regex to identify eggplant tips [^SMEL5_|^Smel_|^SMELG]")
parser$add_argument("--show-bootstrap",
    default = "true",
    help = "Show bootstrap values on nodes [true]")
parser$add_argument("--bootstrap-threshold",
    type = "double", default = 70,
    help = "Minimum bootstrap value to display as node label [70]")
parser$add_argument("--bootstrap-style",
    default = "text",
    help = "Bootstrap display: text (default) or dots (colored circles) [text]")
parser$add_argument("--bootstrap-label-size",
    type = "double", default = 3.0,
    help = "Bootstrap text label font size [3.0]")
parser$add_argument("--bootstrap-color",
    default = "#D32F2F",
    help = "Bootstrap text label color [#D32F2F]")
parser$add_argument("--node-color-high",
    default = "#B71C1C",
    help = "Node dot color for BS >= 95 — dots mode only [#B71C1C]")
parser$add_argument("--node-color-medium",
    default = "#E65100",
    help = "Node dot color for moderate support — dots mode only [#E65100]")
parser$add_argument("--node-size-high",
    type = "double", default = 3.0,
    help = "Node dot size for high-confidence nodes — dots mode only [3.0]")
parser$add_argument("--node-size-medium",
    type = "double", default = 2.0,
    help = "Node dot size for moderate support — dots mode only [2.0]")
parser$add_argument("--tip-label-size",
    type = "double", default = 4.5,
    help = "Tip label font size [4.5]")
parser$add_argument("--tip-label-offset",
    type = "double", default = 0.005,
    help = "Horizontal offset of tip labels [0.005]")
parser$add_argument("--tip-point-size",
    type = "double", default = 2.5,
    help = "Tip point symbol size [2.5]")
parser$add_argument("--branch-width",
    type = "double", default = 0.8,
    help = "Branch line width [0.8]")
parser$add_argument("--branch-color",
    default = "grey20",
    help = "Branch line color [grey20]")
parser$add_argument("--seq-type",
    default = "",
    help = "Sequence type label: Nucleotide, Protein, or empty [empty]")
parser$add_argument("--label-style",
    default = "replace",
    help = "Tip label style: replace (clean short names) or append [replace]")
parser$add_argument("--treescale-fontsize",
    type = "double", default = 3.0,
    help = "Scale bar font size [3.0]")
parser$add_argument("--treescale-offset",
    type = "double", default = 0.3,
    help = "Scale bar vertical offset [0.3]")
parser$add_argument("--treescale-color",
    default = "grey35",
    help = "Scale bar color [grey35]")
parser$add_argument("--color-smeldmp",
    default = "#4527A0",
    help = "Tip color for SmelDMP eggplant genes [#4527A0]")
parser$add_argument("--color-haploid",
    default = "#B71C1C",
    help = "Tip color for validated haploid inducers [#B71C1C]")
parser$add_argument("--color-ortholog",
    default = "#37474F",
    help = "Tip color for other DMP orthologs [#37474F]")
parser$add_argument("--color-outgroup",
    default = "#78909C",
    help = "Tip color for outgroup taxa [#78909C]")

args <- parser$parse_args()

# Build title prefix from seq-type if provided
seq_prefix <- if (nchar(trimws(args$seq_type)) > 0) {
    paste0(args$seq_type, ": ")
} else {
    ""
}

# ======================== Helpers ========================

to_bool <- function(x) tolower(trimws(as.character(x))) == "true"

read_tree_safe <- function(path) {
    if (!file.exists(path)) stop(paste("Tree file not found:", path))
    tree <- tryCatch(ape::read.tree(path), error = function(e) NULL)
    if (is.null(tree)) stop(paste("Cannot parse Newick tree:", path))
    if (!inherits(tree, "phylo")) stop(paste("Unexpected object type for:", path))
    tree
}

root_tree_safe <- function(tree, pattern) {
    if (is.null(pattern) || nchar(trimws(pattern)) == 0) return(tree)
    hits <- grep(pattern, tree$tip.label, value = TRUE)
    if (length(hits) == 0) {
        message("[WARN] No outgroup tips matched pattern: ", pattern, " — skipping re-root")
        return(tree)
    }
    rooted <- tryCatch({
        if (length(hits) == 1L) {
            ape::root(tree, outgroup = hits, resolve.root = TRUE)
        } else {
            # Multiple outgroup taxa: root on their MRCA to avoid ape version ambiguity
            mrca_node <- ape::getMRCA(tree, hits)
            ape::root(tree, node = mrca_node, resolve.root = TRUE)
        }
    }, error = function(e) {
        message("[WARN] Could not root tree: ", e$message, " — keeping original root")
        NULL
    })
    if (!is.null(rooted)) {
        tree <- rooted
        # Clear root node label — no bootstrap support is defined at the new root
        if (!is.null(tree$node.label) && length(tree$node.label) > 0L) {
            tree$node.label[1L] <- ""
        }
    }
    tree
}

# ======================== Package Detection ========================

has_phytools  <- requireNamespace("phytools",  quietly = TRUE)
has_ggtree    <- requireNamespace("ggtree",    quietly = TRUE)
has_ggplot2   <- requireNamespace("ggplot2",   quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

cat(sprintf("[INFO] phytools=%s  ggtree=%s  patchwork=%s\n",
    has_phytools, has_ggtree, has_patchwork))

# ======================== Read Trees ========================

cat("[INFO] Reading tree1:", args$tree1, "\n")
tree1 <- read_tree_safe(args$tree1)
cat("[INFO] Reading tree2:", args$tree2, "\n")
tree2 <- read_tree_safe(args$tree2)

# Exclude tips (comma-separated regex patterns)
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
    tree1 <- drop_matched(tree1, args$label1)
    tree2 <- drop_matched(tree2, args$label2)
}

# Optionally re-root
if (to_bool(args$root_outgroup) && !is.null(args$outgroup_pattern)) {
    tree1 <- root_tree_safe(tree1, args$outgroup_pattern)
    tree2 <- root_tree_safe(tree2, args$outgroup_pattern)
}

# ======================== Tip Statistics ========================

n_t1       <- length(tree1$tip.label)
n_t2       <- length(tree2$tip.label)
shared_tips <- intersect(tree1$tip.label, tree2$tip.label)
n_shared   <- length(shared_tips)

cat(sprintf("[INFO] %s: %d tips | %s: %d tips | Shared: %d tips\n",
    args$label1, n_t1, args$label2, n_t2, n_shared))

if (n_shared < 3) {
    stop(sprintf(
        "Only %d shared tip(s) between trees — need at least 3 for a meaningful comparison.",
        n_shared
    ))
}

# ======================== Auto-scale Canvas ========================
# Side-by-side panels need more vertical space than single-tree figures.
# Scale height so each tip gets ~0.6 in; scale width so long tip labels fit.
# Only grows — never shrinks below the configured value.

n_max <- max(n_t1, n_t2)

auto_h <- n_max * 0.6
if (auto_h > args$height) {
    cat(sprintf("[INFO] Auto-scaled height %.1f → %.1f inches for %d tips\n",
                args$height, auto_h, n_max))
    args$height <- auto_h
}

auto_w <- 10 + n_max * 0.45   # ~10 in base + 0.45 in per tip for label room
if (auto_w > args$width) {
    cat(sprintf("[INFO] Auto-scaled width %.1f → %.1f inches for %d tips\n",
                args$width, auto_w, n_max))
    args$width <- auto_w
}

# Scale tip label size down for dense tanglegrams (panels are narrower than single-tree figures)
if (n_max > 25 && args$tip_label_size == 4.5) {
    args$tip_label_size <- max(1.8, 4.2 - (n_max / 25))
    cat(sprintf("[INFO] Auto-scaled tip_label_size to %.1f for %d tips\n",
                args$tip_label_size, n_max))
}

# ======================== RF Distance ========================

t1_pruned <- ape::drop.tip(tree1, setdiff(tree1$tip.label, shared_tips))
t2_pruned <- ape::drop.tip(tree2, setdiff(tree2$tip.label, shared_tips))

rf_dist <- tryCatch({
    # getFromNamespace accesses RF.dist regardless of export status.
    # Some ape builds (< 5.7) have it internal-only; this handles both cases.
    rf_fn <- utils::getFromNamespace("RF.dist", "ape")
    rf_fn(t1_pruned, t2_pruned, normalize = TRUE)
}, error = function(e) {
    message("[WARN] Could not compute RF distance: ", e$message)
    NA_real_
})

rf_label <- if (!is.na(rf_dist)) {
    sprintf("Normalized RF distance: %.3f  (0=identical, 1=maximally different)", rf_dist)
} else {
    "RF distance: N/A"
}
cat("[INFO]", rf_label, "\n")

# ======================== Output Directory ========================

out_dir <- dirname(args$output)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================== METHOD 1: phytools cophyloplot ========================
# Produces a proper tanglegram with connecting lines between shared tip labels.

if (has_phytools) {
    library(phytools)

    # Build association matrix: each shared tip linked to itself
    assoc <- cbind(shared_tips, shared_tips)

    # Ladderize for a clean layout
    t1l <- ape::ladderize(tree1)
    t2l <- ape::ladderize(tree2)

    # Tip label colours — 4-category scheme consistent with ggtree fallback
    classify_tip_color <- function(lab) {
        if (grepl(args$eggplant_pattern, lab)) return(args$color_smeldmp)
        for (pat in c("AtDMP8", "AtDMP9", "GmDMP")) {
            if (grepl(pat, lab, fixed = TRUE)) return(args$color_haploid)
        }
        if (!is.null(args$outgroup_pattern) && grepl(args$outgroup_pattern, lab)) {
            return(args$color_outgroup)
        }
        args$color_ortholog
    }

    if (to_bool(args$highlight_eggplant)) {
        tip_col1 <- sapply(t1l$tip.label, classify_tip_color, USE.NAMES = FALSE)
        tip_col2 <- sapply(t2l$tip.label, classify_tip_color, USE.NAMES = FALSE)
    } else {
        tip_col1 <- rep("black", length(t1l$tip.label))
        tip_col2 <- rep("black", length(t2l$tip.label))
    }

    # Auto-scale spacing to accommodate tip count
    panel_space <- max(28, ceiling(n_shared * 0.85))
    # Normalize ggplot2-scale tip_label_size (~4.5) to phytools fsize (~0.5) scale
    tip_fsize_base <- args$tip_label_size * 0.11
    tip_fsize      <- max(0.20, min(tip_fsize_base, 18 / max(n_t1, n_t2)))

    png(filename = args$output,
        width    = args$width,
        height   = args$height,
        units    = "in",
        res      = args$dpi,
        bg       = "white")

    # Outer margins: bottom=4 for caption, top=3 for title
    par(oma = c(4.5, 0, 3.5, 0), mar = c(0, 0, 0, 0))

    tryCatch(
        phytools::cophyloplot(
            t1l, t2l,
            assoc       = assoc,
            space       = panel_space,
            length.line = 1,
            gap         = 2,
            lwd         = args$branch_width,
            col         = "grey55",
            link.type   = "curved",
            link.lwd    = 0.7,
            link.lty    = 1,
            fsize       = tip_fsize,
            ftype       = "i",
            tip.color   = list(tip_col1, tip_col2)
        ),
        error = function(e) {
            message("[WARN] tip.color as list not supported; retrying without tip colours")
            phytools::cophyloplot(
                t1l, t2l,
                assoc       = assoc,
                space       = panel_space,
                length.line = 1,
                gap         = 2,
                lwd         = args$branch_width,
                col         = "grey55",
                link.type   = "curved",
                link.lwd    = 0.7,
                link.lty    = 1,
                fsize       = tip_fsize,
                ftype       = "i"
            )
        }
    )

    # Panel labels (left = tree1, right = tree2)
    mtext(args$label1, side = 3, line = 1.5, outer = FALSE,
          at  = par("usr")[1], adj = 0,
          font = 2, cex = 1.1)
    mtext(args$label2, side = 3, line = 1.5, outer = FALSE,
          at  = par("usr")[2], adj = 1,
          font = 2, cex = 1.1)

    # Title
    mtext(
        sprintf("%s%s vs %s — Tanglegram Comparison", seq_prefix, args$label1, args$label2),
        side = 3, outer = TRUE, line = 1.8,
        font = 2, cex  = 1.3
    )

    # RF distance subtitle
    mtext(rf_label, side = 3, outer = TRUE, line = 0.4, cex = 0.95, col = "grey30")

    # Footer
    mtext(
        sprintf("Tips: %s=%d | %s=%d | Shared=%d | Generated %s",
            args$label1, n_t1, args$label2, n_t2, n_shared,
            format(Sys.time(), "%Y-%m-%d")),
        side = 1, outer = TRUE, line = 2.5, cex = 0.78, col = "grey40"
    )

    dev.off()
    cat("[INFO] Tanglegram (phytools) saved:", args$output, "\n")
    quit(status = 0)
}

# ======================== METHOD 2: ggtree + patchwork ========================
# Side-by-side ggtree panels — fallback when phytools is unavailable.

if (!has_ggtree || !has_ggplot2) {
    stop(paste(
        "Neither phytools nor ggtree/ggplot2 is available.",
        "Install phytools or ggtree to generate comparison figures."
    ))
}

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggtree)
})
if (has_patchwork) suppressPackageStartupMessages(library(patchwork))

# ======================== Label Cleaning & Category Classification ========================
# Matches render_tree.R conventions for consistent visual identity across all tree figures.

smeldmp_name_map <- c(
    "SMEL5_01g008730.1" = "SmelDMP01.730",
    "SMEL5_01g026030.1" = "SmelDMP01.990",
    "SMEL5_02g013320.1" = "SmelDMP02",
    "SMEL5_04g005390.1" = "SmelDMP04",
    "SMEL5_10g003660.1" = "SmelDMP10.560",
    "SMEL5_10g017610.1" = "SmelDMP10.200",
    "SMEL5_12g005350.1" = "SmelDMP12"
)

haploid_inducer_patterns <- c(
    "AtDMP8",   # Arabidopsis thaliana — Zhong et al. 2020
    "AtDMP9",   # Arabidopsis thaliana — Zhong et al. 2020
    "GmDMP",    # Glycine max — Zhong et al. 2024 (GmDMP1 + GmDMP2)
    "NtDMP",    # Nicotiana tabacum — X. Zhang et al. 2022 (NtDMP1-3)
    "SlDMP8"    # Solanum lycopersicum (SlDMP8-like, Solyc05g007920) — Deng et al. 2025
    # NOTE: "CsDMP9"/XP_006482605 is Citrus sinensis — NOT Cucumis sativus (Yin 2024 HI = CsaV3_1G028660).
)

classify_tip_orig <- function(label) {
    if (grepl("^SMEL5_|^Smel_|^SMELG", label)) return("SmelDMP")
    for (pat in haploid_inducer_patterns) {
        if (grepl(pat, label, fixed = TRUE)) return("Haploid Inducer")
    }
    if (!is.null(args$outgroup_pattern) && grepl(args$outgroup_pattern, label)) {
        return("Outgroup")
    }
    "Other DMP Ortholog"
}

short_gene_name <- function(lab) {
    lab <- sub("_CDS_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_(XP|XM|XR|WP)_[0-9]+\\.?[0-9]*.*$", "", lab)
    lab <- sub("_AT[0-9]+G[0-9]+\\.?[0-9]*.*$", "", lab)
    lab <- sub("_Glyma\\.[0-9A-Z]+\\.?[0-9]*.*$", "", lab)
    lab <- sub("_G[0-9]+(_Outgroup)?$", "", lab)
    lab <- sub("_(Wang|Zhong|Deng|Liu|Yang|Lin)[_A-Za-z0-9]*$", "", lab)
    lab
}

full_gene_name <- function(lab) {
    lab <- sub("_CDS_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_[0-9]+-[0-9]+$", "", lab)
    lab <- sub("_G[0-9]+(_Outgroup)?$", "", lab)
    lab <- sub("_(Wang|Zhong|Deng|Liu|Yang|Lin)[_A-Za-z0-9]*$", "", lab)
    lab
}

clean_labels <- function(labels) {
    sapply(labels, function(lab) {
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
            paste0(short, " (", full, ")")
        } else {
            for (smel_id in names(smeldmp_name_map)) {
                if (grepl(smel_id, lab, fixed = TRUE)) {
                    return(paste0(lab, "  [", smeldmp_name_map[[smel_id]], "]"))
                }
            }
            lab
        }
    }, USE.NAMES = FALSE)
}

# Classify on ORIGINAL labels, then clean
orig_cats1 <- sapply(tree1$tip.label, classify_tip_orig, USE.NAMES = FALSE)
orig_cats2 <- sapply(tree2$tip.label, classify_tip_orig, USE.NAMES = FALSE)

new_labels1 <- clean_labels(tree1$tip.label)
new_labels2 <- clean_labels(tree2$tip.label)

tip_cats1 <- setNames(orig_cats1, new_labels1)
tip_cats2 <- setNames(orig_cats2, new_labels2)

tree1$tip.label <- new_labels1
tree2$tip.label <- new_labels2

# ======================== Tip Order Alignment ========================
# Rotate tree1 (left) to match tree2 (right) tip ordering so corresponding
# tips appear at the same vertical position for easy visual comparison.

tree2 <- ape::ladderize(tree2)

# Get tree2's tip order from a temporary ggtree layout (bottom -> top by y)
p2_tmp <- ggtree(tree2, layout = "rectangular")
t2_tip_df <- p2_tmp$data[p2_tmp$data$isTip, ]
tree2_tip_order <- t2_tip_df$label[order(t2_tip_df$y)]

tree1 <- tryCatch(
    ape::rotateConstr(tree1, tree2_tip_order),
    error = function(e) {
        message("[WARN] Could not align tip orders: ", e$message, " -- ladderizing instead")
        ape::ladderize(tree1)
    }
)

# ======================== Visual Encoding ========================

cat_levels <- c("SmelDMP", "Haploid Inducer", "Other DMP Ortholog", "Outgroup")

category_colors <- c(
    "SmelDMP"            = args$color_smeldmp,
    "Haploid Inducer"    = args$color_haploid,
    "Other DMP Ortholog" = args$color_ortholog,
    "Outgroup"           = args$color_outgroup
)

category_shapes <- c(
    "SmelDMP"            = 18,   # filled diamond
    "Haploid Inducer"    = 16,   # filled circle
    "Other DMP Ortholog" = 15,   # filled square
    "Outgroup"           = 17    # filled triangle
)

category_labels_map <- c(
    "SmelDMP"            = "SmelDMP (eggplant)",
    "Haploid Inducer"    = "Haploid inducer (validated)",
    "Other DMP Ortholog" = "Other DMP ortholog",
    "Outgroup"           = "Outgroup"
)

# ======================== Panel Builder ========================

parse_bootstrap <- function(labels, sep = "/") {
    sapply(labels, function(x) {
        if (is.na(x) || nchar(trimws(x)) == 0) return(NA_real_)
        parts <- strsplit(x, sep, fixed = TRUE)[[1]]
        suppressWarnings(as.numeric(trimws(parts[length(parts)])))
    })
}

make_ggtree_panel <- function(tree, tip_cats, label, side = "left") {
    tip_size <- args$tip_label_size
    is_right <- (side == "right")

    # Pre-process bootstrap node labels before building ggtree object
    has_node_labels <- !is.null(tree$node.label) && length(tree$node.label) > 0
    if (has_node_labels) {
        bs_vals <- parse_bootstrap(tree$node.label)
        tree$node.label <- ifelse(
            !is.na(bs_vals) & bs_vals >= args$bootstrap_threshold,
            as.character(round(bs_vals)),
            NA_character_
        )
    }

    # Build base tree — branch width and color from config for consistency with single-tree figures
    p <- ggtree(tree, layout = "rectangular",
                linewidth = args$branch_width, color = args$branch_color) +
        theme_tree2() +
        labs(title = label) +
        theme(
            plot.title  = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.margin = margin(5, 15, 5, 15)
        )

    # Inject category and fontface into tree data
    p$data$category <- ifelse(
        p$data$isTip, tip_cats[p$data$label], NA_character_
    )
    p$data$category <- factor(p$data$category, levels = cat_levels)

    p$data$tip_face <- ifelse(
        p$data$isTip & p$data$category == "Outgroup", "italic",
        ifelse(
            p$data$isTip & p$data$category %in% c("SmelDMP", "Haploid Inducer"),
            "bold", "plain"
        )
    )

    # Right panel: last character at branch tip (hjust=1, offset=0)
    # Left panel:  first character near branch tip (hjust=0, offset from config)
    tip_hjust  <- if (is_right) 1 else 0
    tip_offset <- if (is_right) 0 else args$tip_label_offset

    if (to_bool(args$highlight_eggplant)) {
        p <- p +
            geom_tippoint(
                aes(color = category, shape = category),
                size = args$tip_point_size,
                show.legend = !is_right   # single legend from left panel only
            ) +
            geom_tiplab(
                aes(color = category, fontface = tip_face),
                size   = tip_size,
                hjust  = tip_hjust,
                offset = tip_offset,
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
                    override.aes   = list(size = 4),
                    title.position = "top",
                    title.hjust    = 0.5
                ),
                shape = guide_legend(
                    title.position = "top",
                    title.hjust    = 0.5
                )
            )
    } else {
        p <- p + geom_tiplab(
            size = tip_size, hjust = tip_hjust,
            offset = tip_offset, color = "grey20"
        )
    }

    # Bootstrap node labels — gated by show_bootstrap, style from config
    show_bs  <- tolower(args$show_bootstrap) == "true"
    use_dots <- tolower(args$bootstrap_style) == "dots"

    if (show_bs && has_node_labels && use_dots) {
        # Dots mode: two-tier colored circles
        p$data$bs_val <- sapply(p$data$label, function(x) {
            if (is.na(x) || nchar(trimws(as.character(x))) == 0) return(NA_real_)
            suppressWarnings(as.numeric(trimws(x)))
        }, USE.NAMES = FALSE)

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
    } else if (show_bs && has_node_labels) {
        # Text mode: numeric labels from config-driven size and color
        p <- p + geom_nodelab(
            aes(label = label),
            size  = args$bootstrap_label_size,
            hjust = -0.15,
            color = args$bootstrap_color,
            na.rm = TRUE
        )
    }

    # Scale bar — size and color from config
    p <- p + geom_treescale(
        fontsize = args$treescale_fontsize,
        offset   = args$treescale_offset,
        color    = args$treescale_color
    )

    # --- Center tree within panel ---
    # Add root-side padding to balance the label-side extension
    x_range   <- range(p$data$x, na.rm = TRUE)
    x_span    <- diff(x_range)
    max_chars <- max(nchar(tree$tip.label))
    label_ext <- max_chars * tip_size * 0.006 * x_span
    root_pad  <- label_ext   # symmetric padding -> centered tree + labels

    if (is_right) {
        p <- p +
            scale_x_reverse() +
            coord_cartesian(
                xlim = c(x_range[1] - root_pad, x_range[2] + label_ext * 1.1),
                clip = "off"
            )
    } else {
        p <- p +
            coord_cartesian(
                xlim = c(x_range[1] - root_pad, x_range[2] + label_ext * 1.1),
                clip = "off"
            )
    }

    p
}

# ======================== Build Panels ========================

p1 <- make_ggtree_panel(tree1, tip_cats1, args$label1, side = "left")
p2 <- make_ggtree_panel(tree2, tip_cats2, args$label2, side = "right")

caption_text <- sprintf(
    "Tips: %s=%d | %s=%d | Shared=%d | %s | Generated %s",
    args$label1, n_t1, args$label2, n_t2, n_shared,
    rf_label,
    format(Sys.time(), "%Y-%m-%d")
)

if (has_patchwork) {
    combined <- p1 + p2 +
        patchwork::plot_layout(widths = c(1, 1), guides = "collect") +
        patchwork::plot_annotation(
            title    = sprintf("%s%s vs %s — Phylogenetic Tree Comparison",
                               seq_prefix, args$label1, args$label2),
            subtitle = rf_label,
            caption  = caption_text,
            theme = theme(
                plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey35"),
                plot.caption  = element_text(size = 7,  hjust = 0.5, color = "grey50"),
                legend.position = "bottom"
            )
        )

    ggsave(
        args$output,
        plot      = combined,
        width     = args$width,
        height    = args$height,
        dpi       = args$dpi,
        units     = "in",
        limitsize = FALSE,
        bg        = "white"
    )
    cat("[INFO] Comparison (ggtree+patchwork) saved:", args$output, "\n")
} else {
    # No patchwork — save each panel separately
    base_path <- tools::file_path_sans_ext(args$output)
    out1 <- paste0(base_path, "_", args$label1, ".png")
    out2 <- paste0(base_path, "_", args$label2, ".png")
    half_w <- args$width / 2

    ggsave(out1, p1, width = half_w, height = args$height,
           dpi = args$dpi, units = "in", limitsize = FALSE, bg = "white")
    ggsave(out2, p2, width = half_w, height = args$height,
           dpi = args$dpi, units = "in", limitsize = FALSE, bg = "white")

    cat("[WARN] patchwork not available — saved two separate panels:\n")
    cat("  ", out1, "\n")
    cat("  ", out2, "\n")
}
