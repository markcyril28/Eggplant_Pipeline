#!/usr/bin/env Rscript
# ============================================================================
# Root a Newick tree on an outgroup tip and write the rooted tree out.
# ============================================================================
# Usage:
#   Rscript root_tree.R --input <newick_in> --output <newick_out> \
#                       --outgroup <regex_or_exact_label> \
#                       [--resolve true|false]
#
# Behaviour:
#   - The --outgroup argument is matched against tip labels as a regex first
#     (grep). If that produces exactly one tip, the tree is rooted on that
#     tip via ape::root(outgroup=...). If multiple tips match, the script
#     roots on their MRCA via ape::root(node=getMRCA(...)).
#   - Bootstrap support labels stored as node labels are preserved.
#   - --resolve defaults to true (resolve.root = TRUE in ape::root()).
# ============================================================================

suppressPackageStartupMessages({
    library(argparse)
    library(ape)
})

parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE, help = "Input Newick tree file (.treefile, .nwk, .raxml.support, .raxml.bestTree)")
parser$add_argument("--output", required = TRUE, help = "Output Newick file (rooted)")
parser$add_argument("--outgroup", required = TRUE, help = "Outgroup tip label or regex (e.g. 'PpDMP5-like' or '^Pp')")
parser$add_argument("--resolve", default = "true", help = "Pass resolve.root=TRUE to ape::root (default: true)")
args <- parser$parse_args()

resolve_root <- tolower(args$resolve) == "true"

cat(sprintf("[root_tree] input  : %s\n", args$input))
cat(sprintf("[root_tree] output : %s\n", args$output))
cat(sprintf("[root_tree] outgroup pattern: %s\n", args$outgroup))

tree <- ape::read.tree(args$input)
if (is.null(tree) || length(tree$tip.label) == 0) {
    stop(sprintf("Failed to read tree from %s", args$input))
}
cat(sprintf("[root_tree] tips: %d\n", length(tree$tip.label)))

matches <- grep(args$outgroup, tree$tip.label, value = TRUE)
if (length(matches) == 0) {
    stop(sprintf("No tip labels match outgroup pattern '%s'.\nFirst few tips: %s",
                 args$outgroup,
                 paste(head(tree$tip.label, 5), collapse = ", ")))
}
cat(sprintf("[root_tree] %d tip(s) matched outgroup pattern:\n", length(matches)))
for (m in matches) cat(sprintf("    %s\n", m))

if (length(matches) == 1L) {
    rooted <- ape::root(tree, outgroup = matches, resolve.root = resolve_root)
} else {
    mrca <- ape::getMRCA(tree, matches)
    if (is.null(mrca) || mrca <= length(tree$tip.label)) {
        cat(sprintf("[root_tree] outgroup MRCA could not be resolved; rooting on first match: %s\n",
                    matches[1]))
        rooted <- ape::root(tree, outgroup = matches[1], resolve.root = resolve_root)
    } else {
        cat(sprintf("[root_tree] rooting on MRCA node %d of %d outgroup tips\n", mrca, length(matches)))
        rooted <- ape::root(tree, node = mrca, resolve.root = resolve_root)
    }
}

dir.create(dirname(args$output), showWarnings = FALSE, recursive = TRUE)
ape::write.tree(rooted, file = args$output)
cat(sprintf("[root_tree] wrote: %s\n", args$output))
