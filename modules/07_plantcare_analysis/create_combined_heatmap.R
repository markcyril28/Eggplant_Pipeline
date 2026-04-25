#!/usr/bin/env Rscript
# Create combined motif-function heatmap using ComplexHeatmap

suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(grid)
    library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript create_combined_heatmap.R <matrix.tsv> <output_dir> <dpi> [color_palette] [row_font] [col_font] [cell_font] [cell_size] [label_font] [cell_border_color] [cell_border_width] [column_rotation] [legend_height] [gene_label]")
}

matrix_file <- args[1]
output_dir <- args[2]
dpi <- as.numeric(args[3])
color_palette_str <- if (length(args) >= 4 && nchar(args[4]) > 0) args[4] else ""
row_font <- if (length(args) >= 5 && nchar(args[5]) > 0) as.numeric(args[5]) else 10
col_font <- if (length(args) >= 6 && nchar(args[6]) > 0) as.numeric(args[6]) else 8
cell_font <- if (length(args) >= 7 && nchar(args[7]) > 0) as.numeric(args[7]) else 8
cell_size <- if (length(args) >= 8 && nchar(args[8]) > 0) as.numeric(args[8]) else 1
label_font <- if (length(args) >= 9 && nchar(args[9]) > 0) as.numeric(args[9]) else 14
cell_border_color <- if (length(args) >= 10 && nchar(args[10]) > 0) args[10] else "white"
cell_border_width <- if (length(args) >= 11 && nchar(args[11]) > 0) as.numeric(args[11]) else 0.5
column_rotation <- if (length(args) >= 12 && nchar(args[12]) > 0) as.numeric(args[12]) else 22.5
legend_height <- if (length(args) >= 13 && nchar(args[13]) > 0) as.numeric(args[13]) else 4
gene_label <- if (length(args) >= 14 && nchar(args[14]) > 0) args[14] else "Genes of Interest"

safe_read_tsv <- function(file_path, header = FALSE) {
    if (!file.exists(file_path)) {
        message(paste0("Warning: file not found, skipping heatmap: ", file_path))
        return(NULL)
    }
    file_size <- file.info(file_path)$size
    if (is.na(file_size) || file_size == 0) {
        message(paste0("Warning: file is empty, skipping heatmap: ", file_path))
        return(NULL)
    }
    read.table(file_path, sep = "\t", header = header, stringsAsFactors = FALSE, check.names = FALSE)
}

# Load data
mat_df <- safe_read_tsv(matrix_file, header = TRUE)
if (is.null(mat_df)) {
    quit(save = "no", status = 0)
}
mat <- mat_df
if (ncol(mat) == 0) {
    message("Warning: matrix has zero columns, skipping heatmap.")
    quit(save = "no", status = 0)
}

metadata_file <- file.path(dirname(matrix_file), "column_metadata.txt")
metadata <- safe_read_tsv(metadata_file, header = FALSE)
if (is.null(metadata) || nrow(metadata) == 0) {
    message("Warning: column metadata is empty, skipping heatmap.")
    quit(save = "no", status = 0)
}
colnames(metadata) <- c("index", "name", "type", "function")

# Load function groups
function_groups_file <- file.path(dirname(matrix_file), "function_groups.txt")
function_groups <- safe_read_tsv(function_groups_file, header = FALSE)
if (is.null(function_groups) || nrow(function_groups) == 0) {
    message("Warning: function groups are empty, skipping heatmap.")
    quit(save = "no", status = 0)
}
colnames(function_groups) <- c("func_name", "count")

# Convert to matrix
mat <- as.matrix(mat)

# Color palette
if (nchar(color_palette_str) > 0) {
    color_palette <- strsplit(color_palette_str, ",")[[1]]
} else {
    color_palette <- c("#d9e9c8ff","#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
}
col_fun <- colorRamp2(
    seq(min(mat), max(mat), length.out = length(color_palette)),
    color_palette
)

# Identify function columns for spacing - create column splits
col_split <- rep(NA, ncol(mat))
current_pos <- 1
split_id <- 1

for (i in 1:nrow(function_groups)) {
    func_name <- function_groups$func_name[i]
    motif_count <- function_groups$count[i]
    
    for (j in 1:motif_count) {
        if (current_pos <= ncol(mat)) {
            col_split[current_pos] <- split_id
            current_pos <- current_pos + 1
        }
    }
    split_id <- split_id + 1
}
col_split <- factor(col_split, levels = unique(col_split))

# Create labels for functions (one label per group)
function_labels <- rep("", ncol(mat))
current_pos <- 1
for (i in 1:nrow(function_groups)) {
    func_name <- function_groups$func_name[i]
    motif_count <- function_groups$count[i]
    
    # Put function name at the midpoint of its motifs
    mid_pos <- current_pos + floor(motif_count / 2)
    if (mid_pos <= ncol(mat)) {
        function_labels[mid_pos] <- func_name
    }
    current_pos <- current_pos + motif_count
}

# Create custom column labels annotation
col_anno <- HeatmapAnnotation(
    function_label = anno_text(
        function_labels,
        rot = 22.5,
        location = unit(0.01, "npc"),
        just = "left",
        gp = gpar(fontsize = 9),
        height = unit(1, "cm")
    ),
    show_annotation_name = FALSE,
    which = "column"
)

# Column names are all motifs now
# (no need to filter)

# Create heatmap
ht <- Heatmap(
    mat,
    name = "Frequency",
    col = col_fun,
    
    # Cell config
    rect_gp = gpar(col = cell_border_color, lwd = cell_border_width),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(as.character(mat[i, j]), x, y, gp = gpar(fontsize = cell_font))
    },
    show_heatmap_legend = TRUE,
    width = unit(ncol(mat) * cell_size, "cm"),
    height = unit(nrow(mat) * cell_size, "cm"),
    
    # Row config
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_font), #fontface = "bold"),
    row_title = gene_label,
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = label_font, fontface = "bold"),
    cluster_rows = FALSE,
    
    # Column config
    column_names_side = "bottom",
    column_names_rot = column_rotation,
    column_names_gp = gpar(fontsize = col_font),
    column_names_centered = FALSE,
    column_title = "Elements",
    column_title_gp = gpar(fontsize = label_font, fontface = "bold"),
    column_title_side = "bottom",
    cluster_columns = FALSE,
    column_split = col_split,
    column_gap = unit(2, "mm"),
    show_column_dend = FALSE,
    
    # Add function label annotation at top
    top_annotation = col_anno,
    
    # Legend config
    heatmap_legend_param = list(
        title = "Frequency",
        title_position = "topcenter",
        legend_direction = "vertical",
        legend_height = unit(legend_height, "cm"),
        at = seq(min(mat), max(mat), length.out = 5),
        labels_gp = gpar(fontsize = 10),
        title_gp = gpar(fontsize = 12, fontface = "bold")
    )
)

# Save with high DPI and padding
output_file <- file.path(output_dir, "combined_heatmap.png")
png(output_file, width = 10 + ncol(mat) * 0.2, height = 6 + nrow(mat) * 0.4, 
    units = "in", res = dpi)
# padding order: bottom, left, top, right
draw(ht, padding = unit(c(0, 0.25, 0.05, 0.05), "inches"))
dev.off()

cat(paste0("Heatmap saved: ", output_file, "\n"))
