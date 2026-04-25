#!/usr/bin/env Rscript
#
# Visualization script for PlantCARE matrix results using ComplexHeatmap
# Version 1: Heatmap by Motif Names
# Version 2: Heatmap by Function Descriptions
#

# Required libraries
suppressPackageStartupMessages({
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
        BiocManager::install("ComplexHeatmap", ask = FALSE)
    if (!requireNamespace("RColorBrewer", quietly = TRUE))
        install.packages("RColorBrewer", repos='http://cran.us.r-project.org')
    if (!requireNamespace("optparse", quietly = TRUE))
        install.packages("optparse", repos='http://cran.us.r-project.org')
    if (!requireNamespace("circlize", quietly = TRUE))
        install.packages("circlize", repos='http://cran.us.r-project.org')
    if (!requireNamespace("grid", quietly = TRUE))
        install.packages("grid", repos='http://cran.us.r-project.org')
})

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(RColorBrewer)
library(optparse)
library(grid)

# Function to create and save the heatmap
create_heatmap <- function(matrix_file, output_file, gene_label, column_label, row_font, col_font, label_font, cell_size, min_freq, cell_font_size, dpi, color_palette_str = NULL, column_rotation = 22.5, cell_border_color = "white", cell_border_width = 2, legend_height = 3) {
    # Read the matrix
    matrix_data <- read.csv(matrix_file, sep = "\t", row.names = 1, check.names = FALSE)

    # Filter columns (CAREs/Functions) where max frequency across all genes is below threshold
    col_max_freq <- apply(matrix_data, 2, max)
    cols_to_keep <- col_max_freq >= min_freq
    
    if (sum(cols_to_keep) == 0) {
        stop(sprintf("No columns have a maximum frequency >= %d. Cannot create heatmap.", min_freq))
    }
    
    removed_cols <- sum(!cols_to_keep)
    if (removed_cols > 0) {
        cat(sprintf("Filtered out %d column(s) with maximum frequency < %d\n", removed_cols, min_freq))
    }
    
    matrix_data <- matrix_data[, cols_to_keep, drop = FALSE]
    # Convert to numeric matrix; replace any NA (empty TSV cells) with 0
    matrix_data <- as.matrix(matrix_data)
    storage.mode(matrix_data) <- "numeric"
    matrix_data[is.na(matrix_data)] <- 0
    cat(sprintf("Remaining columns: %d\n", ncol(matrix_data)))
    cat(sprintf("Genes (rows): %d\n", nrow(matrix_data)))

    # Define the color palette (violet scale)
    if (!is.null(color_palette_str) && nchar(color_palette_str) > 0) {
        color_palette <- strsplit(color_palette_str, ",")[[1]]
    } else {
        color_palette <- c("#d9e9c8ff","#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
    }

    # Build a proper color-mapping function so ComplexHeatmap can handle any
    # data range, including cases where min == max or the range rounds to fewer
    # than 2 distinct integer breaks.
    data_min <- min(matrix_data)
    data_max <- max(matrix_data)
    if (data_max == data_min) data_max <- data_min + 1   # guard: avoid zero-range
    col_breaks <- seq(data_min, data_max, length.out = length(color_palette))
    col_fun <- colorRamp2(col_breaks, color_palette)

    # Legend break values: up to 5 evenly-spaced integers, always at least 2 distinct
    legend_at <- unique(round(seq(data_min, data_max, length.out = 5)))
    if (length(legend_at) < 2) legend_at <- c(floor(data_min), ceiling(data_max))

    # Create the heatmap
    ht <- Heatmap(as.matrix(matrix_data),
        name = "Frequency",
        col = col_fun,
        rect_gp = gpar(col = cell_border_color, lwd = cell_border_width),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = row_font), #, fontface = "bold"),
        column_names_gp = gpar(fontsize = col_font),
        row_names_side = "left",
        column_names_rot = column_rotation,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%d", matrix_data[i, j]), x, y, gp = gpar(fontsize = cell_font_size))
        },
        heatmap_legend_param = list(
            at = legend_at,
            legend_height = unit(legend_height, "cm"),
            title = "Frequency",
            title_position = "topcenter",
            direction = "vertical"
        ),
        width = unit(ncol(matrix_data) * cell_size, "cm"),
        height = unit(nrow(matrix_data) * cell_size, "cm")
    )

    # Calculate image dimensions (dynamic based on heatmap + 1 inch padding on all sides)
    width_cm <- ncol(matrix_data) * cell_size + 2 * 2.54 + 5  # Add extra space for labels
    height_cm <- nrow(matrix_data) * cell_size + 2 * 2.54 + 5
    width_in <- width_cm / 2.54
    height_in <- height_cm / 2.54
    
    # Cap dimensions to prevent PNG device errors
    max_dim <- 100  # Maximum dimension in inches
    if (width_in > max_dim) {
        cat(sprintf("Warning: Width (%.1f in) exceeds maximum. Capping at %d in.\n", width_in, max_dim))
        width_in <- max_dim
    }
    if (height_in > max_dim) {
        cat(sprintf("Warning: Height (%.1f in) exceeds maximum. Capping at %d in.\n", height_in, max_dim))
        height_in <- max_dim
    }

    cat(sprintf("Generating heatmap: %.1f x %.1f inches at %d dpi\n", width_in, height_in, dpi))

    # Save the heatmap
    png(output_file, width = width_in, height = height_in, units = "in", res = dpi)
    draw(ht,
        row_title = gene_label,
        row_title_side = "left",
        row_title_gp = gpar(fontsize = label_font, fontface = "bold"),
        column_title = column_label,
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize = label_font, fontface = "bold"),
        padding = unit(c(2.54, 2.54, 2.54, 2.54), "cm")
    )
    dev.off()
    
    cat(sprintf("Heatmap saved to: %s\n", output_file))
}

# Main function
main <- function() {
    option_list <- list(
        make_option(c("-i", "--input"),        type = "character", default = NULL,               help = "Input matrix file in TSV format."),
        make_option(c("-o", "--output"),       type = "character", default = NULL,               help = "Output heatmap file (e.g., heatmap.png)."),
        make_option("--gene_label",            type = "character", default = "Genes of Interest",help = "Label for the genes (rows)."),
        make_option("--column_label",          type = "character", default = "CAREs (CARE)",     help = "Label for the columns (CAREs or Functions)."),
        make_option("--row_font",              type = "integer",   default = 12L,                help = "Font size for row names."),
        make_option("--col_font",             type = "integer",   default = 10L,                help = "Font size for column names."),
        make_option("--label_font",            type = "integer",   default = 14L,                help = "Font size for labels."),
        make_option("--cell_size",             type = "double",    default = 1.0,                help = "Cell size in cm."),
        make_option("--cell_font_size",        type = "integer",   default = 10L,                help = "Font size for text inside cells."),
        make_option("--min_freq",              type = "integer",   default = 7L,                 help = "Minimum maximum frequency threshold for columns to be included."),
        make_option("--dpi",                   type = "integer",   default = 600L,               help = "Image DPI (default: 600)."),
        make_option("--color_palette",         type = "character", default = NULL,               help = "Comma-separated hex colors for the heatmap gradient."),
        make_option("--column_rotation",       type = "double",    default = 22.5,               help = "Column name rotation angle (default: 22.5)."),
        make_option("--cell_border_color",     type = "character", default = "white",            help = "Color of cell borders (default: white)."),
        make_option("--cell_border_width",     type = "double",    default = 2.0,                help = "Width of cell borders in pt (default: 2)."),
        make_option("--legend_height",         type = "double",    default = 3.0,                help = "Legend height in cm (default: 3).")
    )

    parser <- OptionParser(
        option_list  = option_list,
        description  = "Generate a heatmap from a PlantCARE matrix."
    )
    args <- parse_args(parser)

    if (is.null(args$input))  stop("--input is required.")
    if (is.null(args$output)) stop("--output is required.")

    create_heatmap(args$input, args$output, args$gene_label, args$column_label,
                   args$row_font, args$col_font, args$label_font, args$cell_size,
                   args$min_freq, args$cell_font_size, args$dpi,
                   args$color_palette, args$column_rotation, args$cell_border_color,
                   args$cell_border_width, args$legend_height)
}

if (!interactive()) {
    main()
}