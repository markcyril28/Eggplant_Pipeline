#!/usr/bin/env Rscript
#
# Visualization script for PlantCARE matrix results
# Creates heatmaps and bar plots for cis-regulatory element analysis
#

# Required libraries
required_packages <- c("ggplot2", "pheatmap", "RColorBrewer", "tidyr", "dplyr")

# Function to install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      install.packages(pkg, repos = "http://cran.us.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
}

# Install and load packages
install_if_missing(required_packages)

# Function to read and visualize count matrix
visualize_count_matrix <- function(count_file, output_prefix) {
  cat("Reading count matrix:", count_file, "\n")
  
  # Read data
  count_data <- read.table(count_file, header = TRUE, row.names = 1, 
                           sep = "\t", check.names = FALSE)
  
  # Remove columns with all zeros
  count_data <- count_data[, colSums(count_data) > 0]
  
  cat("Matrix dimensions:", nrow(count_data), "sequences x", 
      ncol(count_data), "motif types\n")
  
  # 1. Heatmap of motif counts
  if (ncol(count_data) > 1) {
    cat("Creating heatmap...\n")
    
    pdf(paste0(output_prefix, "_heatmap.pdf"), width = 12, height = 8)
    
    # Transpose for better visualization if many motifs
    if (ncol(count_data) > nrow(count_data)) {
      plot_data <- t(as.matrix(count_data))
      main_title <- "Cis-Regulatory Element Counts (Motifs as rows)"
    } else {
      plot_data <- as.matrix(count_data)
      main_title <- "Cis-Regulatory Element Counts"
    }
    
    pheatmap(plot_data,
             color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             display_numbers = TRUE,
             number_format = "%.0f",
             fontsize_number = 8,
             fontsize_row = 8,
             fontsize_col = 10,
             main = main_title,
             cellwidth = 15,
             cellheight = 15,
             border_color = "grey60")
    
    dev.off()
    cat("Saved:", paste0(output_prefix, "_heatmap.pdf\n"))
  }
  
  # 2. Bar plot of total motifs per type
  cat("Creating bar plot...\n")
  
  motif_totals <- data.frame(
    Motif = colnames(count_data),
    Count = colSums(count_data),
    stringsAsFactors = FALSE
  )
  motif_totals <- motif_totals[order(motif_totals$Count, decreasing = TRUE), ]
  motif_totals$Motif <- factor(motif_totals$Motif, levels = motif_totals$Motif)
  
  p <- ggplot(motif_totals, aes(x = Motif, y = Count, fill = Count)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold")) +
    labs(title = "Total Count of Cis-Regulatory Elements",
         x = "Motif Type",
         y = "Count") +
    geom_text(aes(label = Count), vjust = -0.5, size = 3)
  
  ggsave(paste0(output_prefix, "_barplot.pdf"), p, 
         width = 12, height = 6)
  cat("Saved:", paste0(output_prefix, "_barplot.pdf\n"))
  
  return(count_data)
}

# Function to visualize functional categories
visualize_functional_categories <- function(category_file, output_prefix) {
  cat("\nReading functional category matrix:", category_file, "\n")
  
  # Read data
  category_data <- read.table(category_file, header = TRUE, row.names = 1,
                              sep = "\t", check.names = FALSE)
  
  cat("Categories found:", ncol(category_data), "\n")
  
  # Prepare data for plotting
  plot_data <- data.frame(
    Category = colnames(category_data),
    Count = colSums(category_data),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data[order(plot_data$Count, decreasing = TRUE), ]
  plot_data$Category <- factor(plot_data$Category, levels = plot_data$Category)
  
  # Color palette for categories
  category_colors <- c(
    "Light Responsiveness" = "#FFD700",
    "Drought Response" = "#8B4513",
    "ABA Response" = "#4169E1",
    "Anaerobic Response" = "#32CD32",
    "Stress Response" = "#FF4500",
    "Core Promoter" = "#9370DB",
    "Transcription Factor Binding" = "#FF69B4",
    "Other/Unknown" = "#808080"
  )
  
  # Pie chart
  p1 <- ggplot(plot_data, aes(x = "", y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = category_colors) +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right") +
    labs(title = "Distribution of Functional Categories") +
    geom_text(aes(label = paste0(Category, "\n(", Count, ")")),
              position = position_stack(vjust = 0.5), size = 3)
  
  ggsave(paste0(output_prefix, "_functional_pie.pdf"), p1,
         width = 10, height = 6)
  cat("Saved:", paste0(output_prefix, "_functional_pie.pdf\n"))
  
  # Bar plot
  p2 <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = category_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "none") +
    labs(title = "Functional Categories of Cis-Regulatory Elements",
         x = "Functional Category",
         y = "Count") +
    geom_text(aes(label = Count), vjust = -0.5, size = 4)
  
  ggsave(paste0(output_prefix, "_functional_bar.pdf"), p2,
         width = 10, height = 6)
  cat("Saved:", paste0(output_prefix, "_functional_bar.pdf\n"))
  
  return(category_data)
}

# Function to create summary report
create_summary_report <- function(summary_file, output_prefix) {
  cat("\nReading summary statistics:", summary_file, "\n")
  
  summary_data <- read.table(summary_file, header = TRUE, sep = "\t")
  
  # Print summary
  cat("\n")
  cat("=" %R% 60, "\n")
  cat("SUMMARY STATISTICS\n")
  cat("=" %R% 60, "\n")
  print(summary_data)
  cat("=" %R% 60, "\n")
  
  return(summary_data)
}

# Main function
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript visualize_plantCARE_matrix.R <prefix> [output_prefix]\n")
    cat("\nExample:\n")
    cat("  Rscript visualize_plantCARE_matrix.R PlantCARE_5913\n")
    cat("\nThis will look for files:\n")
    cat("  - PlantCARE_5913_count_matrix.tsv\n")
    cat("  - PlantCARE_5913_functional_categories.tsv\n")
    cat("  - PlantCARE_5913_summary.tsv\n")
    quit(status = 1)
  }
  
  prefix <- args[1]
  output_prefix <- ifelse(length(args) >= 2, args[2], paste0(prefix, "_plots"))
  
  cat("\n")
  cat("=" %R% 60, "\n")
  cat("PlantCARE Matrix Visualization\n")
  cat("=" %R% 60, "\n")
  cat("Input prefix:", prefix, "\n")
  cat("Output prefix:", output_prefix, "\n")
  
  # Define input files
  count_file <- paste0(prefix, "_count_matrix.tsv")
  category_file <- paste0(prefix, "_functional_categories.tsv")
  summary_file <- paste0(prefix, "_summary.tsv")
  
  # Visualize count matrix
  if (file.exists(count_file)) {
    count_data <- visualize_count_matrix(count_file, output_prefix)
  } else {
    cat("Warning: Count matrix file not found:", count_file, "\n")
  }
  
  # Visualize functional categories
  if (file.exists(category_file)) {
    category_data <- visualize_functional_categories(category_file, output_prefix)
  } else {
    cat("Warning: Functional category file not found:", category_file, "\n")
  }
  
  # Display summary
  if (file.exists(summary_file)) {
    summary_data <- create_summary_report(summary_file, output_prefix)
  } else {
    cat("Warning: Summary file not found:", summary_file, "\n")
  }
  
  cat("\n")
  cat("=" %R% 60, "\n")
  cat("Visualization complete!\n")
  cat("=" %R% 60, "\n")
}

# Helper function for string repetition
"%R%" <- function(x, n) {
  paste(rep(x, n), collapse = "")
}

# Run main function
if (!interactive()) {
  main()
}
