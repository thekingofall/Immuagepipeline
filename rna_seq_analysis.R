#!/usr/bin/env Rscript

#' RNA-seq Analysis Script
#' 
#' Usage: Rscript rna_seq_analysis.R [options]
#' 
#' Options:
#'   -i, --input     Input count file path
#'   -n, --name      Base name for output directory (default: "New")
#'   -o, --out       Output directory (if not specified, will use name_A05.Rcount)
#'   -m, --min       Minimum count threshold (default: 0)
#'   -w, --width     Plot width in inches (default: 10)
#'   -h, --height    Plot height in inches (default: 6)
#'   -d, --dpi       Plot resolution (default: 300)

# Load required packages
required_packages <- c(
  "data.table", "dplyr", "ggplot2", "pheatmap",
  "clusterProfiler", "org.Hs.eg.db", "FactoMineR", "factoextra",
  "optparse"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", 
              help="Input count file path"),
  make_option(c("-n", "--name"), type="character", default="New",
              help="Base name for output directory [default=%default]"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Output directory (if not specified, will use name_A05.Rcount)"),
  make_option(c("-m", "--min"), type="numeric", default=0,
              help="Minimum count threshold [default=%default]"),
  make_option(c("-w", "--width"), type="numeric", default=10,
              help="Plot width in inches [default=%default]"),
  make_option(c("--height"), type="numeric", default=6,
              help="Plot height in inches [default=%default]"),
  make_option(c("-d", "--dpi"), type="numeric", default=300,
              help="Plot resolution [default=%default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Set output directory based on name or out parameter
if (is.null(opt$out)) {
  opt$output <- paste0(opt$name, "_A05.Rcount")
} else {
  opt$output <- paste0(opt$out, "_A05.Rcount")
}

# Default settings
plot_colors <- list(
  expression = c(expressed = "#377EB8", non_expressed = "#E41A1C"),
  pca = c("#ff7c38", "#e03e36", "#17b978", "#b80d57",
          "#700961", "#2f9296", "#46b7b9", "#87dfd6")
)

metadata_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "SYMBOL")

#' Process RNA-seq count data
process_rnaseq_data <- function(count_file) {
  data <- fread(count_file) %>% as.data.frame()
  rownames(data) <- data[,1]
  data <- data[,-1]
  
  colnames(data) <- gsub(".*A03.Alignment/", "", colnames(data)) %>%
                    gsub(".sorted.bam", "", .)
  
  return(data)
}

#' Convert counts to TPM
countToTpm <- function(counts, effLen) {
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

#' Convert gene IDs
convert_gene_ids <- function(data) {
  # First make a copy of the input data
  data_copy <- data
  
  # Add Geneid column if it doesn't exist
  if (!"Geneid" %in% colnames(data_copy)) {
    data_copy$Geneid <- rownames(data_copy)
  }
  
  # Clean ENSEMBL IDs - handle different possible formats
  data_copy$Geneid <- gsub("\\..*$", "", data_copy$Geneid)  # Remove version numbers if present
  
  # Print some diagnostic information
  message("Number of input genes: ", nrow(data_copy))
  message("Sample of gene IDs: ", paste(head(data_copy$Geneid), collapse=", "))
  
  # Perform the conversion
  result <- tryCatch({
    bitr(data_copy$Geneid,
         fromType = "ENSEMBL",
         toType = "SYMBOL",
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("Error in gene ID conversion: ", e$message)
    return(NULL)
  })
  
  # Check if conversion was successful
  if (is.null(result) || nrow(result) == 0) {
    message("No genes could be converted. Returning original data with ENSEMBL IDs as symbols")
    data_copy$SYMBOL <- data_copy$Geneid
    return(data_copy)
  }
  
  # Join the conversion results with original data
  data_with_symbols <- data_copy %>%
    left_join(result, by = c("Geneid" = "ENSEMBL")) %>%
    filter(!is.na(SYMBOL)) %>%
    filter(!duplicated(SYMBOL))
  
  # Check if we have any results
  if (nrow(data_with_symbols) == 0) {
    message("No matches found after joining. Returning original data with ENSEMBL IDs as symbols")
    data_copy$SYMBOL <- data_copy$Geneid
    return(data_copy)
  }
  
  message("Successfully converted ", nrow(data_with_symbols), " genes")
  
  rownames(data_with_symbols) <- data_with_symbols$SYMBOL
  return(data_with_symbols)
}

#' Create expression plots
create_expression_plots <- function(data, output_dir) {
  stats_df <- data.frame(
    Sample = colnames(data),
    Expressed = colSums(data > opt$min),
    Non_expressed = colSums(data == opt$min)
  )
  
  # Save gene count statistics to CSV
  write.csv(stats_df, 
            file.path(output_dir, "gene_count_statistics.csv"), 
            row.names = FALSE)
  
  p1 <- ggplot(stats_df, aes(x = Sample, y = Expressed)) +
    geom_bar(stat = "identity", fill = plot_colors$expression["expressed"]) +
    geom_text(aes(label = Expressed), 
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4,
              fontface = "bold") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sample", 
         y = "Number of Expressed Genes", 
         title = "Number of Expressed Genes per Sample")
  
  ggsave(file.path(output_dir, "expressed_genes.pdf"), 
         p1, 
         width = opt$width, 
         height = opt$height,
         dpi = opt$dpi)
  
  return(p1)
}

#' Perform PCA analysis
perform_pca_analysis <- function(data, output_dir) {
  pca.res <- prcomp(log10(t(data)+1))
  eigen_values <- get_eigenvalue(pca.res)
  
  p <- fviz_pca_ind(pca.res,
                    habillage = as.factor(colnames(data)),
                    palette = plot_colors$pca,
                    addEllipses = TRUE,
                    label = FALSE) +
    theme_classic() +
    ggtitle("PCA Analysis")
  
  ggsave(file.path(output_dir, "pca_plot.pdf"), 
         p, 
         width = opt$width, 
         height = opt$height,
         dpi = opt$dpi)
  
  write.csv(eigen_values, file.path(output_dir, "pca_eigenvalues.csv"))
  
  return(list(pca_result = pca.res, pca_plot = p))
}

#' Main analysis function
main <- function() {
  # Check if input file is provided
  if (is.null(opt$input)) {
    stop("Please provide input count file path using -i or --input option")
  }
  
  # Create output directory
  dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
  
  # Process data
  message("Processing count data...")
  data <- process_rnaseq_data(opt$input)
  message("Raw data dimensions: ", nrow(data), " x ", ncol(data))
  message("Number of genes with any expression: ", 
          sum(rowSums(data > 0) > 0))
  
  # Convert gene IDs
  message("Converting gene IDs...")
  data_with_symbols <- convert_gene_ids(data)
  message("After ID conversion dimensions: ", 
          nrow(data_with_symbols), " x ", ncol(data_with_symbols))
  
  # Filter for expressed genes
  message("Filtering expressed genes...")
  data_cols <- setdiff(colnames(data_with_symbols), metadata_cols)
  expressed_data <- data_with_symbols[
    rowSums(data_with_symbols[, data_cols] > opt$min) > 0, 
    data_cols
  ]
  message("Final expressed genes: ", nrow(expressed_data))
  
  # Calculate TPM
  message("Calculating TPM values...")
  tpm_data <- apply(expressed_data, 2, 
                    function(x) countToTpm(x, data_with_symbols$Length)) %>% 
    as.data.frame()
  
  # Create plots
  message("Creating expression plots...")
  expr_plots <- create_expression_plots(expressed_data, opt$output)
  
  message("Performing PCA analysis...")
  pca_results <- perform_pca_analysis(tpm_data, opt$output)
  
  # Create heatmap
  message("Creating heatmap...")
  pdf(file.path(opt$output, "heatmap.pdf"),
      width = opt$width,
      height = opt$height)
  pheatmap(tpm_data, scale = "row", 
           show_rownames = FALSE, 
           cluster_cols = FALSE,
           main = "Gene Expression Heatmap")
  dev.off()
  
  # Save results
  message("Saving results...")
  write.csv(tpm_data, file.path(opt$output, "tpm_data.csv"))
  write.csv(expressed_data, file.path(opt$output, "raw_counts_filtered.csv"))
  
  # Save session info
  writeLines(capture.output(sessionInfo()), 
            file.path(opt$output, "session_info.txt"))
  
  message("Analysis complete! Results saved in: ", opt$output)
}

# Run the analysis
if (!interactive()) {
  main()
}
