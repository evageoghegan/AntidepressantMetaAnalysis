#Taking top differentially expressed genes and the datasets where they are leading edge in fgsea results
#nov 20th 2024
#old version from 2024

# Set working directory
setwd("/Users/evageoghegan/Brain Data Alchemy Project")

# Load required library
library(dplyr)

# Input file paths (update with your file paths)
top_genes_file <- "Top_Genes_FullMetaOutput.csv"       # CSV file containing a single column of top genes
fgsea_results_file <- "fgsea_significant_datasets_fullmetaoutput.csv"  # CSV file containing FGSEA results

# Output file path
filtered_results_file <- "filtered_fgsea_results.csv"

# Read the list of top genes
# Assumes a single column named "Gene" in the CSV
top_genes <- read.csv(top_genes_file)$Gene

# Read FGSEA results
# Assumes the FGSEA results CSV has a "pathway" column and a "leadingEdge" column
# The "leadingEdge" column should be a comma-separated string of genes (e.g., "GeneA,GeneB,GeneC")
fgsea_results <- read.csv(fgsea_results_file)

# Split the leadingEdge strings into lists of genes
fgsea_results$leadingEdge <- strsplit(as.character(fgsea_results$leadingEdge), ",")

# Filter to only include rows where leadingEdge contains any top gene
filtered_results <- fgsea_results %>%
  rowwise() %>%
  filter(any(top_genes %in% leadingEdge)) %>%
  ungroup()

# Deduplicate datasets (remove rows with identical pathways or leading edges)
filtered_results <- filtered_results %>%
  distinct(pathway, .keep_all = TRUE)

# Add a column showing the matching top genes for each dataset
filtered_results <- filtered_results %>%
  mutate(matching_genes = sapply(leadingEdge, function(genes) {
    paste(intersect(genes, top_genes), collapse = ",")
  }))

# Convert the leadingEdge column back to a comma-separated string
filtered_results$leadingEdge <- sapply(filtered_results$leadingEdge, paste, collapse = ",")

# Write the filtered results to a CSV file
write.csv(filtered_results, filtered_results_file, row.names = FALSE)

# Message to indicate completion
cat("Filtered FGSEA results saved to", filtered_results_file, "\n")

