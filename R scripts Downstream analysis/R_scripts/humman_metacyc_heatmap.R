# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# --- Load Data ---

# MetaCyc functional annotation file from HUMAnN3 output
metacyc_file <- "humann_combined_MetaCyc_named.tsv"
metadata_file <- "metadata.txt"

# Read data
metacyc_data <- fread(metacyc_file, header = TRUE, sep = "\t")
metadata <- fread(metadata_file)

# --- Clean and Reshape Data ---

# Remove unwanted suffixes from sample names
metacyc_data <- metacyc_data %>%
  rename_with(~ gsub("_1_kneaddata_paired_2_Abundance-RPKs", "", .), -`# Gene Family`)

# Pivot to long format
metacyc_long <- metacyc_data %>%
  pivot_longer(cols = -`# Gene Family`, names_to = "Sample", values_to = "Abundance")

# Clean metadata
metadata <- metadata %>%
  rename(Sample = `Sample Names`) %>%
  filter(Sample %in% unique(metacyc_long$Sample))

# Filter out unwanted gene families
metacyc_filtered <- metacyc_long %>%
  filter(!`# Gene Family` %in% c("UNMAPPED", "UNGROUPED", "NO_NAME"))

# Compute relative abundance
metacyc_abundance <- metacyc_filtered %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance, na.rm = TRUE) * 100) %>%
  ungroup()

# --- Prepare Heatmap Data ---

prepare_heatmap_data <- function(data_abundance, top_n = 30) {
  top_genes <- data_abundance %>%
    group_by(`# Gene Family`) %>%
    summarise(Total_Abundance = sum(Relative_Abundance, na.rm = TRUE)) %>%
    arrange(desc(Total_Abundance)) %>%
    slice_head(n = top_n) %>%
    pull(`# Gene Family`)
  
  data_filtered <- data_abundance %>%
    filter(`# Gene Family` %in% top_genes) %>%
    group_by(`# Gene Family`, Sample) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Sample, values_from = Relative_Abundance, values_fill = 0)

  rownames(data_filtered) <- make.unique(data_filtered$`# Gene Family`)
  heatmap_data <- data_filtered[, -1]
  heatmap_data <- as.data.frame(lapply(heatmap_data, as.numeric))
  rownames(heatmap_data) <- rownames(data_filtered)
  heatmap_data[is.na(heatmap_data)] <- 0

  log10(heatmap_data + 1e-6)
}

# --- Generate Heatmap ---

generate_heatmap <- function(heatmap_data, metadata, output_file, title) {
  metadata_clean <- metadata %>%
    select(Sample, `Sample Type`, `Sample Source`) %>%
    column_to_rownames("Sample")

  png(paste0(output_file, ".png"), width = 2200, height = 1000, res = 150)
  pheatmap(
    heatmap_data,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = metadata_clean,
    color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
    main = title,
    fontsize_row = 8,
    fontsize_col = 10,
    scale = "row"
  )
  dev.off()
}

# --- Run Analysis ---

metacyc_heatmap_data <- prepare_heatmap_data(metacyc_abundance, top_n = 30)
generate_heatmap(metacyc_heatmap_data, metadata, "MetaCyc_Heatmap_Corrected", "MetaCyc Pathway Abundance")

# Optionally print summarized table
print(head(metacyc_abundance))
