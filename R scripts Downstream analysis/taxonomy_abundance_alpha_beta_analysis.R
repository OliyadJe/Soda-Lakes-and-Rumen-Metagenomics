# Load libraries
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(stringr)
library(grid)

# Load metaphlan data and metadata
metaphlan_file <- "merged_abundance_table_latest.txt"
metadata_file <- "metadata.txt"
metaphlan_data <- fread(metaphlan_file, header = TRUE)
metadata <- fread(metadata_file, header = TRUE)

# Clean sample names
clean_sample_names <- function(names_vector) {
  str_replace(names_vector, "_1_metaphlan$", "")
}
colnames(metaphlan_data)[2:ncol(metaphlan_data)] <- clean_sample_names(colnames(metaphlan_data)[2:ncol(metaphlan_data)])

# Split clade_name
taxonomy_split <- as.data.table(tstrsplit(metaphlan_data$clade_name, "\\|", fill = NA))
colnames(taxonomy_split) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
metaphlan_data <- cbind(taxonomy_split, metaphlan_data[, -1])

# Summarize Phylum and Genus
phylum_abundance <- metaphlan_data %>%
  filter(!is.na(Phylum), is.na(Genus)) %>%
  mutate(Phylum = gsub("p__", "", Phylum)) %>%
  group_by(Phylum) %>%
  summarize(across(where(is.numeric), sum))

genus_abundance <- metaphlan_data %>%
  filter(!is.na(Genus)) %>%
  mutate(Genus = gsub("g__", "", Genus)) %>%
  group_by(Genus) %>%
  summarize(across(where(is.numeric), sum))

# Top 10 Phylum and Genus
phylum_abundance$Total_Abundance <- rowSums(phylum_abundance[, 2:(ncol(phylum_abundance))])
top_phylum <- phylum_abundance %>%
  arrange(desc(Total_Abundance)) %>%
  slice(1:10)
phylum_matrix <- as.data.frame(top_phylum)
rownames(phylum_matrix) <- phylum_matrix$Phylum
phylum_matrix <- phylum_matrix[, -c(1, ncol(phylum_matrix))]

genus_abundance$Total_Abundance <- rowSums(genus_abundance[, 2:(ncol(genus_abundance))])
top_genus <- genus_abundance %>%
  arrange(desc(Total_Abundance)) %>%
  slice(1:10)
genus_matrix <- as.data.frame(top_genus)
rownames(genus_matrix) <- genus_matrix$Genus
genus_matrix <- genus_matrix[, -c(1, ncol(genus_matrix))]

# Transform to log10
phylum_matrix_log <- log10(phylum_matrix + 1e-5)
genus_matrix_log <- log10(genus_matrix + 1e-5)

# Arrange metadata and samples
colnames(metadata) <- make.names(colnames(metadata))
metadata <- metadata %>%
  mutate(Sample.Names = as.character(Sample.Names)) %>%
  arrange(Sample.Type, Sample.Source)
ordered_samples <- metadata$Sample.Names

phylum_matrix_log <- phylum_matrix_log[, ordered_samples]
genus_matrix_log <- genus_matrix_log[, ordered_samples]

# Define different color palettes
phylum_palette <- brewer.pal(9, "YlGnBu")  # Yellow-Green-Blue for Phylum
genus_palette <- brewer.pal(9, "RdPu")     # Red-Purple for Genus

col_fun_phylum <- colorRamp2(
  seq(min(phylum_matrix_log), max(phylum_matrix_log), length.out = 9),
  phylum_palette
)
col_fun_genus <- colorRamp2(
  seq(min(genus_matrix_log), max(genus_matrix_log), length.out = 9),
  genus_palette
)

# Metadata annotation
column_annot <- HeatmapAnnotation(
  SampleType = metadata$Sample.Type,
  SampleSource = metadata$Sample.Source,
  pH = metadata$pH,
  Salinity = metadata$Salinity....,
  col = list(
    SampleType = c("Soda Lake" = "skyblue", "Rumen" = "orange"),
    SampleSource = c("Abijata" = "lightgreen", "Chitu" = "gold", "Shala" = "lightpink",
                     "Cattle" = "tan", "Goat" = "salmon", "Sheep" = "mediumpurple")
  ),
  annotation_legend_param = list(
    SampleType = list(title = "Sample Type"),
    SampleSource = list(title = "Sample Source"),
    pH = list(title = "pH"),
    Salinity = list(title = "Salinity (%)")
  )
)

# Create Phylum Heatmap
ht_phylum <- Heatmap(as.matrix(phylum_matrix_log),
                     name = "Phylum log10 Abundance",
                     col = col_fun_phylum,
                     border = TRUE,                # Thin white grid
                     cluster_rows = TRUE,
                     cluster_columns = TRUE,
                     show_column_dend = FALSE,      # Remove dotted dendrogram line
                     column_order = ordered_samples,
                     column_title = "",
                     row_title = "Top 10 Phylum ",
                     show_row_names = TRUE,
                     row_names_gp = gpar(fontsize = 8),
                     column_names_gp = gpar(fontsize = 8),
                     column_names_rot = 45,
                     top_annotation = column_annot,
                     column_split = metadata$Sample.Type,
                     column_gap = unit(2, "mm"))

# Create Genus Heatmap
ht_genus <- Heatmap(as.matrix(genus_matrix_log),
                    name = "Genus log10 Abundance",
                    col = col_fun_genus,
                    border = TRUE,                 # Thin white grid
                    cluster_rows = TRUE,
                    cluster_columns = TRUE,
                    show_column_dend = FALSE,       # Remove dotted dendrogram line
                    column_order = ordered_samples,
                    column_title = "",
                    row_title = "Top 10 Genus ",
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 8),
                    column_names_rot = 45,
                    top_annotation = NULL,
                    column_split = metadata$Sample.Type,
                    column_gap = unit(2, "mm"))

# Draw final plot
png("final_phylum_genus_polished3.png", width = 5200, height = 3000, res = 300)
draw(ht_phylum %v% ht_genus,
     merge_legends = TRUE,
     gap = unit(5, "mm"),         # Extra gap between Phylum and Genus
     heatmap_legend_side = "right")
dev.off()


# Load libraries
library(data.table)
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(tibble)
library(patchwork)
library(ggpubr)

# 1. Load your data
metaphlan_file <- "merged_abundance_table_latest.txt"
metadata_file <- "metadata.txt"
metaphlan_data <- fread(metaphlan_file, header = TRUE)
metadata <- fread(metadata_file, header = TRUE)

# 2. Clean sample names
clean_sample_names <- function(names_vector) {
  stringr::str_replace(names_vector, "_1_metaphlan$", "")
}
colnames(metaphlan_data)[2:ncol(metaphlan_data)] <- clean_sample_names(colnames(metaphlan_data)[2:ncol(metaphlan_data)])

# 3. Split taxonomy
taxonomy_split <- as.data.table(tstrsplit(metaphlan_data$clade_name, "\\|", fill = NA))
colnames(taxonomy_split) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
metaphlan_data <- cbind(taxonomy_split, metaphlan_data[, -1])

# 4. Extract species-level matrix
species_abundance <- metaphlan_data %>%
  filter(!is.na(Species)) %>%
  mutate(Species = gsub("s__", "", Species)) %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), sum))

# Set Species as rownames
rownames(species_abundance) <- species_abundance$Species
species_abundance <- species_abundance[, !(colnames(species_abundance) %in% "Species")]

# 5. Transpose and normalize
species_matrix_t <- as.data.frame(t(species_abundance))
species_matrix_t <- species_matrix_t / rowSums(species_matrix_t)

# 6. Prepare metadata
colnames(metadata) <- make.names(colnames(metadata))
metadata <- metadata %>%
  mutate(Sample.Names = as.character(Sample.Names)) %>%
  filter(Sample.Names %in% rownames(species_matrix_t))

# 7. Calculate Shannon index
alpha_shannon <- data.frame(
  Sample = rownames(species_matrix_t),
  Shannon = diversity(species_matrix_t, index = "shannon")
)

# Merge metadata
alpha_shannon <- left_join(alpha_shannon, metadata, by = c("Sample" = "Sample.Names"))

# 8. Order Sample Source
alpha_shannon$Sample.Source <- factor(alpha_shannon$Sample.Source,
                                      levels = c("Abijata", "Chitu", "Shala", "Cattle", "Goat", "Sheep"))

# 9. Set color palette
sample_source_palette <- c(
  "Abijata" = "lightgreen",
  "Chitu" = "gold",
  "Shala" = "lightpink",
  "Cattle" = "tan",
  "Goat" = "salmon",
  "Sheep" = "mediumpurple"
)

sample_type_palette <- c(
  "Soda Lake" = "skyblue",
  "Rumen" = "orange"
)

# 10. Create Alpha Diversity plot
p_shannon_combined <- ggplot(alpha_shannon, aes(x = Sample.Source, y = Shannon, fill = Sample.Source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, shape = 21) +
  scale_fill_manual(values = sample_source_palette) +
  facet_wrap(~ Sample.Type, scales = "free_x", nrow = 1) +
  theme_bw() +
  labs(title = "Alpha Diversity (Shannon Index)", 
       y = "Shannon Index", x = "Sample Source") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))

# 11. Beta Diversity (PCoA)

# Bray-Curtis distance
bray_dist <- vegdist(species_matrix_t, method = "bray")

# PCoA
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)

# Variance explained
var_explained <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 1)

# Prepare PCoA data
pcoa_df <- data.frame(
  Sample = rownames(pcoa_result$points),
  Axis1 = pcoa_result$points[,1],
  Axis2 = pcoa_result$points[,2]
)

pcoa_df <- left_join(pcoa_df, metadata, by = c("Sample" = "Sample.Names"))

# PCoA plot
p_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Sample.Type, shape = Sample.Source)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = sample_type_palette) +
  theme_bw() +
  labs(title = "Beta Diversity (PCoA Bray-Curtis)",
       x = paste0("PCoA1 (", var_explained[1], "%)"),
       y = paste0("PCoA2 (", var_explained[2], "%)"),
       color = "Sample Type", shape = "Sample Source") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14))

# 12. PERMANOVA
permanova_result <- adonis2(bray_dist ~ Sample.Type, data = metadata, permutations = 999, method = "bray")
print(permanova_result)

# OPTIONAL: Add PERMANOVA p-value to PCoA title
pcoa_pval <- signif(permanova_result$`Pr(>F)`[1], 2)
p_pcoa <- p_pcoa + labs(title = paste0("Beta Diversity (PERMANOVA p = ", pcoa_pval, ")"))

# 13. Combine Alpha and Beta plots
combined_plot <- p_shannon_combined + p_pcoa + plot_layout(ncol = 2)

# Save combined plot
ggsave("alpha_beta_combined2.png", plot = combined_plot, width = 14, height = 6, dpi = 300)

# Swap colors and shapes
p_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Sample.Source, shape = Sample.Type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = sample_source_palette) +
  scale_shape_manual(values = c("Soda Lake" = 16, "Rumen" = 17)) +
  theme_bw() +
  labs(title = paste0("Beta Diversity (PERMANOVA p = ", signif(permanova_result$`Pr(>F)`[1], 2),
                      ", R² = ", round(permanova_result$R2[1], 2), ")"),
       x = paste0("PCoA1 (", var_explained[1], "%)"),
       y = paste0("PCoA2 (", var_explained[2], "%)"),
       color = "Sample Source", shape = "Sample Type") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14))

# Combine with more space between plots
combined_plot <- p_shannon_combined + p_pcoa + 
  plot_layout(ncol = 2, widths = c(1, 1.2)) & 
  theme(legend.position = "bottom")

# Save final plot
ggsave("alpha_beta_combined_final2.png", plot = combined_plot, width = 16, height = 7, dpi = 300)
