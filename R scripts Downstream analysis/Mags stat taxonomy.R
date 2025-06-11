# Load libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)

# 1. Load Data
gtdbk_data <- fread("gtdbtk.bac120.summary_latest.tsv", header = TRUE)
metadata <- fread("metadata.txt", header = TRUE)
checkm_data <- fread("checkm_summary.tsv", header = TRUE)

# 2. Merge Data
merged_data <- gtdbk_data %>%
  rename(Bin_Id = user_genome) %>%
  left_join(checkm_data, by = c("Bin_Id" = "Bin Id")) %>%
  mutate(Sample = sub("_.*", "", Bin_Id)) %>%
  left_join(metadata, by = c("Sample" = "Sample.Names")) %>%
  filter(!is.na(Completeness), !is.na(Contamination))

# 3. Define custom consistent color palette
sample_source_palette <- c(
  "Abijata" = "lightgreen",
  "Chitu" = "gold",
  "Shala" = "lightpink",
  "Cattle" = "tan",
  "Goat" = "salmon",
  "Sheep" = "mediumpurple"
)

# Make sure Sample Source follows correct order
merged_data$`Sample Source` <- factor(merged_data$`Sample Source`,
                                      levels = c("Abijata", "Chitu", "Shala", "Cattle", "Goat", "Sheep"))

# 4. Completeness vs Contamination Scatter Plot
plot_completeness_contamination <- ggplot(merged_data, aes(x = Contamination, y = Completeness, color = `Sample Source`)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "red") +
  scale_color_manual(values = sample_source_palette) +
  theme_bw() +
  labs(title = "Completeness vs Contamination",
       x = "Contamination (%)",
       y = "Completeness (%)",
       color = "Sample Source") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA)
  )

# 5. Genome Size Distribution Boxplot
plot_genome_size <- ggplot(merged_data, aes(x = `Sample Source`, y = `Genome size (bp)`, fill = `Sample Source`)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, shape = 21) +
  scale_fill_manual(values = sample_source_palette) +
  theme_bw() +
  labs(title = "Genome Size Distribution",
       y = "Genome Size (bp)", x = "Sample Source", fill = "Sample Source") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "none",
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA)
  )

# 6. RED Value Density Plot
merged_data$red_value <- as.numeric(merged_data$red_value)

plot_red_density <- ggplot(merged_data %>% filter(!is.na(red_value)), 
                           aes(x = red_value, fill = `Sample Source`)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = sample_source_palette) +
  theme_bw() +
  labs(title = "RED Value Density",
       x = "RED Value", y = "Density", fill = "Sample Source") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA)
  )

# 7. Arrange all three plots HORIZONTALLY
final_combined_plot_horizontal <- ggarrange(plot_completeness_contamination,
                                            plot_genome_size,
                                            plot_red_density,
                                            ncol = 3, nrow = 1,
                                            labels = c("A", "B", "C"),
                                            align = "hv",
                                            common.legend = TRUE, legend = "bottom")

# 8. Save the Combined Horizontal Plot
ggsave("MAGs_Quality_GenomeSize_RED_Final_Horizontal_Clean.png",
       plot = final_combined_plot_horizontal,
       width = 18, height = 6, dpi = 300)


# Load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# 1. Load Data
gtdbtk_data <- fread("formatted_gtdbtk_data_latest.tsv", header = TRUE, sep = "\t")
metadata <- fread("metadata.txt", header = TRUE)

# 2. Clean missing taxa
gtdbtk_data$Phylum[is.na(gtdbtk_data$Phylum) | gtdbtk_data$Phylum == ""] <- "Unclassified"
gtdbtk_data$Genus[is.na(gtdbtk_data$Genus) | gtdbtk_data$Genus == ""] <- "Unclassified"

# 3. Extract Sample Names
gtdbtk_data$Sample <- sub("_.*", "", gtdbtk_data$user_genome)

# Fix all names
metadata <- metadata %>%
  rename(Sample.Names = `Sample Names`,
         Sample.Type = `Sample Type`,
         Sample.Source = `Sample Source`)

# 4. Merge metadata
merged_data <- gtdbtk_data %>%
  left_join(metadata, by = c("Sample" = "Sample.Names"))

# 5. Prepare Phylum and Genus count tables
phylum_table <- merged_data %>%
  group_by(Phylum, Sample.Source) %>%
  summarise(MAG_Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Sample.Source, values_from = MAG_Count, values_fill = 0)

genus_table <- merged_data %>%
  group_by(Genus, Sample.Source) %>%
  summarise(MAG_Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Sample.Source, values_from = MAG_Count, values_fill = 0)

# 6. Keep Top 20 Phylum and Genus
top_phyla <- phylum_table %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(Total)) %>%
  slice(1:20) %>%
  select(-Total)

top_genera <- genus_table %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(Total)) %>%
  slice(1:20) %>%
  select(-Total)

# 7. Create matrices
phylum_matrix <- as.data.frame(top_phyla)
rownames(phylum_matrix) <- phylum_matrix$Phylum
phylum_matrix <- phylum_matrix[,-1]

genus_matrix <- as.data.frame(top_genera)
rownames(genus_matrix) <- genus_matrix$Genus
genus_matrix <- genus_matrix[,-1]

# 8. Align samples
common_samples <- intersect(colnames(phylum_matrix), colnames(genus_matrix))
phylum_matrix <- phylum_matrix[, common_samples]
genus_matrix <- genus_matrix[, common_samples]

# 9. Prepare metadata for annotation
sample_meta <- metadata %>%
  filter(Sample.Source %in% common_samples) %>%
  distinct(Sample.Source, Sample.Type)

sample_meta <- as.data.frame(sample_meta)
rownames(sample_meta) <- sample_meta$Sample.Source








# --- Fix sample order manually ---
sample_meta$Sample.Type <- factor(sample_meta$Sample.Type, levels = c("Soda Lake", "Rumen"))
sample_meta <- sample_meta %>% arrange(Sample.Type, Sample.Source)
ordered_samples <- rownames(sample_meta)

phylum_matrix <- phylum_matrix[, ordered_samples]
genus_matrix <- genus_matrix[, ordered_samples]

# --- Define Colors ---
library(RColorBrewer)

sample_type_colors <- c("Soda Lake" = "skyblue", "Rumen" = "orange")
sample_source_colors <- setNames(
  viridis(length(unique(sample_meta$Sample.Source)), option = "turbo"),
  unique(sample_meta$Sample.Source)
)

# --- Metadata Annotation ---
column_annot <- HeatmapAnnotation(
  SampleType = sample_meta$Sample.Type,
  SampleSource = sample_meta$Sample.Source,
  col = list(
    SampleType = sample_type_colors,
    SampleSource = sample_source_colors
  ),
  show_annotation_name = FALSE,
  annotation_height = unit(c(4, 4), "mm")
)

# --- Color Palettes for Heatmaps (YOUR request)
phylum_palette <- brewer.pal(9, "YlGnBu")
genus_palette  <- brewer.pal(9, "RdPu")

# --- Convert to pure matrices ---
phylum_matrix <- as.matrix(phylum_matrix)
genus_matrix  <- as.matrix(genus_matrix)

# --- Now safely calculate quantiles
phylum_upper_limit <- quantile(as.vector(phylum_matrix), 0.9)
genus_upper_limit  <- quantile(as.vector(genus_matrix), 0.9)

# --- Then continue normal coloring
col_fun_phylum <- colorRamp2(c(0, phylum_upper_limit), c(phylum_palette[1], phylum_palette[9]))
col_fun_genus  <- colorRamp2(c(0, genus_upper_limit), c(genus_palette[1], genus_palette[9]))




# --- Create Heatmaps with Borders ---
ht_phylum <- Heatmap(
  as.matrix(phylum_matrix),
  name = "Phylum MAG Count",
  col = col_fun_phylum,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  border = TRUE,  # Black cell borders
  top_annotation = column_annot,
  column_split = sample_meta$Sample.Type,
  show_column_names = TRUE,
  show_row_names = TRUE,
  show_column_dend = FALSE,
  row_title = "Top 20 Phyla",
  heatmap_legend_param = list(title_position = "topcenter")
)

ht_genus <- Heatmap(
  as.matrix(genus_matrix),
  name = "Genus MAG Count",
  col = col_fun_genus,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  border = TRUE,  # Black cell borders
  top_annotation = NULL,
  column_split = sample_meta$Sample.Type,
  show_column_names = TRUE,
  show_row_names = TRUE,
  show_column_dend = FALSE,
  row_title = "Top 20 Genera",
  heatmap_legend_param = list(title_position = "topcenter")
)

# --- Save as High Quality PNG ---
png("MAGs_Top20_Phylum_Genus_Heatmap_Final_Bordered_BrewerColors_v2.png", width = 5400, height = 3600, res = 300)
draw(ht_phylum %v% ht_genus,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     gap = unit(4, "mm"))
dev.off()
