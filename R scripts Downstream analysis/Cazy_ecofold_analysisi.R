# --- Load Required Libraries ---
library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

# -------------------------------
# 1. Load and Merge Metadata
# -------------------------------
metadata <- fread("metadata.txt")
annotations <- fread("merged_ecofolddb_annotations.tsv")

# Clean column names and merge
colnames(metadata) <- c("Sample.Names", "Sample.Type", "Sample.Source", "pH", "Salinity")
annotations[, Sample := tstrsplit(MAG_ID, "_", fixed = TRUE)[1]]
merged_data <- merge(annotations, metadata, by.x = "Sample", by.y = "Sample.Names", all.x = TRUE)

# -------------------------------
# 2. Define Color Palettes
# -------------------------------
sample_source_palette <- c(
  "Abijata" = "lightgreen", "Chitu" = "gold", "Shala" = "lightpink",
  "Cattle" = "tan", "Goat" = "salmon", "Sheep" = "mediumpurple"
)
sample_type_palette <- c("Soda Lake" = "skyblue", "Rumen" = "orange")

# -------------------------------
# 3. Panel A: Functional Category Composition
# -------------------------------
category_df <- merged_data[, .N, by = .(`Sample Type`, `Sample Source`, Category)]
category_df[, Total := sum(N), by = .(`Sample Type`, `Sample Source`)]
category_df[, RelAbundance := 100 * N / Total]
category_df[, `Sample Type` := factor(`Sample Type`, levels = c("Soda Lake", "Rumen"))]

num_categories <- length(unique(category_df$Category))
category_palette <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(num_categories),
  unique(category_df$Category)
)

plot_A <- ggplot(category_df, aes(x = `Sample Source`, y = RelAbundance, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  facet_wrap(~ `Sample Type`, scales = "free_x") +
  scale_fill_manual(values = category_palette) +
  labs(x = "Sample Source", y = "Relative Abundance (%)", fill = "Functional Category") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))

# -------------------------------
# 4. Panel B: Carbon Pathways Dot Plot
# -------------------------------
carbon_df <- merged_data[Category == "Carbon cycling"]
pathway_df <- carbon_df[, .N, by = .(`Sample Type`, `Sample Source`, `Pathway/Activity`)]
pathway_df[, Total := sum(N), by = .(`Sample Type`, `Sample Source`)]
pathway_df[, RelAbundance := 100 * N / Total]

top_pathways <- pathway_df[, .(TotalCount = sum(N)), by = `Pathway/Activity`][order(-TotalCount)][1:10, `Pathway/Activity`]
pathway_df_top <- pathway_df[`Pathway/Activity` %in% top_pathways]

plot_B_dot <- ggplot(pathway_df_top, aes(x = RelAbundance, y = reorder(`Pathway/Activity`, RelAbundance))) +
  geom_point(aes(color = `Sample Type`, shape = `Sample Source`), size = 4, alpha = 0.85) +
  scale_color_manual(values = sample_type_palette) +
  scale_shape_manual(values = c("Abijata" = 16, "Chitu" = 17, "Shala" = 15, "Cattle" = 3, "Goat" = 7, "Sheep" = 8)) +
  labs(x = "Relative Abundance (%)", y = "Carbon Cycle Pathway") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 12, face = "bold"))

# -------------------------------
# 5. Panel C: CAZy GH Heatmap
# -------------------------------
cazy_metadata <- fread("merged_cazy_results_with_samples.tsv")
cazy_metadata <- cazy_metadata %>%
  mutate(Parent_Sample = sub("_.*", "", Sample_ID)) %>%
  left_join(metadata, by = c("Parent_Sample" = "Sample.Names")) %>%
  separate(CAZy_Class, into = c("Protein_ID", "CAZy_Info", "Extra_Info"), sep = "\\|", fill = "right", remove = FALSE) %>%
  mutate(CAZy_Class = sub("([A-Za-z]+).*", "\\1", CAZy_Info),
         Gene_Family = sub("_.*", "", CAZy_Info))

gh_data <- cazy_metadata %>% filter(CAZy_Class == "GH")

gh_summary <- gh_data %>%
  group_by(Gene_Family, Sample.Source) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Sample.Source, values_from = Count, values_fill = 0)

top20_gh <- gh_summary %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(Total)) %>%
  slice(1:20) %>%
  select(-Total)

gh_matrix <- as.matrix(top20_gh[, -1])
rownames(gh_matrix) <- top20_gh$Gene_Family
gh_matrix <- gh_matrix[, c("Abijata", "Chitu", "Shala", "Cattle", "Goat", "Sheep")]

# GH labels (manual or programmatic, your choice)
rownames(gh_matrix) <- c(
  "GH3 (β-glucosidase, arabinofuranosidase)",
  "GH43 (arabinofuranosidase, xylosidase)",
  "GH2 (β-galactosidase, mannosidase)",
  "GH13 (α-amylase, pullulanase)",
  "GH0 (Unclassified GH)",
  "GH23 (peptidoglycan hydrolase)",
  "GH5 (cellulase, xylanase, mannanase)",
  "GH1 (β-glucosidase, β-galactosidase)",
  "GH53 (β-1,4-galactanase)",
  "GH73 (β-N-acetylglucosaminidase)",
  "GH65 (trehalose/maltose phosphorylase)",
  "GH16 (xyloglucanase, β-glucanase)",
  "GH31 (α-glucosidase, α-xylosidase)",
  "GH20 (β-N-acetylhexosaminidase)",
  "GH171 (glucan hydrolase, rare)",
  "GH109 (α-N-acetylgalactosaminidase)",
  "GH28 (polygalacturonase, pectinase)",
  "GH51 (α-L-arabinofuranosidase)",
  "GH78 (α-L-rhamnosidase)",
  "GH97 (α-glucosidase, α-galactosidase)"
)

# Sample annotations
sample_meta_df <- tibble(
  SampleSource = colnames(gh_matrix),
  SampleType = c(rep("Soda Lake", 3), rep("Rumen", 3))
)

column_annot <- HeatmapAnnotation(
  SampleType = sample_meta_df$SampleType,
  SampleSource = sample_meta_df$SampleSource,
  col = list(
    SampleType = sample_type_palette,
    SampleSource = sample_source_palette
  ),
  show_annotation_name = FALSE,
  annotation_height = unit(c(4, 4), "mm")
)

gh_col_fun <- colorRamp2(c(0, max(gh_matrix)), c("white", "navy"))

ht_gh <- Heatmap(
  gh_matrix,
  name = "GH Gene Count",
  col = gh_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = column_annot,
  show_column_names = TRUE,
  show_row_names = TRUE,
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  width = unit(8, "cm"),
  heatmap_legend_param = list(title_position = "topcenter", legend_direction = "horizontal")
)

# Save heatmap
png("GH_Top20_Functional_Heatmap_Final_BottomLegend.png", width = 3000, height = 2200, res = 300)
draw(ht_gh, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# -------------------------------
# Save Combined Bar + Dot Panels
# -------------------------------
combined_horizontal <- plot_A + plot_B_dot + plot_C +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

ggsave("Combined_Functional_Panels.png", combined_horizontal, width = 18, height = 6, dpi = 600)
