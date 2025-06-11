# Load required libraries
library(ape)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(dplyr)
library(tidyr)
library(data.table)
library(ggnewscale)
library(RColorBrewer)

# Load tree
tree <- read.tree("RAxML_bestTree.high_quality_bins_refined.tre")

# Load GTDB-Tk taxonomy and clean names
gtdb <- fread("gtdbtk.bac120.summary_latest.tsv") %>%
  separate(classification,
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right") %>%
  mutate(across(everything(), ~ sub("^[a-z]__", "", .)))

# Match GTDB to tree tips
tips <- tree$tip.label
gtdb <- gtdb %>% filter(user_genome %in% tips)
tree <- drop.tip(tree, setdiff(tips, gtdb$user_genome))

# Merge with metadata
gtdb <- gtdb %>%
  mutate(Sample = sub("_.*", "", user_genome)) %>%
  left_join(metadata, by = c("Sample" = "Sample.Names")) %>%
  mutate(
    Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum),
    Sample.Type = ifelse(is.na(Sample.Type), "Unknown", Sample.Type),
    Sample.Source = ifelse(is.na(Sample.Source), "Unknown", Sample.Source),
    Sample.Type = as.factor(Sample.Type),
    Sample.Source = as.factor(Sample.Source)
  )

phylum_colors <- c(
  Actinomycetota    = "#FF0000",  # Red
  Bacillota         = "#000000",  # Black
  Bacillota_A       = "#0000FF",  # Blue
  Bacillota_D       = "#FFFF00",  # Yellow
  Bacillota_I       = "#00FF00",  # Green
  Bacteroidota      = "#FF00FF",  # Magenta
  Bdellovibrionota  = "#A52A2A",  # Brown
  Chloroflexota     = "#808080",  # Gray
  Cyanobacteriota   = "#00FFFF",  # Cyan
  Desulfobacterota  = "#FFA500",  # Orange
  Elusimicrobiota   = "#800000",  # Maroon
  Fibrobacterota    = "#008000",  # Dark green
  Pseudomonadota    = "#800080",  # Purple
  Spirochaetota     = "#FFC0CB",  # Pink
  UBP6              = "#4682B4",  # Steel blue
  Verrucomicrobiota = "#B8860B"   # Dark goldenrod
)


scale_color_manual(values = phylum_colors, name = "Phylum") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5)))


# Plot circular tree
p <- ggtree(tree, layout = "circular") %<+% gtdb +
  geom_tree(size = 0.6, color = "black") +
  geom_tippoint(aes(color = Phylum), size = 1.2, alpha = 0.9) +
  scale_color_manual(values = phylum_colors, name = "Phylum") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme_tree2() +
  ggtitle("Phylogenetic Tree with Sample Metadata Rings") +
  theme(legend.position = "right")

# Highlight phylum clades (no labels)
phyla <- names(phylum_colors)
for (phylum in phyla) {
  tips_phylum <- gtdb$user_genome[gtdb$Phylum == phylum]
  if (length(tips_phylum) >= 2) {
    node <- getMRCA(tree, tips_phylum)
    p <- p + geom_hilight(node = node, fill = phylum_colors[phylum], alpha = 0.35, extendto = 1.15)
  }
}

# Add metadata rings (tight to tree)
p <- p +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    aes(y = label, fill = Sample.Type),
    width = 0.07, offset = 0.5, color = NA
  ) +
  scale_fill_brewer(palette = "Set2", name = "Sample Type") +
  
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    aes(y = label, fill = Sample.Source),
    width = 0.07, offset = 0.1, color = NA
  ) +
  scale_fill_brewer(palette = "Pastel2", name = "Sample Source")

# Save high-resolution PNG and PDF
ggsave("Final_Circular_Tree_SampleType_Source.png", plot = p, width = 16, height = 16, dpi = 600)
ggsave("Final_Circular_Tree_SampleType_Source.pdf", plot = p, width = 16, height = 16)
