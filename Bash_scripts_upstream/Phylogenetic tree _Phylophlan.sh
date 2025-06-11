#!/bin/bash

# ----------------------------------------------
# Run PhyloPhlAn to generate a phylogenetic tree
# from high-quality metagenome-assembled genomes (MAGs)
# ----------------------------------------------

# Load Conda and activate PhyloPhlAn environment
module load anaconda3
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate /path/to/phylophlan_env

# Set paths
INPUT_DIR="/path/to/high_quality_MAGs"          # Folder containing MAGs (*.fa)
OUTPUT_DIR="${INPUT_DIR}/phylophlan_results"
CONFIG_FILE="/path/to/phylophlan_config/supermatrix_aa.cfg"
DATABASES_FOLDER="/path/to/phylophlan_databases"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "🧬 Starting PhyloPhlAn analysis at $(date)"

# Run PhyloPhlAn
phylophlan \
  -i "$INPUT_DIR" \
  -d phylophlan \
  -o "$OUTPUT_DIR" \
  --diversity high \
  --fast \
  --nproc 32 \
  -f "$CONFIG_FILE" \
  --genome_extension .fa \
  --databases_folder "$DATABASES_FOLDER"

echo "✅ PhyloPhlAn analysis completed at $(date)"

# Deactivate Conda environment
conda deactivate
