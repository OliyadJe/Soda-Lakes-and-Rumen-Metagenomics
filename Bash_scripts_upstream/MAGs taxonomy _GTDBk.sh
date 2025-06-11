#!/bin/bash

# ----------------------------------------------
# GTDB-Tk: Download DB and run taxonomic classification
# ----------------------------------------------

# Load necessary modules or conda environment
module load anaconda3
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate /path/to/gtdbtk_env

# ----------------------------
# Part 1: Database Download
# ----------------------------
DB_DIR="/path/to/gtdbtk_db"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

# Download and extract GTDB-Tk Release 220 (adjust release if needed)
wget -c https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
tar -xvzf gtdbtk_r220_data.tar.gz --strip 1
rm gtdbtk_r220_data.tar.gz

# Set the environment variable for GTDBTK
export GTDBTK_DATA_PATH="$DB_DIR"

# ----------------------------
# Part 2: Run GTDB-Tk classification
# ----------------------------
# Define input and output paths
INPUT_MAG_DIR="/path/to/high_quality_MAGs"           # Folder with MAGs (*.fa)
OUTPUT_DIR="/path/to/gtdbtk_output"
mkdir -p "$OUTPUT_DIR"

# Run GTDB-Tk classification workflow
gtdbtk classify_wf \
  --genome_dir "$INPUT_MAG_DIR" \
  --out_dir "$OUTPUT_DIR" \
  --cpus 16 \
  --extension fa

# Optional: summarize results
echo "🧾 Summary of classifications:"
cat "$OUTPUT_DIR/gtdbtk.bac120.summary.tsv" | head -n 5
echo "..."

# Final message
echo "✅ GTDB-Tk classification completed."
