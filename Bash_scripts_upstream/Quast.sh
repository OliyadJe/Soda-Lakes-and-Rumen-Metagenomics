#!/bin/bash

# ----------------------------------------------
# Run QUAST to evaluate metagenome assemblies
# across multiple samples (e.g., from MEGAHIT)
# ----------------------------------------------

# Load QUAST module or activate Conda environment
# Option 1 (module):
module load bioinfo-tools
module load quast

# Option 2 (Conda):
# source /path/to/miniconda3/etc/profile.d/conda.sh
# conda activate quast_env

# Set paths
INPUT_DIR="/path/to/megahit_output"         # Directory containing per-sample subfolders with final.contigs.fa
OUTPUT_DIR="/path/to/quast_results"
mkdir -p "$OUTPUT_DIR"

# Loop through each MEGAHIT output folder
for SAMPLE_DIR in "$INPUT_DIR"/*/; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    CONTIG_FILE="${SAMPLE_DIR}/final.contigs.fa"

    if [[ -f "$CONTIG_FILE" ]]; then
        echo "▶️ Running QUAST on $SAMPLE_NAME"

        quast.py "$CONTIG_FILE" \
            -o "${OUTPUT_DIR}/${SAMPLE_NAME}" \
            --threads 8 \
            --min-contig 500

        echo "✔️ QUAST complete: $SAMPLE_NAME"
    else
        echo "⚠️ Contig file not found for $SAMPLE_NAME — skipping."
    fi
done

echo "✅ All QUAST evaluations completed."
