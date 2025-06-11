#!/bin/bash

# ----------------------------------------------
# Batch script to run HUMAnN 3 on paired-end
# kneaddata-cleaned metagenomic reads
# ----------------------------------------------

# Load conda and activate HUMAnN3 environment (adjust to your system)
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate /path/to/humann3_env

# Define input, output, and database directories
INPUT_DIR="/path/to/kneaddata_output"
OUTPUT_DIR="/path/to/humann3_output"
DATABASE_DIR="/path/to/humann3_databases"  # Should include 'chocophlan' and 'uniref'

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through paired files
for FORWARD_FILE in "$INPUT_DIR"/*_kneaddata_paired_1.fastq; do
    SAMPLE=$(basename "$FORWARD_FILE" "_kneaddata_paired_1.fastq")
    REVERSE_FILE="${INPUT_DIR}/${SAMPLE}_kneaddata_paired_2.fastq"
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"

    # Check if both paired reads exist and if output already exists
    if [[ -f "$FORWARD_FILE" && -f "$REVERSE_FILE" ]]; then
        if [[ -d "$SAMPLE_OUTPUT_DIR" && -f "$SAMPLE_OUTPUT_DIR/${SAMPLE}_genefamilies.tsv" ]]; then
            echo "✔️ Output for $SAMPLE already exists. Skipping."
        else
            echo "▶️ Processing $SAMPLE"

            humann \
                --input "$FORWARD_FILE" \
                --input "$REVERSE_FILE" \
                --output "$SAMPLE_OUTPUT_DIR" \
                --nucleotide-database "${DATABASE_DIR}/chocophlan" \
                --protein-database "${DATABASE_DIR}/uniref" \
                --threads 32
        fi
    else
        echo "⚠️ Warning: Missing paired files for $SAMPLE. Skipping."
    fi
done

echo "✅ HUMAnN3 analysis completed."
