#!/bin/bash


# ----------------------------------------------
# Script to run MEGAHIT for metagenome assembly
# on paired-end kneaddata-filtered reads
# ----------------------------------------------

# Load Conda (adjust for your system)
module load anaconda3
conda activate myenv  # Replace with your actual environment name

# Define input and output directories
INPUT_DIR="/path/to/kneaddata_output"
OUTPUT_DIR="/path/to/megahit_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Automatically detect sample prefixes from paired-end files
SAMPLES=($(ls "$INPUT_DIR"/*_kneaddata_paired_1.fastq | sed 's/_kneaddata_paired_1.fastq//' | xargs -n 1 basename))

# Loop through all samples
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${INPUT_DIR}/${SAMPLE}_kneaddata_paired_1.fastq"
    R2="${INPUT_DIR}/${SAMPLE}_kneaddata_paired_2.fastq"
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"

    # Skip if one of the pair is missing
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "⚠️ Skipping $SAMPLE — one or both paired files not found."
        continue
    fi

    # Remove and recreate output dir
    rm -rf "$SAMPLE_OUTPUT_DIR"
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    echo "▶️ Running MEGAHIT on $SAMPLE..."
    megahit \
        -1 "$R1" \
        -2 "$R2" \
        --min-contig-len 500 \
        --num-cpu-threads 4 \
        -o "$SAMPLE_OUTPUT_DIR"

    echo "✔️ Assembly complete: $SAMPLE"
done

# Deactivate Conda environment
conda deactivate

echo "✅ All assemblies completed."
