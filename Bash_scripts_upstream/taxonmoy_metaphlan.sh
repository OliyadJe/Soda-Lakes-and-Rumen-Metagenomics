#!/bin/bash

# ----------------------------------------------
# MetaPhlAn 4 batch script for paired-end reads
# Automatically processes all *_1.fq.gz + *_2.fq.gz pairs
# ----------------------------------------------

# Load required modules (adjust as needed)
module load bioinfo-tools
module load MetaPhlAn4
module load bowtie2

# Define input/output directories and database path
INPUT_DIR="/path/to/host_removed_reads"
OUTPUT_DIR="/path/to/metaphlan_output"
DB_DIR="/path/to/metaphlan_db"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Automatically find sample base names (without _1/_2 suffix)
samples=($(ls ${INPUT_DIR}/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//' | sort -u))

# Loop through all samples
for SAMPLE in "${samples[@]}"; do
    R1="${INPUT_DIR}/${SAMPLE}_1.fq.gz"
    R2="${INPUT_DIR}/${SAMPLE}_2.fq.gz"

    metaphlan "${R1},${R2}" \
      --input_type fastq \
      --bowtie2db "${DB_DIR}" \
      --bowtie2out "${OUTPUT_DIR}/${SAMPLE}_bowtie2.bz2" \
      --nproc 8 \
      -o "${OUTPUT_DIR}/${SAMPLE}_profile.txt"

    echo "✔️ Completed: ${SAMPLE}"
done

echo "✅ All MetaPhlAn taxonomic profiling complete."
