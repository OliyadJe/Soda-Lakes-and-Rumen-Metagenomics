
# Load modules
module load anaconda3
module load PDC  # Optional, if your HPC requires it

# Activate conda environment (contains samtools, metabat2, parallel)
conda activate myenv

# Set base directory containing MEGAHIT results
MEGAHIT_OUTPUT_DIR="./megahit_output"

# Function: Convert SAM to sorted BAM
convert_sam_to_bam() {
    SAMPLE_DIR=$1
    SAMPLE=$(basename "$SAMPLE_DIR")

    if [ -f "${SAMPLE_DIR}/mapped_reads.sam" ]; then
        echo "[$SAMPLE] Converting SAM to BAM..."
        samtools view -bS "${SAMPLE_DIR}/mapped_reads.sam" | samtools sort -o "${SAMPLE_DIR}/mapped_reads_sorted.bam"
        samtools index "${SAMPLE_DIR}/mapped_reads_sorted.bam"
    else
        echo "[$SAMPLE] SAM not found. Skipping..."
    fi
}

# Function: Generate depth file
generate_depth_file() {
    SAMPLE_DIR=$1
    if [ -f "${SAMPLE_DIR}/mapped_reads_sorted.bam" ]; then
        echo "[$SAMPLE_DIR] Generating depth file..."
        jgi_summarize_bam_contig_depths --outputDepth "${SAMPLE_DIR}/depth.txt" "${SAMPLE_DIR}/mapped_reads_sorted.bam"
    else
        echo "[$SAMPLE_DIR] BAM not found. Skipping..."
    fi
}

# Function: Run MetaBAT2
run_metabat2() {
    SAMPLE_DIR=$1
    CONTIGS="${SAMPLE_DIR}/final.contigs.fa"
    DEPTH="${SAMPLE_DIR}/depth.txt"
    BINS="${SAMPLE_DIR}/metabat_bins"

    if [ -f "$CONTIGS" ] && [ -f "$DEPTH" ]; then
        echo "[$SAMPLE_DIR] Running MetaBAT2..."
        mkdir -p "$BINS"
        metabat2 -i "$CONTIGS" -a "$DEPTH" -o "${BINS}/bin"
    else
        echo "[$SAMPLE_DIR] Missing input files. Skipping..."
    fi
}

# Export functions for GNU Parallel
export -f convert_sam_to_bam generate_depth_file run_metabat2

# Run all steps in parallel
find "${MEGAHIT_OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -type d | parallel --jobs 4 convert_sam_to_bam
find "${MEGAHIT_OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -type d | parallel --jobs 4 generate_depth_file
find "${MEGAHIT_OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -type d | parallel --jobs 4 run_metabat2
