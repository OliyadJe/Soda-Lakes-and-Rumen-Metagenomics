#!/bin/bash

# -------------------------------------------------------
# Batch script to annotate proteins using EcoFoldDB
# Requires ProstT5, and FoldSeek binaries
# -------------------------------------------------------

# Load Conda and activate environment (with GPU support)
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate nvcc  # Replace with your conda env name

# Set input and output paths
INPUT_DIR="/path/to/faa_files"                            # Directory with *.faa input protein files
OUTDIR="${INPUT_DIR}/EcoFoldDB_Outputs"
ECOFOLDB_DIR="/path/to/EcoFoldDB/EcoFoldDB_v1.3"
PROSTT5_DIR="/path/to/EcoFoldDB/ProstT5_dir"
FOLDSEEK_BIN="/path/to/EcoFoldDB/foldseek/bin"

mkdir -p "$OUTDIR"
cd "$INPUT_DIR"

# Loop through each protein fasta (.faa) file
for faa in *.faa; do
    base=$(basename "$faa" .faa)
    sample_outdir="${OUTDIR}/${base}_out"

    echo "▶️ Running EcoFoldDB annotation for: $base"
    rm -rf "$sample_outdir"

    "${ECOFOLDB_DIR}/EcoFoldDB-annotate.sh" \
        --EcoFoldDB_dir "$ECOFOLDB_DIR" \
        --ProstT5_dir "$PROSTT5_DIR" \
        --gpu 0 \
        --prefilter-mode 0 \
        -e 1e-19 \
        --qcov 0.8 \
        --tcov 0.8 \
        --foldseek_bin "$FOLDSEEK_BIN" \
        -o "$sample_outdir" \
        "$faa"
done

echo "✅ EcoFoldDB annotation completed for all samples."
