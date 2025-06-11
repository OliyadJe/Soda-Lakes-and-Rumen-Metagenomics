#!/bin/bash

# -------------------------------------------------------------
# Batch DIAMOND blastp search against a CAZy database
# For protein files (*.faa) annotated from MAGs
# -------------------------------------------------------------


# Load and activate conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate /path/to/myenv  # Contains DIAMOND

# Set directories
WORKDIR="/path/to/cazy_inputs"                         # Folder with input .faa files
DB_DIR="/path/to/dbcan_db"                             # Folder with CAZyDB FASTA and DIAMOND DB
DB_FASTA="${DB_DIR}/CAZyDB.07142024.fa"
DB_DMND="${DB_DIR}/CAZy.dmnd"

cd "$WORKDIR"

# Check DIAMOND version
echo "🔍 DIAMOND version:"
diamond version

# Create DIAMOND database if not already present
if [[ ! -f "$DB_DMND" ]]; then
    echo "🔧 Building DIAMOND database from: $DB_FASTA"
    diamond makedb --in "$DB_FASTA" -d "$DB_DMND"
else
    echo "✅ DIAMOND database found: $DB_DMND"
fi

# Run DIAMOND blastp for each .faa file
for faa in *.faa; do
    output="${faa%.faa}_cazy_diamond_results.tsv"
    if [[ -f "$output" ]]; then
        echo "⏩ Skipping $faa — output already exists."
        continue
    fi

    echo "🚀 Processing $faa..."
    diamond blastp \
        -q "$faa" \
        -d "$DB_DMND" \
        -o "$output" \
        --outfmt 6 \
        --threads 32 || echo "❌ Error processing $faa" >> error.log
done

echo "✅ All DIAMOND CAZy analyses completed."
