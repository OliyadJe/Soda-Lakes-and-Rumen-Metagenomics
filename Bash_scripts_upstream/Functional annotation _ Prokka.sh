# Load necessary modules
module load anaconda3/2023.02-1         # Load Anaconda module (adjust version as needed)

# Activate the Prokka environment
source activate /cfs/klemming/projects/snic/naiss2023-23-323/tools/prokka_env

# Export the number of threads for Prokka
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Define the project path
PROJECT_DIR="/cfs/klemming/projects/snic/naiss2023-23-323"
BIN_DIR="${PROJECT_DIR}/rumen_sodalak/megahit_output/min/high_quality_bins"

# Create Prokka output directory if it doesn't exist
mkdir -p ${BIN_DIR}/prokka_annotations

# Run Prokka on all bins, skipping those that are already processed
for bin in ${BIN_DIR}/*.fa; do
    bin_name=$(basename "$bin" .fa)
    output_dir="${BIN_DIR}/prokka_annotations/$bin_name"
    # Check if the annotation is already completed
    if [ -f "$output_dir/$bin_name.gff" ]; then
        echo "Prokka annotation for $bin_name already exists. Skipping."
        continue
    fi
    mkdir -p "$output_dir"
    echo "Running Prokka annotation for $bin_name..."
    prokka --outdir "$output_dir" --prefix "$bin_name" --cpus $SLURM_CPUS_PER_TASK "$bin"
done

# Deactivate the environment
conda deactivate
