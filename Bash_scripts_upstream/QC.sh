#!/bin/bash

# -------------------------------------------------------
# Batch script to run KneadData on multiple paired-end 
# metagenomic samples using Conda and SLURM-compatible clusters
# -------------------------------------------------------

# SLURM directives (adjust for your own cluster if needed)
#SBATCH --job-name=kneaddata_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=64G

# Load required modules (customize to your HPC module system)
module load bioinfo-tools
module load bowtie2
module load java
module load python
module load trimmomatic
module load FastQC
module load MultiQC

# Activate KneadData environment
# Replace with your actual environment path or use `conda activate`
source /path/to/kneaddata_env/bin/activate

# Define paths (customize for your dataset)
input_dir="/path/to/raw_reads"
output_dir="/path/to/kneaddata_output"
database_dir="/path/to/kneaddata_db"

# Create output directory (reset if it already exists)
rm -rf $output_dir
mkdir -p $output_dir

# Get list of sample base names (assumes *_1.fq.gz and *_2.fq.gz)
samples=($(ls ${input_dir}/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//' | sort -u))

# Process samples in batches
batch_size=5

for ((start_idx = 0; start_idx < ${#samples[@]}; start_idx += batch_size)); do
  end_idx=$((start_idx + batch_size - 1))
  if [[ $end_idx -ge ${#samples[@]} ]]; then
    end_idx=$((${#samples[@]} - 1))
  fi

  for ((i = start_idx; i <= end_idx; i++)); do
    sample=${samples[$i]}
    input1="${input_dir}/${sample}_1.fq.gz"
    input2="${input_dir}/${sample}_2.fq.gz"
    kneaddata_output="${output_dir}/${sample}"

    mkdir -p $kneaddata_output

    kneaddata \
      --input1 $input1 \
      --input2 $input2 \
      --reference-db $database_dir \
      --output $kneaddata_output \
      --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
      --sequencer-source TruSeq3 \
      --run-trim-repetitive \
      --run-fastqc-end \
      --bypass-trf \
      --log-level INFO &
  done

  wait  # Wait for current batch to finish
done

done

echo "Analysis complete."
