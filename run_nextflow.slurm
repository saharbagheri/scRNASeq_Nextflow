#!/usr/bin/env bash
#SBATCH --job-name=nextflow_driver
#SBATCH --partition=synergy
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=nextflow_runner_%j.out
#SBATCH --error=nextflow_runner_%j.err

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow_env

cd /work/IMC_binf/sbagheri/sc_RNA_seq_nextflow

echo "Pipeline started at: $(date)"

nextflow run main.nf \
    -profile arc \
    -resume \
    > pipeline_run.log 2>&1

EXIT_CODE=$?

echo "Pipeline finished at: $(date)"
echo "Exit code: $EXIT_CODE"

exit "$EXIT_CODE"




