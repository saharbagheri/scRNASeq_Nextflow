echo "=== FastQC ==="
apptainer exec container/gse174609-tools.sif fastqc --version

echo "=== MultiQC ==="
apptainer exec container/gse174609-tools.sif multiqc --version

echo "=== samtools ==="
apptainer exec container/gse174609-tools.sif samtools --version | head -1

echo "=== fasterq-dump ==="
apptainer exec container/gse174609-tools.sif fasterq-dump --version

echo "=== prefetch ==="
apptainer exec container/gse174609-tools.sif prefetch --version

echo "=== pigz ==="
apptainer exec container/gse174609-tools.sif pigz --version

echo "=== Cell Ranger ==="
apptainer exec container/gse174609-tools.sif cellranger --version

