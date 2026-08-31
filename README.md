# GSE174609 scRNA-seq Nextflow pipeline

## Workflow

```text
samples.tsv
    |
    v
SRA_TO_FASTQ
    |--------------------> FASTQC (R1 + R2) ----\
    |                                            \
    +--------------------> CELLRANGER_COUNT       >----> MULTIQC
                                 |                 /
                                 v                /
                              BAM_QC -------------/
                         (flagstat + stats)
```

The workflow now performs QC at two levels:

1. **Raw reads:** FastQC on both R1 and R2.
2. **Cell Ranger BAM:** `samtools flagstat` and `samtools stats`.

All of these results are aggregated into one final MultiQC report.

Cell Ranger is run with BAM creation enabled (`--create-bam=true`).

## Files

- `main.nf` — Nextflow DSL2 workflow
- `nextflow.config` — paths, SLURM profile and execution provenance settings
- `samples.tsv` — biological sample / SRA mapping
- `envs/scrna.yml` — Conda environment containing SRA Tools, FastQC, MultiQC, pigz and samtools

Cell Ranger is intentionally not installed through Conda. The workflow points to the existing Cell Ranger 10 installation.

## Main outputs

```text
results/
├── fastq/
│   ├── Healthy_1_S1_L001_R1_001.fastq.gz
│   ├── Healthy_1_S1_L001_R2_001.fastq.gz
│   └── ...
│
├── fastqc/
│   ├── Healthy_1_S1_L001_R1_001_fastqc.html
│   ├── Healthy_1_S1_L001_R1_001_fastqc.zip
│   ├── Healthy_1_S1_L001_R2_001_fastqc.html
│   ├── Healthy_1_S1_L001_R2_001_fastqc.zip
│   └── ...
│
├── cellranger/
│   ├── Healthy_1/
│   │   └── outs/
│   │       ├── possorted_genome_bam.bam
│   │       ├── possorted_genome_bam.bam.bai
│   │       ├── filtered_feature_bc_matrix.h5
│   │       ├── web_summary.html
│   │       └── ...
│   └── ...
│
├── bam_qc/
│   ├── Healthy_1.flagstat.txt
│   ├── Healthy_1.samtools.stats.txt
│   └── ...
│
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
│
└── pipeline_info/
    ├── execution_trace.txt
    ├── execution_report.html
    ├── execution_timeline.html
    └── pipeline_dag.html
```

## Why FastQC is run on both reads

For 10x 3' single-cell RNA-seq, R1 and R2 serve different purposes:

- **R1:** cell barcode and UMI sequence.
- **R2:** cDNA/transcript sequence.

Their FastQC profiles are therefore expected to differ substantially. In particular, ordinary assumptions about sequence composition and duplication for genomic FASTQ files should not automatically be applied to R1.

## Check configuration

```bash
nextflow config -profile arc
```

## Run on ARC / SLURM

From the pipeline directory:

```bash
nextflow run main.nf \
    -profile arc \
    -resume
```

`-resume` allows Nextflow to reuse successfully completed process outputs after a failed or interrupted run. Do not delete the Nextflow `work/` directory if you want `-resume` to work.

## Test one sample first

For an HPC workflow containing SRA downloads and Cell Ranger, testing one sample end-to-end before launching the complete cohort is recommended. You can temporarily create a one-row `samples_test.tsv` and run:

```bash
nextflow run main.nf \
    -profile arc \
    --samples samples_test.tsv \
    --outdir results_test \
    -resume
```

## Override parameters

```bash
nextflow run main.nf \
    -profile arc \
    --samples samples.tsv \
    --outdir results \
    --expected_cells 5000 \
    --chemistry SC3Pv3 \
    -resume
```

## Apptainer / Singularity container

This version of the workflow includes a container definition for the open-source tools used by the pipeline:

- SRA Toolkit (`prefetch`, `fasterq-dump`)
- FastQC
- MultiQC
- samtools
- pigz

Cell Ranger is **not** bundled in the image. The pipeline continues to use the existing 10x Genomics installation specified by `params.cellranger`, for example:

```text
/work/IMC_binf/sbagheri/sc_RNA_seq/software/cellranger-10.0.0/cellranger
```

### Build the image on ARC

```bash
cd container
./build_image.sh
```

This creates:

```text
container/gse174609-tools.sif
```

You can test it with:

```bash
./test_image.sh
```

If ARC does not allow `apptainer build --fakeroot`, build the image on a system where Apptainer/Singularity builds are permitted and copy the resulting `.sif` file into `container/`.

### Run the workflow

From the project directory:

```bash
nextflow run main.nf \
    -profile arc \
    -resume
```

The `arc` run uses SLURM and the containerized open-source processes. `CELLRANGER_COUNT` intentionally runs outside the container so that it can use the existing Cell Ranger installation directly.

The configuration binds `/work`, `/home`, and `/tmp` into the container so the workflow can access ARC project data, references, Nextflow work directories, and temporary files.

### Run without the container

For debugging on a machine where all tools are installed directly:

```bash
nextflow run main.nf \
    -profile host \
    -resume
```
