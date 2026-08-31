#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * GSE174609 scRNA-seq preprocessing + QC workflow
 *
 * samples.tsv -> SRA download -> FASTQ.gz -> FastQC (R1 + R2) ---\
 *                                      \-> Cell Ranger count ------> BAM QC -> MultiQC
 *
 * Cell Ranger BAM creation is enabled.
 */

params.samples         = params.samples         ?: "${projectDir}/samples.tsv"
params.outdir          = params.outdir          ?: "${projectDir}/results"
params.transcriptome = params.transcriptome ?: "/opt/refdata-gex-GRCh38-2024-A"
params.cellranger    = params.cellranger    ?: "cellranger"
params.expected_cells  = params.expected_cells  ?: 5000
params.chemistry       = params.chemistry       ?: "SC3Pv3"
params.tools_container = params.tools_container ?: "${projectDir}/container/gse174609-tools.sif"


/* ---------------------------------------------------------
 * Download SRA, convert to FASTQ, compress, and rename
 * directly to the Cell Ranger-compatible naming convention.
 * --------------------------------------------------------- */
process SRA_TO_FASTQ {
    tag "${meta.sample} (${meta.srr})"

    container params.tools_container

    publishDir "${params.outdir}/fastq", mode: 'copy'

    cpus 4
    memory '8 GB'
    time '6h'

    input:
    val meta

    output:
    tuple val(meta),
          path("${meta.sample}_S${meta.sample_number}_L001_R1_001.fastq.gz"),
          path("${meta.sample}_S${meta.sample_number}_L001_R2_001.fastq.gz"),
          emit: reads

    script:
    def prefix = "${meta.sample}_S${meta.sample_number}_L001"
    """
    set -euo pipefail

    echo "Sample: ${meta.sample}"
    echo "SRA accession: ${meta.srr}"

    prefetch \
        --output-directory . \
        ${meta.srr}

    SRA_FILE="${meta.srr}/${meta.srr}.sra"

    if [[ ! -s \"\$SRA_FILE\" ]]; then
        echo "ERROR: SRA file was not downloaded: \$SRA_FILE" >&2
        exit 1
    fi

    fasterq-dump \
        --split-files \
        --include-technical \
        --threads ${task.cpus} \
        --outdir . \
        \"\$SRA_FILE\"

    if [[ ! -s "${meta.srr}_1.fastq" || ! -s "${meta.srr}_2.fastq" ]]; then
        echo "ERROR: Expected R1/R2 FASTQ files were not produced for ${meta.srr}" >&2
        exit 1
    fi

    # Compress R1 and R2. pigz uses all CPUs for each file sequentially,
    # avoiding two simultaneous compressors competing for the same cores.
    pigz -p ${task.cpus} "${meta.srr}_1.fastq"
    pigz -p ${task.cpus} "${meta.srr}_2.fastq"

    mv "${meta.srr}_1.fastq.gz" "${prefix}_R1_001.fastq.gz"
    mv "${meta.srr}_2.fastq.gz" "${prefix}_R2_001.fastq.gz"

    # Remove unused technical FASTQs and the downloaded SRA directory only
    # after R1/R2 have been produced successfully.
    rm -f "${meta.srr}"_*.fastq "${meta.srr}"_*.fastq.gz || true
    rm -rf "${meta.srr}"

    test -s "${prefix}_R1_001.fastq.gz"
    test -s "${prefix}_R2_001.fastq.gz"
    """
}


/* ---------------------------------------------------------
 * FastQC on BOTH R1 and R2.
 * R1 contains 10x cell barcode/UMI sequence, while R2 contains
 * the cDNA read, so their FastQC profiles are expected to differ.
 * --------------------------------------------------------- */
process FASTQC {
    tag "${meta.sample}"

    container params.tools_container

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    cpus 2
    memory '4 GB'
    time '1h'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: zip
    tuple val(meta), path("*_fastqc.html"), emit: html

    script:
    """
    set -euo pipefail

    fastqc \
        --threads ${task.cpus} \
        --quiet \
        ${r1} \
        ${r2}
    """
}


/* ---------------------------------------------------------
 * Run Cell Ranger count independently for every sample.
 * BAM creation is enabled so possorted_genome_bam.bam and its
 * index are generated under each sample's outs/ directory.
 * --------------------------------------------------------- */
process CELLRANGER_COUNT {
    tag "${meta.sample}"

    container params.tools_container

    publishDir "${params.outdir}/cellranger", mode: 'copy'

    cpus 12
    memory '64 GB'
    time '6h'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample}"), emit: results

    script:
    """
    set -euo pipefail

    ${params.cellranger} count \
        --id="${meta.sample}" \
        --fastqs="\$PWD" \
        --sample="${meta.sample}" \
        --transcriptome="${params.transcriptome}" \
        --expect-cells=${params.expected_cells} \
        --localcores=${task.cpus} \
        --localmem=64 \
        --chemistry="${params.chemistry}" \
        --create-bam=true

    test -s "${meta.sample}/outs/possorted_genome_bam.bam"
    """
}


/* ---------------------------------------------------------
 * BAM QC after Cell Ranger.
 * samtools flagstat: compact alignment summary.
 * samtools stats: detailed BAM/alignment statistics.
 * Both formats are recognized by MultiQC.
 * --------------------------------------------------------- */
process BAM_QC {
    tag "${meta.sample}"

    container params.tools_container

    publishDir "${params.outdir}/bam_qc", mode: 'copy'

    cpus 1
    memory '2 GB'
    time '1h'

    input:
    tuple val(meta), path(cellranger_dir)

    output:
    tuple val(meta),
          path("${meta.sample}.flagstat.txt"),
          path("${meta.sample}.samtools.stats.txt"),
          emit: qc

    script:
    def bam = "${cellranger_dir}/outs/possorted_genome_bam.bam"
    """
    set -euo pipefail

    test -s "${bam}"

    samtools flagstat \
        -@ ${task.cpus} \
        "${bam}" \
        > "${meta.sample}.flagstat.txt"

    samtools stats \
        -@ ${task.cpus} \
        "${bam}" \
        > "${meta.sample}.samtools.stats.txt"
    """
}


/* ---------------------------------------------------------
 * Aggregate raw-read FastQC and BAM-level samtools QC into
 * one final MultiQC report.
 * --------------------------------------------------------- */
process MULTIQC {
    tag "all samples"

    container params.tools_container

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    cpus 1
    memory '1 GB'
    time '30m'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    set -euo pipefail

    multiqc . \
        --filename multiqc_report.html \
        --title "GSE174609 Raw Read and BAM Quality Control" \
        --force
    """
}


/* ---------------------------------------------------------
 * Workflow
 * --------------------------------------------------------- */
workflow {

    /*
     * Convert each TSV row into a metadata map, for example:
     * [sample: Healthy_1, srr: SRR14575500, sample_number: 1]
     */
    samples_ch = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            [
                sample       : row.sample,
                srr          : row.srr,
                sample_number: row.sample_number as Integer
            ]
        }

    SRA_TO_FASTQ(samples_ch)

    /*
     * Raw-read QC and Cell Ranger can proceed independently once a
     * sample's R1/R2 FASTQs are ready.
     */
    FASTQC(SRA_TO_FASTQ.out.reads)
    CELLRANGER_COUNT(SRA_TO_FASTQ.out.reads)

    /* BAM QC starts only after the corresponding Cell Ranger run succeeds. */
    BAM_QC(CELLRANGER_COUNT.out.results)

    /*
     * FASTQC produces two ZIP files per sample (R1 and R2), so flatten
     * those file lists. Then add both samtools outputs for every sample.
     * MULTIQC receives the complete set only after all QC jobs finish.
     */
    fastqc_files_ch = FASTQC.out.zip
        .map { meta, files -> files }
        .flatten()

    bam_qc_files_ch = BAM_QC.out.qc
        .flatMap { meta, flagstat, stats -> [flagstat, stats] }

    all_qc_files_ch = fastqc_files_ch
        .mix(bam_qc_files_ch)
        .collect()

    MULTIQC(all_qc_files_ch)
}
