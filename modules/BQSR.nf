process baseRecalibrator {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'broadinstitute/gatk:4.1.4.0'

    tag "$bamFile"

    publishDir("$params.outdir/BAM", mode: "copy")
    publishDir("$params.outdir/BQSR_QC", mode: "copy", pattern: "*.recal_data.table")
    publishDir("$params.outdir/BQSR_QC", mode: "copy", pattern: "*.post_recal_data.table")
    publishDir("$params.outdir/BQSR_QC", mode: "copy", pattern: "*.BQSR.before_after.pdf")

    input:
    tuple val(sample_id), file(bamFile), file(baiFile)
    val knownSites
    path indexFiles
    path qsrcVcfFiles

    output:
    tuple val(sample_id),
          file("${bamFile.baseName}_recalibrated.bam"),
          file("${bamFile.baseName}_recalibrated.bai"),
          file("${bamFile.baseName}.recal_data.table"),
          file("${bamFile.baseName}.post_recal_data.table"),
          file("${bamFile.baseName}.BQSR.before_after.pdf")

    script:
    def knownSitesArgs = knownSites.join(' ')
    """
    echo "Running BQSR"

    if [[ -n "${params.genome_file}" ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta' | head -1)
    fi

    echo "Genome File: \${genomeFasta}"

    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    # 1. First pass: generate pre-BQSR recalibration table
    gatk --java-options "-Xmx8G" BaseRecalibrator \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        ${knownSitesArgs} \
        -O ${bamFile.baseName}.recal_data.table

    # 2. Apply BQSR
    gatk --java-options "-Xmx8G" ApplyBQSR \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        --bqsr-recal-file ${bamFile.baseName}.recal_data.table \
        -O ${bamFile.baseName}_recalibrated.bam

    # 3. Index recalibrated BAM
    samtools index ${bamFile.baseName}_recalibrated.bam ${bamFile.baseName}_recalibrated.bai

    # 4. Second pass: generate post-BQSR recalibration table
    gatk --java-options "-Xmx8G" BaseRecalibrator \
        -R "\${genomeFasta}" \
        -I ${bamFile.baseName}_recalibrated.bam \
        ${knownSitesArgs} \
        -O ${bamFile.baseName}.post_recal_data.table

    # 5. Create before/after diagnostic plots
    gatk AnalyzeCovariates \
        -before ${bamFile.baseName}.recal_data.table \
        -after ${bamFile.baseName}.post_recal_data.table \
        -plots ${bamFile.baseName}.BQSR.before_after.pdf

    echo "BQSR Complete"
    """
}