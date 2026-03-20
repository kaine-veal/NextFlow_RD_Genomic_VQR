process alignReadsBowtie2 {
    container 'staphb/bowtie2:2.4.4'
    
    input:
    tuple val(sample_id), path(reads)
    path requiredIndexFiles

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    INDEX=\$(find -L ./ -name "*.bt2" | sed 's/\\.[0-9]*\\.bt2\$//' | sed 's/\\.rev//' | sort -u | head -1)
    
    bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} > ${sample_id}.sam
    """
}

process samToBam {

    container 'variantvalidator/indexgenome:1.1.0'

    tag "$sample_id"

    publishDir("$params.outdir/aligned_reads", mode: "copy")

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    samtools view -b ${sam} |
    samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
    """
}