process alignReadsBowtie2 {
    label 'process_low'
    container 'staphb/bowtie2:2.4.4'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path requiredIndexFiles

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    INDEX=\$(find -L ./ -name "*.bt2" | sed 's/\\.[0-9]*\\.bt2\$//' | sed 's/\\.rev//' | sort -u | head -1)

    bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} | \
    samtools view -b - | \
    samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
    """
}