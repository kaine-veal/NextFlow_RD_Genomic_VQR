/*
 * Trim reads with fastp
 */
process TrimReads {

    label 'process_low'

    container 'staphb/fastp:1.1.0'


    // Add a tag to identify the process
    tag "$sample_id"

    // Specify the output directory for trimmed reads and fastp reports
    publishDir("$params.outdir/trimmed_reads", mode: "copy")

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: html_report
    path "${sample_id}_fastp.json", emit: json_report

script:
"""
echo "Running fastp on sample: ${sample_id}"

fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${sample_id}_R1_trimmed.fastq.gz \
    -O ${sample_id}_R2_trimmed.fastq.gz \
    --detect_adapter_for_pe \
    --length_required 36 \
    -j ${sample_id}_fastp.json \
    -h ${sample_id}_fastp.html \
    -w ${task.cpus}

echo "fastp complete for ${sample_id}"
"""
}
