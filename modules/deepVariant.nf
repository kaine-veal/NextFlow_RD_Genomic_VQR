process deepVariant {

    label 'process_high'

    container 'google/deepvariant:1.6.0'

    tag "$sample_id"

    publishDir("$params.outdir/VCF/deepvariant", mode: "copy")

    input:
    tuple val(sample_id), path(bam), path(bai)
    path indexFiles

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.g.vcf.gz")

    script:
    """
    run_deepvariant \
    --model_type=WES \
    --ref=\$(find -L . -name '*.fasta') \
    --reads=${bam} \
    --output_vcf=${sample_id}.vcf.gz \
    --output_gvcf=${sample_id}.g.vcf.gz \
    --num_shards=${task.cpus}
    """
}