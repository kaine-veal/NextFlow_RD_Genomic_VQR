process filterdeepvariantVCF {
    label 'process_low'
    container 'staphb/bcftools:1.23'
    tag "$sample_id"

    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("dv_filtered.vcf.gz"), path("dv_filtered.vcf.gz.tbi")

    script:
    def min_qual = params.degraded_dna ? 5  : 10
    def min_dp   = params.degraded_dna ? 2  : 2
    """
    bcftools filter \
        -e "FMT/GQ < ${min_qual} || FMT/DP < ${min_dp} || (GT='het' && FMT/DP > 0 && (FMT/AD[0:1]/FMT/DP) < 0.2)" \
        -s "LowConf" \
        -m + \
        -O z \
        -o dv_filtered.vcf.gz \
        ${vcf}

    bcftools index -t dv_filtered.vcf.gz
    """
}