process variantRecalibrator {

    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'
    tag "$vcf"

    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcf)
    val knownSitesArgs
    path genome
    path qsrc_vcf

    output:
    tuple val(sample_id), file("${vcf.baseName}.recalibrated.final.vcf")

    script:
    def knownSitesArgsStr = knownSitesArgs.join(' ')
    def degradedDna = params.degraded_dna == "true"

    """
    echo "Running VQSR"

    genomeFasta=\$(find -L . -name '*.fasta' | head -1)

    # Create dict if not present
    if [ ! -f "\${genomeFasta%.fasta}.dict" ]; then
        gatk CreateSequenceDictionary -R "\${genomeFasta}"
    fi

    # Index the VCF if index not present
    if [ ! -f "${vcf}.tbi" ]; then
        gatk IndexFeatureFile -I "${vcf}"
    fi

    echo "Genome File: \${genomeFasta}"

    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    if ${degradedDna}; then
        echo "Running VQSR for degraded DNA (1x coverage)"
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an FS -an SOR \
            -mode SNP \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -O ${vcf.baseName}.recalibrated_SNP.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -recal-file ${vcf.baseName}.recalibrated_SNP.recal \
            -mode SNP \
            -O ${vcf.baseName}.output_SNP.vcf
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an FS -an SOR \
            -mode INDEL \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -O ${vcf.baseName}.recalibrated_INDEL.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -recal-file ${vcf.baseName}.recalibrated_INDEL.recal \
            -mode INDEL \
            -O ${vcf.baseName}.output_INDEL.vcf
    else
        echo "Running VQSR for standard DNA (10x+ coverage)"
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -O ${vcf.baseName}.recalibrated_SNP.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -recal-file ${vcf.baseName}.recalibrated_SNP.recal \
            -mode SNP \
            -O ${vcf.baseName}.output_SNP.vcf
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
            -mode INDEL \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -O ${vcf.baseName}.recalibrated_INDEL.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -recal-file ${vcf.baseName}.recalibrated_INDEL.recal \
            -mode INDEL \
            -O ${vcf.baseName}.output_INDEL.vcf
    fi

    gatk MergeVcfs \
        -I ${vcf.baseName}.output_SNP.vcf \
        -I ${vcf.baseName}.output_INDEL.vcf \
        -O ${vcf.baseName}.recalibrated.vcf

    # Remove unevaluated variants (FILTER='.') before hap.py
    bcftools view -i 'FILTER!="."' \
        ${vcf.baseName}.recalibrated.vcf \
        -O v -o ${vcf.baseName}.recalibrated.final.vcf

    echo "VQSR Complete"
    """
}