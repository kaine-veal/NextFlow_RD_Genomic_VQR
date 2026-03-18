process filterVCF {
    label 'process_low'
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$vcfFile"

    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcfFile), file(vcfIndex)
    path indexFiles

    output:
    tuple val(sample_id), file("*_filtered.vcf")

    script:
    def isDegradedDNA = params.degraded_dna ? 'true' : 'false'
    """
    genomeFasta=\$(find -L . -name '*.fasta' | head -1)

    # Create dict if not staged
    if [ ! -f "\${genomeFasta%.fasta}.dict" ]; then
        gatk CreateSequenceDictionary -R "\${genomeFasta}"
    fi

    # Index the VCF if index not staged
    if [ ! -f "${vcfFile}.tbi" ]; then
        gatk IndexFeatureFile -I "${vcfFile}"
    fi

    baseName=\$(basename ${vcfFile} .gz)
    baseName=\$(basename \${baseName} .vcf)
    outputVcf="\${baseName}_filtered.vcf"

    if [ "$isDegradedDNA" == "true" ]; then
        echo "Running variant filtration for degraded DNA (2 x coverage)"
        gatk VariantFiltration -R "\${genomeFasta}" -V "${vcfFile}" -O "\${outputVcf}" \
            --filter-name "LowQUAL" --filter-expression 'QUAL < 150.0' \
            --filter-name "LowQD" --filter-expression 'QD < 2.0' \
            --filter-name "LowCoverage" --filter-expression 'DP < 5' \
            --filter-name "HighFS" --filter-expression 'FS > 60.0' \
            --filter-name "HighSOR" --filter-expression 'SOR > 3.0' \
            --filter-name "LowMQ" --filter-expression 'MQ < 60.0' \
            --genotype-filter-name "LowGQ" --genotype-filter-expression 'GQ < 30' \
            --set-filtered-genotype-to-no-call
    else
        echo "Running variant filtration for standard DNA (10x+ coverage)"
        gatk VariantFiltration -R "\${genomeFasta}" -V "${vcfFile}" -O "\${outputVcf}" \
            --filter-name "LowQUAL" --filter-expression 'QUAL < 500.0' \
            --filter-name "LowQD" --filter-expression 'QD < 2.0' \
            --filter-name "LowCoverage" --filter-expression 'DP < 30' \
            --filter-name "HighFS" --filter-expression 'FS > 60.0' \
            --filter-name "HighSOR" --filter-expression 'SOR > 3.0' \
            --filter-name "LowMQ" --filter-expression 'MQ < 60.0' \
            --genotype-filter-name "LowGQ" --genotype-filter-expression 'GQ < 30' \
            --set-filtered-genotype-to-no-call
    fi

    echo "Variant Filtering for Sample: ${vcfFile} Complete"
    """
}