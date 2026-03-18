process haplotypeCaller {

    label 'process_high'

    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$bamFile"

    input:
    tuple val(sample_id), file(bamFile), file(bamIndex)
    path indexFiles

    output:
    tuple val(sample_id), file("*.vcf"), file("*.vcf.idx")

    script:
    """
    echo "Running HaplotypeCaller for Sample: ${bamFile}"

    # Find the fasta file from staged index files
    genomeFasta=\$(find -L . -name "*.fasta")

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    outputVcf="\$(basename ${bamFile} _sorted_dedup_recalibrated.bam).vcf"

    gatk HaplotypeCaller -R "\${genomeFasta}" -I ${bamFile} -O "\${outputVcf}" -ERC GVCF \
        -A BaseQuality -A DepthPerSampleHC -A MappingQuality -A QualByDepth \
        -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A StrandOddsRatio \
        -A MappingQualityZero -A InbreedingCoeff -A BaseQualityRankSumTest -A HaplotypeFilteringAnnotation

    gatk IndexFeatureFile -I \${outputVcf}

    echo "Variant Calling for Sample: ${sample_id} Complete"
    """
}