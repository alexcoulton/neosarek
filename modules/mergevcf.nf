process MERGE_VCF {
    input:
    tuple val(patient), val(tumour_samples), val(normal_samples), path(vcf_paths), path(vcf_tbi_paths)
    path genome
    path genome_index
    path genome_dict

    publishDir "${params.outputdir}/${patient}/mergedvcf/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(tumour_samples), path("${patient}.tumour.merged.indel200.norm.vcf.gz"), path("${patient}.tumour.merged.indel200.norm.vcf.gz.tbi"), emit: merged_vcf
    
    script:
    """
    for vcf in ./*.vcf.gz; do
        output="\${vcf%.vcf.gz}.norm.vcf.gz"
        bcftools norm -f ${genome} -m-any "\$vcf" -o "\$output" -O z
        bcftools view -i 'FILTER="PASS" || FILTER="clustered_events"' -o \$output.passonly.vcf.gz \$output -O z
        tabix -p vcf \$output.passonly.vcf.gz
    done

    bcftools merge ./*.norm.vcf.gz.passonly.vcf.gz -o ${patient}.tumour.merged.vcf.gz -O z --force-samples

    #mutect2 force calling will not work with very large indels
    #limit these to 200 in length:
    bcftools view -i '(ILEN >= -200 && ILEN <= 200) || TYPE!="INDEL"' ${patient}.tumour.merged.vcf.gz -Oz -o ${patient}.tumour.merged.indel200.vcf.gz
    bcftools norm -f ${genome} -m-any ${patient}.tumour.merged.indel200.vcf.gz -o ${patient}.tumour.merged.indel200.norm.vcf.gz -O z
    tabix -p vcf ${patient}.tumour.merged.indel200.norm.vcf.gz
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/mergedvcf/* ./
    """
}
