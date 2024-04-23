process NORMALIZE_FILTER_VCF {
    input:
    tuple val(patient), val(tumour_sample), val(ffpe), val(normal_sample), path(vcf), path(tbi)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(patient), val(tumour_sample), val(ffpe), val(normal_sample), path("${tumour_sample}.norm.passonly.vcf.gz"), path("${tumour_sample}.norm.passonly.vcf.gz.tbi"), emit: normalized_filtered_vcf

    script:
    """
    bcftools norm -f ${genome} -m-any ${vcf} -o ${tumour_sample}.norm.vcf.gz -O z
    bcftools view -i 'FILTER="PASS" || FILTER="clustered_events"' -o ${tumour_sample}.norm.passonly.vcf.gz ${tumour_sample}.norm.vcf.gz -O z
    tabix -p vcf ${tumour_sample}.norm.passonly.vcf.gz
    """
}

process MERGE_VCF {
    input:
    tuple val(patient), val(tumour_samples), val(ffpe_status), val(normal_samples), path(vcf), path(tbi)
    path genome
    path genome_index
    path genome_dict

    publishDir "${params.outputdir}/${patient}/mergedvcf/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(tumour_samples), path("${patient}.tumour.merged.indel200.norm.vcf.gz"), path("${patient}.tumour.merged.indel200.norm.vcf.gz.tbi"), emit: merged_vcf
    
    script:
    """
    file_count=\$(ls -1 ./*.norm.passonly.vcf.gz 2>/dev/null | wc -l)

    echo \$file_count

    #if there is more than one sample for this patient, perform merge
    if [ "\$file_count" -gt 1 ]; then 
        bcftools merge ./*.norm.passonly.vcf.gz -o ${patient}.tumour.merged.vcf.gz -O z --force-samples
    else
        cp ./*.norm.passonly.vcf.gz ${patient}.tumour.merged.vcf.gz
    fi

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

process FUNCOTATOR {
    input:
    tuple val(patient), val(tumour_samples), path(vcf), path(tbi)
    path genome
    path genome_index
    path genome_dict
    val genome_version
    path funcotator_data

    publishDir "${params.outputdir}/${patient}/funcotator/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(tumour_samples), path("${patient}.tumour.merged.indel200.norm.anno.vcf.gz"), path("${patient}.tumour.merged.indel200.norm.anno.vcf.gz.tbi"), emit: annotated_vcf

    script:
    """
    gatk Funcotator \
        --output ./${patient}.tumour.merged.indel200.norm.anno.vcf.gz \
        --ref-version ${genome_version} \
        --data-sources-path ${funcotator_data} \
        --output-file-format VCF \
        --variant ${vcf} \
        --reference ${genome}
    """
}
