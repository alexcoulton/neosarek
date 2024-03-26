process MUTECT2_FORCECALL {
    input:
    tuple val(patient), val(normal_sample), path(normal_cram), path(tumour_merged_vcf), path(tumour_merged_tbi), val(tumour_sample), path(tumour_cram)
    path germline_resource
    path germline_resource_index
    path panel_of_normals
    path panel_of_normals_index
    path genome
    path genome_index
    path genome_dict

    publishDir "${params.outputdir}/${patient}/tumour_samples/${tumour_sample}/mutect2/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite


    output:
    tuple val(patient), val(tumour_sample), path("${tumour_sample}.forcecall.norm.vcf.gz"), emit: mutect2_vcf
    
    script:
    """
    gatk Mutect2 \
        -R ${genome} \
        -I ${tumour_cram} \
        -I ${normal_cram} \
        -normal ${patient}_${normal_sample} \
        --germline-resource ${germline_resource} \
        --panel-of-normals ${panel_of_normals} \
        --f1r2-tar-gz f1r2.tar.gz \
        -alleles ${tumour_merged_vcf} \
        --force-call-filtered-alleles \
        -L ${tumour_merged_vcf} \
        --sequence-dictionary ${genome_dict} \
        -O ${tumour_sample}.forcecall.vcf.gz

    bcftools norm -f ${genome} -m-any ${tumour_sample}.forcecall.vcf.gz -o ${tumour_sample}.forcecall.norm.vcf.gz -O z
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/tumour_samples/${tumour_sample}/mutect2/* ./
    """
}
