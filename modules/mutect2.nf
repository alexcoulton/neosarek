process MUTECT2_FORCECALL {
    input:
    tuple val(patient), val(normal_sample), path(normal_cram), path(tumour_merged_vcf), path(tumour_merged_tbi), val(tumour_sample), path(tumour_cram), path(orig_mutect2_vcf), path(orig_mutect2_tbi)
    path germline_resource
    path germline_resource_index
    path panel_of_normals
    path panel_of_normals_index
    path genome
    path genome_index
    path genome_dict

    publishDir "${params.outputdir}/${patient}/tumour_samples/${tumour_sample}/mutect2/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(tumour_sample), path("${tumour_sample}.forcecall.roq.norm.isec.rescue_annot.vcf.gz"), emit: mutect2_vcf
    
    script:
    """
    gatk Mutect2 \
        -R ${genome} \
        -I ${tumour_cram} \
        -I ${normal_cram} \
        -normal ${patient}_${normal_sample} \
        --germline-resource ${germline_resource} \
        --panel-of-normals ${panel_of_normals} \
        --f1r2-tar-gz ${tumour_sample}.f1r2.tar.gz \
        -alleles ${tumour_merged_vcf} \
        --force-call-filtered-alleles \
        -L ${tumour_merged_vcf} \
        --sequence-dictionary ${genome_dict} \
        -O ${tumour_sample}.forcecall.vcf.gz

    gatk LearnReadOrientationModel \
        -I ${tumour_sample}.f1r2.tar.gz \
        -O ${tumour_sample}.readorientation.model.tar.gz

    gatk FilterMutectCalls \
        -R ${genome} \
        -V ${tumour_sample}.forcecall.vcf.gz \
        --ob-priors ${tumour_sample}.readorientation.model.tar.gz \
        -O ${tumour_sample}.forcecall.roq.vcf.gz


    bcftools norm -f ${genome} -m-any ${tumour_sample}.forcecall.roq.vcf.gz -o ${tumour_sample}.forcecall.roq.norm.vcf.gz -O z
    bcftools index ${tumour_sample}.forcecall.roq.norm.vcf.gz

    #intersect the force-called vcf with the master vcf to make sure no new variants are called
    bcftools isec -n=2 -w1 ${tumour_sample}.forcecall.roq.norm.vcf.gz ${tumour_merged_vcf} -Oz -o ${tumour_sample}.forcecall.roq.norm.isec.vcf.gz

    tabix ${tumour_sample}.forcecall.roq.norm.isec.vcf.gz

    bcftools annotate \
        -a ${orig_mutect2_vcf} \
        -m DENOVO_CALLED \
        ${tumour_sample}.forcecall.roq.norm.isec.vcf.gz \
        -o ${tumour_sample}.forcecall.roq.norm.isec.rescue_annot.vcf.gz 
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/tumour_samples/${tumour_sample}/mutect2/* ./
    """
}
