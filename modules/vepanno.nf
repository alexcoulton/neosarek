process VEP_ANNO {
    input:
    path vepcache
    path vepfasta
    path vepplugins
    tuple val(patient), val(tumour_samples), path(merged_vcf), path(merged_vcf_index)

    publishDir "${params.outputdir}/${patient}/vep_annotated_merged_vcf/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite


    output:
    tuple val(patient), val(tumour_samples), path("${patient}.vep.vcf"), emit: vep_vcf

    script:
    """
    vep \
        --input_file ${merged_vcf} \
        --output_file ${patient}.vep.vcf \
        --format vcf \
        --vcf \
        --symbol \
        --terms SO \
        --tsl \
        --biotype \
        --hgvs \
        --fasta ${vepfasta} \
        --offline \
        --cache ${vepcache} \
        --plugin Frameshift \
        --plugin Wildtype \
        --dir_plugins ${vepplugins}
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/vep_annotated_merged_vcf/${patient}.vep.vcf ./
    """
}
