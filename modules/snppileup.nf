process FACETS_SNP_PILEUP {
    input:
    tuple val(patient), val(normal_sample), path(normal_cram), val(tumour_sample), path(tumour_cram), path(normal_vcf)

    publishDir "${params.outputdir}/${patient}/facets/snppileup/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(normal_sample), val(tumour_sample), path("${tumour_sample}.facets.pileup.csv"), emit: pileup_file

    script:
    """
    /camp/home/coultoa/work/external.repos/facets/inst/extcode/snp-pileup ${normal_vcf} ./${tumour_sample}.facets.pileup.csv ${normal_cram} ${tumour_cram}
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/facets/snppileup/${tumour_sample}.facets.pileup.csv ./
    """
}

process FACETS_RUN {
    input:
    tuple val(patient), val(normal_sample), val(tumour_sample), path(snp_pileup_file)

    publishDir "${params.outputdir}/${patient}/facets/output/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(tumour_sample), path("${tumour_sample}.facets.rds"), path("${tumour_sample}.facets.plot.pdf"), path("${tumour_sample}.facets.diag.plot.pdf"), emit: facets_output

    script:
    """
    ~/work/cpi.nextflow/bin/run.facets.R ${tumour_sample} ${snp_pileup_file} ${tumour_sample}.facets.rds ${tumour_sample}.facets.plot.pdf ${tumour_sample}.facets.diag.plot.pdf
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/facets/output/${tumour_sample}* ./
    """
}

