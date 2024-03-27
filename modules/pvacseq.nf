process PVACSEQ {
    input:
    tuple val(patient), val(normal_sample), val(hlatype), val(tumour_sample), path(vepvcf)

    publishDir "${params.outputdir}/${patient}/pvacseq/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite


    output:
    path "./output"

    script:
    """
    grep "#CHROM" ${vepvcf}

    pvacseq run \
        ${vepvcf} \
        ${patient}_${tumour_sample} \
        ${hlatype} \
        NetMHC \
        ./output \
        -e1 8,9,10 \
        -e2 15 \
        -t 8 \
        --normal-sample-name ${patient}_${normal_sample} \
        --iedb-install-directory /opt/iedb
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/pvacseq/* ./
    """
}
