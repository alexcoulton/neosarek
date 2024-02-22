process pvacseq {
    input:
    tuple val(patient), val(normal_sample), val(hlatype), val(tumour_samples), path(vepvcf)

    output:
    path "./output"

    script:
    """
    grep "#CHROM" ${vepvcf}

    pvacseq run \
        ${vepvcf} \
        ${patient}_${sample} \
        ${hlatype} \
        NetMHC \
        ./output \
        -e1 8,9,10 \
        -e2 15 \
        --normal-sample-name PEA020_SPA659A22_normal
    """
}
