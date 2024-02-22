process mergefastq {
    input:
    tuple val(patient), val(meta), path(fastq_r1_paths)
    tuple val(patient), val(meta), path(fastq_r2_paths)

    publishDir "${params.outputdir}/${patient}/mergedfastq/", mode: 'copy'

    output:
    tuple val(patient), val(meta), path("$patient.${meta['sample']}.merged.r1.fastq.gz"), path("$patient.${meta['sample']}.merged.r2.fastq.gz"), emit: merged_fastq
    
    script:
    """
    cat ${fastq_r1_paths.join(" ")} > ${patient}.${meta['sample']}.merged.r1.fastq.gz
    cat ${fastq_r2_paths.join(" ")} > ${patient}.${meta['sample']}.merged.r2.fastq.gz
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/PEA020/mergedfastq/* ./


    #touch ${patient}.${meta['sample']}.merged.r1.fastq.gz
    #touch ${patient}.${meta['sample']}.merged.r2.fastq.gz
    """
}
