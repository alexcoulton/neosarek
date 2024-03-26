process MERGE_FASTQ {
    input:
    tuple val(patient), val(sample), path(fastq_r1_paths)
    tuple val(patient), val(sample), path(fastq_r2_paths)

    publishDir "${params.outputdir}/${patient}/mergedfastq/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(sample), path("$patient.${sample}.merged.r1.fastq.gz"), path("$patient.${sample}.merged.r2.fastq.gz"), emit: merged_fastq
    
    script:
    """
    cat ${fastq_r1_paths.join(" ")} > ${patient}.${sample}.merged.r1.fastq.gz
    cat ${fastq_r2_paths.join(" ")} > ${patient}.${sample}.merged.r2.fastq.gz
    """

    stub:
    """
    ln -s ${params.stub_data_dir}/${patient}/mergedfastq/* ./
    """
}
