process PREPARE_PYCLONE_INPUT {
    input:
    tuple val(patient), val(sample), path(mutect2_vcf), path(facets_output)

    publishDir "${params.outputdir}/${patient}/pyclone/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(sample), path("${patient}.${sample}.pyclone.input.csv"), emit: pyclone_input_per_sample

    script:
    """
    /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/bin/parse.mut.cn.for.pyclone.R \
        ${patient} \
        ${sample} \
        ${mutect2_vcf} \
        ${facets_output}
    """

    stub:
    """
    ln -sf ${params.stub_data_dir}/${patient}/pyclone/* ./
    """
}

process MERGE_PYCLONE_INPUT_PER_PATIENT {
    input:
    tuple val(patient), val(samples), path(csvs_to_merge)

    publishDir "${params.outputdir}/${patient}/pyclone/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.pyclone.input.merged.tsv"), emit: merged_pyclone_input

    script:
    """
    head -n 1 ${csvs_to_merge[0]} > ${patient}.pyclone.input.csv
    tail -n +2 -q ${csvs_to_merge} >> ${patient}.pyclone.input.csv
    #convert to TSV
    sed 's/,/\\t/g' ${patient}.pyclone.input.csv > ${patient}.pyclone.input.merged.tsv
    """

    stub:
    """
    ln -sf ${params.stub_data_dir}/${patient}/pyclone/* ./
    """
}

process RUN_PYCLONE {
    input:
    tuple val(patient), path(pyclone_input)

    publishDir "${params.outputdir}/${patient}/pyclone/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.pyclone.results.tsv"), emit: pyclone_results

    script:
    """
    pyclone-vi fit -i ${pyclone_input} -o ${patient}.pyclone.h5 -c 40 -d beta-binomial -r 10
    pyclone-vi write-results-file -i ${patient}.pyclone.h5 -o ${patient}.pyclone.results.tsv
    """

    stub:
    """
    ln -sf ${params.stub_data_dir}/${patient}/pyclone/* ./
    """
}
