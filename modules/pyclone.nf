process PREPARE_PYCLONE_INPUT {
    input:
    tuple val(patient), val(sample), path(mutect2_vcf), path(facets_output)
    val(conipher_prefix)

    publishDir "${params.outputdir}/${patient}/pyclone/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(sample), path("${patient}.${sample}.pyclone.input.tsv"), emit: pyclone_input_per_sample
    tuple val(patient), val(sample), path("${patient}.${sample}.conipher.input.tsv"), emit: conipher_input_per_sample

    script:
    """
    /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/bin/parse.mut.cn.for.pyclone.R \
        ${patient} \
        ${sample} \
        ${mutect2_vcf} \
        ${facets_output} \
        ${conipher_prefix}
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
    head -n 1 ${csvs_to_merge[0]} > ${patient}.pyclone.input.merged.tsv
    tail -n +2 -q ${csvs_to_merge} >> ${patient}.pyclone.input.merged.tsv
    #convert to TSV
    #sed 's/,/\\t/g' ${patient}.pyclone.input.csv > ${patient}.pyclone.input.merged.tsv
    """
}

process MERGE_CONIPHER_INPUT_PER_PATIENT {
    input:
    tuple val(patient), val(samples), path(csvs_to_merge)

    publishDir "${params.outputdir}/${patient}/conipher/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.conipher.input.merged.tsv"), emit: merged_conipher_input

    script:
    """
    head -n 1 ${csvs_to_merge[0]} > ${patient}.conipher.input.merged.tsv
    tail -n +2 -q ${csvs_to_merge} >> ${patient}.conipher.input.merged.tsv
    #convert to TSV
    #sed 's/,/\\t/g' ${patient}.conipher.input.csv > ${patient}.conipher.input.merged.tsv
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

process PREPARE_PAIRTREE_INPUT {
    input:
    tuple val(patient), val(sample), path(mutect_vcfs), path(facets_rds), path(pyclone_results)

    publishDir "${params.outputdir}/${patient}/pairtree/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.pairtree.ssm"), path("${patient}.pairtree.input.json"), path("${patient}.combined.mutations.tsv"), emit: pairtree_input_files

    script:
    """
    echo ""
    /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/bin/pairtree.input.prep.R \
        ${patient} \
        ${sample} \
        ${mutect_vcfs} \
        ${facets_rds} \
        ${pyclone_results}
    """
}

process RUN_PAIRTREE {
    input:
    tuple val(patient), path(pairtree_ssm), path(pairtree_json)

    publishDir "${params.outputdir}/${patient}/pairtree/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple path("${patient}.pairtree.results.html"), path("${patient}.pairtree.results.npz"), emit: pairtree_results

    script:
    """
    pairtree \
        --params ${pairtree_json} \
        ${pairtree_ssm} \
        ${patient}.pairtree.results.npz

    plottree \
        --runid ${patient} \
        ${pairtree_ssm} \
        ${pairtree_json} \
        ${patient}.pairtree.results.npz \
        ${patient}.pairtree.results.html
    """
}

process RUN_CONIPHER {
    input:
    tuple val(patient), path(conipher_input)
    val(conipher_prefix)
    path(CONIPHER_DIR)

    publishDir "${params.outputdir}/${patient}/conipher/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}_conipher_results"), emit: pyclone_results

    script:
    """
    mkdir ./numba_cache
    export NUMBA_CACHE_DIR=./numba_cache
    /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/bin/run.conipher.R \
        ${patient} \
        ${conipher_prefix} \
        ${conipher_input}
    """
}
