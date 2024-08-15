process PREPARE_PYCLONE_INPUT {
    input:
    tuple val(patient), val(sample), path(mutect2_vcf), path(facets_output), path(annotations)
    val(conipher_prefix)

    publishDir "${params.outputdir}/${patient}/pyclone/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), val(sample), env(PURITY), path("${patient}.${sample}.pyclone.input.tsv"), emit: pyclone_input_per_sample
    tuple val(patient), val(sample), path("${patient}.${sample}.conipher.input.tsv"), emit: conipher_input_per_sample

    script:
    """
    echo ""
    echo ""
    echo ""
    echo ""
    echo ""
    echo ""
    /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/bin/parse.mut.cn.for.pyclone.R \
        ${patient} \
        ${sample} \
        ${mutect2_vcf} \
        ${facets_output} \
        ${conipher_prefix}

    PURITY=`cat "${patient}.${sample}.purity.txt"`
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

process RUN_PYCLONEVI {
    input:
    tuple val(patient), path(pyclone_input)

    publishDir "${params.outputdir}/${patient}/pyclonevi/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.pyclonevi.results.tsv"), emit: pyclonevi_results

    script:
    """
    echo "Run pyclone-vi fit"
    pyclone-vi
    pyclone-vi fit -i ${pyclone_input} -o ${patient}.pyclone.h5 -c 40 -d beta-binomial -r 10
    echo "Run pyclone-vi write-results-file"
    pyclone-vi write-results-file -i ${patient}.pyclone.h5 -o ${patient}.pyclone.results.tsv
    """

    stub:
    """
    ln -sf ${params.stub_data_dir}/${patient}/pyclonevi/* ./
    """
}

process RUN_PYCLONE {
    debug true
    input:
    tuple val(patient), val(samples), val(purities), path(sample_tsvs)

    publishDir "${params.outputdir}/${patient}/pyclone/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.pyclone.results.tsv"), emit: pyclonevi_results

    script:
    def purity_list = purities.join(' ') //convert to space delimited list
    def sample_list = samples.join(' ')
    """
    #mkdir ./numba_cache
    #export NUMBA_CACHE_DIR=./numba_cache

    echo "PyClone binary:"
    echo \$(which PyClone)

    echo "PyClone setup_analysis"
    PyClone setup_analysis \
        --working_dir ./pyclone_output \
        --in_files ${sample_tsvs} \
        --tumour_contents ${purity_list} \
        --samples ${sample_list}



    echo "PyClone run_analysis"
    /nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/envs/pyclone_custom2/bin/python \
        -m trace \
        --trace /nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/envs/pyclone_custom2/bin/PyClone \
        run_analysis \
        --config_file ./pyclone_output/config.yaml




    #PyClone run_analysis \
    #    --config_file ./pyclone_output/config.yaml


    #PyClone run_analysis_pipeline \
    #    --working_dir ./pyclone_output \
    #    --in_files ${sample_tsvs} \
    #    --tumour_contents ${purity_list} \
    #    --samples ${sample_list} &> pyclone.log
    """

    stub:
    """
    ln -sf ${params.stub_data_dir}/${patient}/pyclonevi/* ./
    """
}

process COLLATE_MUTATION_CN_DATA {
    input:
    tuple val(patient), val(sample), path(mutect_vcfs), path(facets_rds), path(annotations), path(pyclone_results)

    publishDir "${params.outputdir}/${patient}/mutation_cn/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.combined.mutations.tsv"), path("${patient}.combined.seg.tsv"), emit: combined_mutation_seg_data

    script:
    """
    /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/bin/collate.mut.cn.R \
        ${patient} \
        ${sample} \
        ${mutect_vcfs} \
        ${facets_rds} \
        ${pyclone_results} \
        ${annotations}
    """
}

process PREPARE_PAIRTREE_INPUT {
    input:
    tuple val(patient), path(mut), path(cn)

    publishDir "${params.outputdir}/${patient}/pairtree/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite

    output:
    tuple val(patient), path("${patient}.pairtree.ssm"), path("${patient}.pairtree.input.json"), path("${patient}.combined.mutations.tsv"), emit: pairtree_input_files

    script:
    """
    echo ""
    /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/bin/pairtree.input.prep.R \
        ${patient} \
        ${mut}
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
