include { PREPARE_PYCLONE_INPUT } from '../modules/pyclone.nf'
include { MERGE_PYCLONE_INPUT_PER_PATIENT } from '../modules/pyclone.nf'
include { MERGE_CONIPHER_INPUT_PER_PATIENT } from '../modules/pyclone.nf'
include { RUN_PYCLONE } from '../modules/pyclone.nf'
include { RUN_PAIRTREE } from '../modules/pyclone.nf'
include { PREPARE_PAIRTREE_INPUT } from '../modules/pyclone.nf'
include { RUN_CONIPHER } from '../modules/pyclone.nf'

workflow PYCLONE {
    take:
    ch_pyclone_input

    main:
    //[patient, sample, mutect2_vcf, facets_rds]
    //pyclone_input = ch_mutect2_out
        //.map { ["${it[0]}_${it[1]}", it[0], it[1], it[2]] }
        //.combine(
            //ch_facets_out
                //.map { ["${it[0]}_${it[1]}", it[0], it[1], it[2], it[3]] }, by: 0
            //)
        //.map{ it.unique() }
        //.map { [it[1], it[2], it[3], it[4]]}

    PREPARE_PYCLONE_INPUT(ch_pyclone_input, params.conipher_prefix)

    MERGE_PYCLONE_INPUT_PER_PATIENT(
        PREPARE_PYCLONE_INPUT.out.pyclone_input_per_sample
            .groupTuple()
    )

    RUN_PYCLONE(MERGE_PYCLONE_INPUT_PER_PATIENT.out.merged_pyclone_input)

    pairtree_input_prep = ch_pyclone_input
        .groupTuple()

    pairtree_input_prep = pairtree_input_prep
        .combine(RUN_PYCLONE.out.pyclone_results, by: 0)

    PREPARE_PAIRTREE_INPUT(pairtree_input_prep)
    //RUN_PAIRTREE(PREPARE_PAIRTREE_INPUT.out.pairtree_input_files)

    //MERGE_CONIPHER_INPUT_PER_PATIENT(
        //PREPARE_PYCLONE_INPUT.out.conipher_input_per_sample
            //.groupTuple()
    //)

    //RUN_CONIPHER(
        //MERGE_CONIPHER_INPUT_PER_PATIENT.out.merged_conipher_input,
        //params.conipher_prefix,
        //'/nemo/project/proj-tracerX/working/CMELA/alex/work/external.repos/CONIPHER'
    //)
}
