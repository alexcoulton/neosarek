include { PREPARE_PYCLONE_INPUT } from '../modules/pyclone.nf'
include { MERGE_PYCLONE_INPUT_PER_PATIENT } from '../modules/pyclone.nf'
include { RUN_PYCLONE } from '../modules/pyclone.nf'

workflow PYCLONE {
    take:
    ch_mutect2_out
    ch_facets_out

    main:
    //[patient, sample, mutect2_vcf, facets_rds]
    pyclone_input = ch_mutect2_out
        .map { ["${it[0]}_${it[1]}", it[0], it[1], it[2]] }
        .combine(
            ch_facets_out
                .map { ["${it[0]}_${it[1]}", it[0], it[1], it[2], it[3]] }, by: 0
            )
        .map{ it.unique() }
        .map { [it[1], it[2], it[3], it[4]]}

    PREPARE_PYCLONE_INPUT(pyclone_input)

    MERGE_PYCLONE_INPUT_PER_PATIENT(
        PREPARE_PYCLONE_INPUT.out.pyclone_input_per_sample
            .groupTuple()
    )

    RUN_PYCLONE(MERGE_PYCLONE_INPUT_PER_PATIENT.out.merged_pyclone_input)
}
