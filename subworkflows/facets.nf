include { FACETS_SNP_PILEUP } from '../modules/snppileup.nf'
include { FACETS_RUN } from '../modules/snppileup.nf'

workflow FACETS_CN_CALLING {
    take:
    ch_tumour_crams
    ch_normal_crams
    ch_normal_vcf

    main:
    //[patient, normal_sample_id, normal_sample_cram_path, tumour_sample_id, tumour_sample_cram_path, normal_haplotypecaller_vcf]
    facets_cram_input = ch_normal_crams
        .cross(ch_tumour_crams)
        .map { it.flatten().unique() }
        .combine(ch_normal_vcf, by: 0)
        .map { it.flatten().unique() }
        .dump(tag: 'facets_cram_input', pretty: true)

    FACETS_SNP_PILEUP(facets_cram_input)
    FACETS_RUN(FACETS_SNP_PILEUP.out.pileup_file.dump(tag: 'FACETS_SNP_PILEUP.out.pileup_file', pretty: true))

    emit:
    facets_output = FACETS_RUN.out.facets_output
}
