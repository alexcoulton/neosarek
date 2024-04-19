nextflow.enable.dsl = 2

//////////////////////////////
//STAGE PARAMETER
//////////////////////////////

validStages = ['initial', 'review', 'post_review', 'facets_only']
if(!validStages.contains(params.stage)){
    error "Invalid parameter stage. Possible values are ${validStages.join(', ')}"
}

//////////////////////////////
//IMPORT SUBWORKFLOWS
//////////////////////////////

include { FACETS_CN_CALLING } from './subworkflows/facets.nf'
include { HLA_NEOANTIGEN } from './subworkflows/hla.neoantigen.nf'
include { PYCLONE } from './subworkflows/pyclone.nf'

//////////////////////////////
//IMPORT PROCESSES
//////////////////////////////

include { NORMALIZE_FILTER_VCF } from './modules/mergevcf.nf'
include { MERGE_VCF } from './modules/mergevcf.nf'
include { MUTECT2_FORCECALL } from './modules/mutect2.nf'

//////////////////////////////
//WORKFLOW
//////////////////////////////

workflow {
    //////////////////////////////
    //CHANNEL INITIATION
    //////////////////////////////
    sample_channel = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> [row.patient, [sample: row.sample, status: row.status, ffpe: row.ffpe], file(row.fastq_1), file(row.fastq_2)] }
        .dump(tag: 'sample_channel')

    //[patient, [sample, status, ffpe], fastq_1, fastq_2]
    normal_samples = sample_channel.filter{ it[1].status == '0' }
        .dump(tag: 'normal_samples')

    //[patient, [sample, status, ffpe], fastq_1, fastq_2]
    tumour_samples = sample_channel.filter{ it[1].status == '1' }

    //[patient, sample]
    tumour_samples_trimmed = tumour_samples //remove lane / fastq information
        .map { [it[0], it[1].sample, it[1].ffpe] }
        .unique()
        .dump(tag: 'tumour_samples')

    //[patient, sample]
    normal_samples_trimmed = normal_samples //get rid of lane / fastq information
        .map { [it[0], it[1].sample] }
        .unique()
        .dump(tag: 'normal_samples_trimmed')

    //[patient, tumour_sample, tumour_cram]
    tumour_crams = tumour_samples_trimmed
        .map { [it[0], it[1], params.sarek_output_dir + '/preprocessing/recalibrated/' + it[1] + '/' + it[1] + '.recal.cram'] }
        .dump(tag: 'tumour_crams')

    //[patient, normal_sample, normal_cram]
    normal_crams = normal_samples_trimmed
        .map { [it[0], it[1], params.sarek_output_dir + '/preprocessing/recalibrated/' + it[1] + '/' + it[1] + '.recal.cram'] }
        .dump(tag: 'normal_crams')

    //[patient, normal_sample, normal_vcf]
    normal_vcf = normal_samples_trimmed
        .map { [it[0], it[1], params.sarek_output_dir + '/variant_calling/haplotypecaller/' + it[1] + '/' + it[1] + '.haplotypecaller.filtered.vcf.gz'] }
        .dump(tag: 'normal_vcf')


    //[patient, [tumour_samples], [non_ffpe_tumour_samples], [normal_sample(s)], [tumour_sample_mut2_vcfs], [tumour_sample_mut2_vcf_indexes]]
    tumour_samples_vcf = tumour_samples_trimmed
        .combine(normal_samples_trimmed, by: 0)
        .map { [patient: it[0], tumour_samp: it[1], ffpe: it[2], normal_samp: it[3]]  }
        .map { [
            it.patient,
            it.tumour_samp,
            it.ffpe,
            it.normal_samp,
            params.sarek_output_dir + '/variant_calling/mutect2/' + it.tumour_samp + '_vs_' + it.normal_samp + '/' + it.tumour_samp + '_vs_' + it.normal_samp + '.mutect2.filtered.vcf.gz',
            params.sarek_output_dir + '/variant_calling/mutect2/' + it.tumour_samp + '_vs_' + it.normal_samp + '/' + it.tumour_samp + '_vs_' + it.normal_samp + '.mutect2.filtered.vcf.gz.tbi'
            ]}

    //////////////////////////////
    //EXECUTION STAGES
    //////////////////////////////

    if(params.stage == 'initial') {
        NORMALIZE_FILTER_VCF(
            tumour_samples_vcf,
            params.genome_reference,
            params.genome_index,
            params.genome_dict
        )

        //[patient, tumour_samples, ffpe_status, normal_samples, vcfs, tbis]
        merge_vcf_input = NORMALIZE_FILTER_VCF.out.normalized_filtered_vcf
            .filter { it[2] != '1' } //filter out FFPE samples
            .groupTuple()

        MERGE_VCF(
            merge_vcf_input,
            params.genome_reference,
            params.genome_index,
            params.genome_dict
        )

        //HLA_NEOANTIGEN(normal_samples, MERGE_VCF.out.merged_vcf.dump(tag: 'MERGE_VCF.out.merged_vcf'))
        //FACETS_CN_CALLING(tumour_crams, normal_crams, normal_vcf)

        ////for every tumour .cram, create a list containing the merged tumour VCF and the normal .cram
        //[patient, normal_sample, normal_cram, master_tumour_vcf, master_tumour_vcf_index, tumour_sample, tumour_cram]
        mutect2_files = normal_crams
            .cross(
                MERGE_VCF.out.merged_vcf
                    .map { [it[0], it[2], it[3]] }
                    .cross(tumour_crams)
                    .map { it.flatten().unique() }
            )
            .map { it.flatten().unique() }
            .map { [
                it[0],
                it[1],
                it[2],
                it[3],
                it[4],
                it[5],
                it[6],
                params.sarek_output_dir + '/variant_calling/mutect2/' + it[5] + '_vs_' + it[1] + '/' + it[5] + '_vs_' + it[1] + '.mutect2.filtered.vcf.gz',
                params.sarek_output_dir + '/variant_calling/mutect2/' + it[5] + '_vs_' + it[1] + '/' + it[5] + '_vs_' + it[1] + '.mutect2.filtered.vcf.gz.tbi'
                ] }
            .dump(tag: 'mutect2_files')

        MUTECT2_FORCECALL(
            mutect2_files,
            params.germline_resource,
            params.germline_resource_index,
            params.panel_of_normals,
            params.panel_of_normals_index,
            params.genome_reference,
            params.genome_index,
            params.genome_dict
        )
    }

    if(params.stage == 'facets_only'){
        FACETS_CN_CALLING(tumour_crams, normal_crams, normal_vcf)
    }

    if(params.stage == 'review'){
        //[patient, sample, mutect2_vcf, reviewed_facets_rds]
        tumour_samples_pyclone = tumour_samples
            .map { [it[0], it[1].sample, it[1].ffpe] }
            .unique()
            .dump(tag: 'tumour_samples')
            .map {
                    [
                        it[0],
                        it[1],
                        params.outputdir + '/' + it[0] + '/tumour_samples/' + it[1] + '/mutect2/' + it[1] + '.forcecall.norm.isec.vcf.gz',
                        params.outputdir + '/' + it[0] + '/facets/review/' + it[1] + '.facets.rds', //reviewed solutions
                    ]
                }

        PYCLONE(tumour_samples_pyclone)
    }
}
