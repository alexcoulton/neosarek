//////////////////////////////
//FUNCTIONS
//////////////////////////////

// Function to process normal samples
def processNormalSamples(normalSamples, int readIndex) {
    return normalSamples
        .map { [it[0], it[1].sample, it[readIndex]] }
        .groupTuple()
        .map { [it[0], it[1].unique()[0], it[2]]}
}

include { MERGE_FASTQ } from '../modules/merge.fastq.nf'
include { HLAHD } from '../modules/hlahd.nf'
include { VEP_ANNO } from '../modules/vepanno.nf'
include { PVACSEQ } from '../modules/pvacseq.nf'

workflow HLA_NEOANTIGEN {
    take:
    ch_normal_samples
    ch_merged_vcf

    main:
    // should only have one normal sample in the sample sheet per patient.
    // need to add the necessary validation to ensure this is the case
    // OR - specify one GERMLINE sample per patient. have a distinction 
    // between normal and GERMLINE samples in the pipeline samplesheet

    //for now I will assume that normal == GERMLINE

    // Processing normal samples for read1 and read2 lanes

    //[patient, [normal_sample, status], fastq_1, fastq_2]
    normal_read1_lanes = processNormalSamples(ch_normal_samples, 2) // for read1
        .dump(tag: 'normal_read1_lanes', pretty: true)
    //[patient, [normal_sample, status], fastq_1, fastq_2]
    normal_read2_lanes = processNormalSamples(ch_normal_samples, 3) // for read2
        .dump(tag: 'normal_read2_lanes', pretty: true)


    MERGE_FASTQ(normal_read1_lanes, normal_read2_lanes)
    HLAHD(MERGE_FASTQ.out.merged_fastq.dump(tag: 'MERGE_FASTQ.out.merged_fastq', pretty: true))

    //Sarek VEP does not have the required pVACseq plugins
    //Need to run our own VEP annotation process
    VEP_ANNO(params.vepcache, params.genome_reference, params.vepplugins, ch_merged_vcf)

    //Here we cross the hlahd results from the normal samples
    //with the VEP-annotated tumour sample VCFs
    //The cross is performed using the patient ID as the key

    PVACSEQ(
        HLAHD.out.hlaresult
            .cross(
                VEP_ANNO.out.vep_vcf
                    .map { [it[0], it[1][0], it[2]] }
            )
            .map{ it.flatten() }
            .map{ it.unique() }
            .dump(tag: 'HLAHD.out.hlaresult.cross', pretty: true)
    )
}
