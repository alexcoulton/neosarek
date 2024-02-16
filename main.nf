nextflow.enable.dsl = 2
params.samplesheet = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/mela.samplesheet.csv'
params.hlahd_install_directory = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/singularity/hlahd.1.7.0'

params.vepcache = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/sarek/work/stage-196db640-7c93-4a0d-8583-26529571c7d8/bd/911b8f40e0b07793cc188da2b38c92/110_GRCh38'
params.vepfasta = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/sarek/work/stage-22201135-817d-4102-b5f9-4718d011a38a/33/4a7f91b4d1de78398e241f657a139d/Homo_sapiens_assembly38.fasta'
params.vepplugins = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/singularity/pvac_vep_plugins'

include { hlahd } from './modules/hlahd.nf'
include { pvacseq } from './modules/pvacseq.nf'
include { vepanno } from './modules/vepanno.nf'
include { mergefastq } from './modules/merge.fastq.nf'

// Function to process normal samples
def processNormalSamples(normalSamples, int readIndex) {
    return normalSamples
        .map { [it[0], it[1], it[readIndex]] }
        .collect()
        .map { it.unique() }
        .map { [it[0], it[1], it.subList(2, it.size())] } //account for varying numbers of lanes
}


workflow {
    sample_channel = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> [row.patient, [sample: row.sample, status: row.status], file(row.merged_fastq_1), file(row.merged_fastq_2), file(row.variants_vcf)] }

    normal_samples = sample_channel.filter{ it[1].status == '0'}
    tumour_samples = sample_channel.filter{ it[1].status == '1'}

    // should only have one normal sample in the sample sheet per patient.
    // need to add the necessary validation to ensure this is the case
    // OR - specify one GERMLINE sample per patient. have a distinction 
    // between normal and GERMLINE samples in the pipeline samplesheet

    //for now I will assume that normal == GERMLINE

    // Processing normal samples for read1 and read2 lanes
    normal_read1_lanes = processNormalSamples(normal_samples, 2) // for read1
    normal_read2_lanes = processNormalSamples(normal_samples, 3) // for read2

    mergefastq(normal_read1_lanes, normal_read2_lanes)
    hlahd(mergefastq.out.merged_fastq)

    //Sarek VEP does not have the required pVACseq plugins
    //Need to run our own VEP annotation process
    //vepanno(params.vepcache, params.vepfasta, params.vepplugins, tumour_samples)

    //Here we cross the hlahd results from the normal samples
    //with the VEP-annotated tumour sample VCFs
    //The cross is performed using the patient ID as the key
    //pvacseq(hlahd.out.hlaresult.cross(vepanno.out.vep_vcf))
    //pvacseq(
        //hlahd.out.hlaresult
            //.cross(vepanno.out.vep_vcf)
            //.map{ it.flatten() }
            //.map{ it.unique() }
            //.view()
    //)
}
