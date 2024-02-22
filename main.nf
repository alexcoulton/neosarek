//////////////////////////////
//PARAMETERS
//////////////////////////////
nextflow.enable.dsl = 2

// mela sample sheet
//params.samplesheet = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/mela.samplesheet.csv'
// sarek sample sheet
params.samplesheet = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/sarek/newsamplesheetpea020.csv'
params.sarek_output_dir = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/sarek/output_comb'

params.outputdir = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/output/'
params.stub_data_dir = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/stub_data/'

params.hlahd_install_directory = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/singularity/hlahd.1.7.0'

params.vepcache = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/vep_cache/'
params.vepplugins = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/singularity/pvac_vep_plugins'

params.genome_reference = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/igenome/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'
params.genome_index = params.genome_reference + '.fai'
params.genome_dict = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/igenome/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict'

params.germline_resource = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/igenome/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz'
params.germline_resource_index = params.germline_resource + '.tbi'

params.panel_of_normals = '/nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/igenome/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz'
params.panel_of_normals_index = params.panel_of_normals + '.tbi'


//////////////////////////////
//IMPORT PROCESSES
//////////////////////////////

include { mergefastq } from './modules/merge.fastq.nf'
include { mergevcf } from './modules/mergevcf.nf'
include { vepanno } from './modules/vepanno.nf'
include { hlahd } from './modules/hlahd.nf'
include { mutect2_forcecall } from './modules/mutect2.nf'
include { pvacseq } from './modules/pvacseq.nf'

//////////////////////////////
//FUNCTIONS
//////////////////////////////

// Function to process normal samples
def processNormalSamples(normalSamples, int readIndex) {
    return normalSamples
        .map { [it[0], it[1], it[readIndex]] }
        .collect()
        .map { it.unique() }
        .map { [it[0], it[1], it.subList(2, it.size())] } //account for varying numbers of lanes
}

//////////////////////////////
//WORKFLOW
//////////////////////////////

workflow {
    sample_channel = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> [row.patient, [sample: row.sample, status: row.status], file(row.fastq_1), file(row.fastq_2)] }

    normal_samples = sample_channel.filter{ it[1].status == '0' }
    tumour_samples = sample_channel.filter{ it[1].status == '1' }

    //Sarek outputs tumour VCFs in a tumour_vs_normal path
    //Here we join by patient and get the mutect2 output VCFs

    //won't need the fastq columns for the tumour samples so I scrub them here
    tumour_samples = tumour_samples
        .map { [it[0], it[1].sample] }
        .unique()

    normal_samples_trimmed = normal_samples //get rid of lane / fastq information
        .map { [it[0], it[1].sample] }
        .unique()

    //group tumour VCFs by patient for mutect2 forced calling
    tumour_samples_vcf = tumour_samples
        .join(
            normal_samples
                .map { [it[0], it[1].sample] }
        )
        .map { [patient: it[0], tumour_samp: it[1], normal_samp: it[2]]  }
        .map { [it.patient, it.tumour_samp, it.normal_samp, params.sarek_output_dir + '/variant_calling/mutect2/' + it.tumour_samp + '_vs_' + it.normal_samp + '/' + it.tumour_samp + '_vs_' + it.normal_samp + '.mutect2.filtered.vcf.gz', params.sarek_output_dir + '/variant_calling/mutect2/' + it.tumour_samp + '_vs_' + it.normal_samp + '/' + it.tumour_samp + '_vs_' + it.normal_samp + '.mutect2.filtered.vcf.gz.tbi']}
        .groupTuple() //groupTuple throws error when list elements are named

    mergevcf(
        tumour_samples_vcf,
        params.genome_reference,
        params.genome_index,
        params.genome_dict
    )

    tumour_crams = tumour_samples
        .map { [it[0], it[1], params.sarek_output_dir + '/preprocessing/recalibrated/' + it[1] + '/' + it[1] + '.recal.cram'] }

    normal_crams = normal_samples_trimmed
        .map { [it[0], it[1], params.sarek_output_dir + '/preprocessing/recalibrated/' + it[1] + '/' + it[1] + '.recal.cram'] }

    ////for every tumour .cram, create a list containing the merged tumour VCF and the normal .cram
    mutect2_files = normal_crams
        .cross(
            mergevcf.out.merged_vcf
                .cross(tumour_crams)
                .map { it.flatten().unique() }
        )
        .map { it.flatten().unique() }

    //mutect2_forcecall(
        //mutect2_files,
        //params.germline_resource,
        //params.germline_resource_index,
        //params.panel_of_normals,
        //params.panel_of_normals_index,
        //params.genome_reference,
        //params.genome_index,
        //params.genome_dict
    //)

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
    vepanno(params.vepcache, params.genome_reference, params.vepplugins, mergevcf.out.merged_vcf)

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
