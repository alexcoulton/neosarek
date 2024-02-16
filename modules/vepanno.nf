process vepanno {
    input:
    path vepcache
    path vepfasta
    path vepplugins
    tuple val(patient), val(meta), path(merged_fastq_1), path(merged_fastq_2), path(variants_vcf)

    output:
    tuple val(patient), val("${meta['sample']}"), path("${patient}.${meta['sample']}.vep.vcf"), emit: vep_vcf

    script:
    """
    vep \
        --input_file ${variants_vcf} \
        --output_file ${patient}.${meta['sample']}.vep.vcf \
        --format vcf \
        --vcf \
        --symbol \
        --terms SO \
        --tsl \
        --biotype \
        --hgvs \
        --fasta ${vepfasta} \
        --offline \
        --cache ${vepcache} \
        --plugin Frameshift \
        --plugin Wildtype \
        --dir_plugins ${vepplugins}
    """

    stub:
    """
    cp -r \
        /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/stub_data/vepanno/* ./
    """
}
