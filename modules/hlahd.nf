process hlahd {
    input:
    tuple val(patient), val(meta), path(merged_fastq_1), path(merged_fastq_2) 

    output:
    path output, emit: hlahd_output
    tuple val(patient), val("${meta['sample']}"), env(HLA), emit: hlaresult

    script:
    """
    mkdir output

    hlahd.sh \
        -t 8 \
        -m 100 \
        -c 0.95 \
        -f ${params.hlahd_install_directory}/freq_data/ \
        ${merged_fastq_1} \
        ${merged_fastq_2} \
        ${params.hlahd_install_directory}/HLA_gene.split.txt \
        ${params.hlahd_install_directory}/dictionary/ \
        ${meta['sample']} \
        output

    grep -E '^A|^B|^C' ./output/${meta['sample']}/result/${meta['sample']}_final.result.txt \
        | sed 's/^.\t//' \
        | sed 's/\t/,/g' \
        | tr ',' '\n' \
        | sed 's/:[^:]*\$//' \
        | tr '\n', ',' \
        | sed 's/.\$//g' > ./parsed.alleles.csv

    HLA=`cat ./parsed.alleles.csv`
    """

    stub:
    """
    #copy example output instead of running hlahd
    cp -r \
        /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/stub_data/hlahd/output \
        ./output

    grep -E '^A|^B|^C' ./output/${meta['sample']}/result/${meta['sample']}_final.result.txt \
        | sed 's/^.\t//' \
        | sed 's/\t/,/g' \
        | tr ',' '\n' \
        | sed 's/:[^:]*\$//' \
        | tr '\n', ',' \
        | sed 's/.\$//g' > ./parsed.alleles.csv

    HLA=`cat ./parsed.alleles.csv`
    """
}
