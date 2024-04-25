process HLAHD {
    input:
    tuple val(patient), val(sample), path(merged_fastq_1), path(merged_fastq_2) 

    publishDir "${params.outputdir}/${patient}/hlahd/", mode: params.publish_dir_mode, overwrite: params.publish_dir_overwrite


    output:
    path output, emit: hlahd_output
    tuple val(patient), val(sample), env(HLA), emit: hlaresult

    script:
    """
    mkdir output
    echo ""

    hlahd.sh \
        -t 8 \
        -m 100 \
        -c 0.95 \
        -f ${params.hlahd_install_directory}/freq_data/ \
        ${merged_fastq_1} \
        ${merged_fastq_2} \
        ${params.hlahd_install_directory}/HLA_gene.split.txt \
        ${params.hlahd_install_directory}/dictionary/ \
        ${sample} \
        output

    grep -E '^A|^B|^C' ./output/${sample}/result/${sample}_final.result.txt \
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
    ln -s ${params.stub_data_dir}/${patient}/hlahd/* ./

    grep -E '^A|^B|^C' ./output/${sample}/result/${sample}_final.result.txt \
        | sed 's/^.\t//' \
        | sed 's/\t/,/g' \
        | tr ',' '\n' \
        | sed 's/:[^:]*\$//' \
        | tr '\n', ',' \
        | sed 's/.\$//g' > ./parsed.alleles.csv

    HLA=`cat ./parsed.alleles.csv`
    """
}
