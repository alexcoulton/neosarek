#!/bin/bash
ml GATK/4.4.0.0-GCCcore-12.3.0-Java-17

gatk Funcotator \
    --output ./PEA004.tumour.merged.indel200.norm.anno.vcf.gz \
    --ref-version hg38 \
    --data-sources-path ~/work/cpi.nextflow/funcotator/funcotator_dataSources.v1.7.20200521s \
    --output-file-format VCF \
    --variant /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/output/pea004_ffpe_roq/PEA004/mergedvcf/PEA004.tumour.merged.indel200.norm.vcf.gz \
    --reference ~/work/cpi.nextflow/igenome/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta


