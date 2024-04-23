# Subworkflow documentation

## FACETS_CN_CALLING

This subworkflow performs copy number calling with Facets.

## HLA_NEOANTIGEN

This subworkflow performs HLA typing with HLA-HD and neoantigen calling with
pVACseq

## PYCLONE

This subworkflow performs clustering of mutations with Pyclone-VI and tree
inference with Pairtree

# Process documentation

## NORMALIZE_FILTER_VCF
## MERGE_VCF
## MUTECT2_FORCECALL
## PREPARE_PYCLONE_INPUT
## MERGE_PYCLONE_INPUT_PER_PATIENT
## MERGE_CONIPHER_INPUT_PER_PATIENT
## RUN_PYCLONE
## RUN_PAIRTREE
## PREPARE_PAIRTREE_INPUT
## RUN_CONIPHER
## MERGE_FASTQ
## HLAHD
## VEP_ANNO
## PVACSEQ
## FACETS_SNP_PILEUP
## FACETS_RUN
## FACETS_RUN_ALTERNATE_SOLUTION
Gererates an alternate Facets solution. Typically, Facets determines the value
of the dipLogR parameter algorithmically, and uses this to generate copy number
profiles. However, this value is not always appropriate. Here we generate an
alternate solution based on an alternate dipLogR value, which is calculated as follows:

1. Remove lowest quartile of segments by length
2. Get the top 20% of most balanced segments
3. Set the dipLogR as the cnlr.median of the segment where this value is lowest

## FACETS_DIAGNOSTICS
Generates diagnostic plots for the original Facets solution
## FACETS_DIAGNOSTICS_ALT
Generates diagnostic plots for the alternate Facets solution





