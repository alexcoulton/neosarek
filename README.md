# NeoSarek

A Nextflow pipeline designed for use subsequent to the nf-core Sarek pipeline
that performs multi-sample / multi-regional mutation calling, copy number
calling, HLA-typing, Neoantigen calling, and phylogenetic reconstruction.

# Usage

As copy number solutions need to be inspected manually, the pipeline is run in two stages.

Example first stage run:

```
#standard run
cd ~/work/cpi.nextflow/pipelines/mela/
nextflow ./main.nf \
    --samplesheet ~/work/ucl/data/mela.peace.new/pea004.ffpe.sarek.samplesheet.csv \
    --outputdir /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/output/pea004_ffpe_roq \
    --sarek_output_dir ~/work/cpi.nextflow/pipelines/sarek/output_pea004/ \
    -with-timeline run_report/pea004.ffpe.hla.roq.timeline.txt \
    -with-trace run_report/pea004.ffpe.hla.roq.trace.txt \
    -with-report run_report/pea004.ffpe.hla.roq.rep.html \
    -profile cluster \
    --stage 'initial' \
```

Example second stage run:

```
#standard run
cd ~/work/cpi.nextflow/pipelines/mela/
nextflow ./main.nf \
    --samplesheet ~/work/ucl/data/mela.peace.new/pea004.ffpe.sarek.samplesheet.csv \
    --outputdir /nemo/project/proj-tracerX/working/CMELA/alex/work/cpi.nextflow/pipelines/mela/output/pea004_ffpe_roq \
    --sarek_output_dir ~/work/cpi.nextflow/pipelines/sarek/output_pea004/ \
    -with-timeline run_report/pea004.ffpe.hla.roq.timeline.txt \
    -with-trace run_report/pea004.ffpe.hla.roq.trace.txt \
    -with-report run_report/pea004.ffpe.hla.roq.rep.html \
    -profile cluster \
    --stage 'review' \
```




# Detailed information


# Pipeline parameters

## Basic input parameters

### samplesheet
Samplesheet for NeoSarek - path to CSV file
### sarek_output_dir
Path of Sarek output directory
### outputdir
Path for NeoSarek output directory
### stub_data_dir
Path of stub data directory; a development parameter
### stage
String, specifies which stage of the pipeline to execute. Possible values are
`initial`, `review` or `facets_only`
### conipher_prefix
String. CONIPHER adds a prefix to all sample names, here we supply this prefix.
Example: `PEA`

## Tool resource parameters

### hlahd_install_directory
Path to HLAHD installation directory
### vepcache
Path to VEP cache directory, required for pVACseq
### vepplugins
Path to VEP plugin directory, required for pVACseq

## Genome resource parameters

### genome_reference
Path to reference genome fasta file
### genome_index
Path to reference genome index
### genome_dict
Path to reference genome dictionary
### germline_resource
Path to germline resource VCF file
### germline_resource_index
Path to germline resource index
### panel_of_normals
Path to panel of normals VCF file
### panel_of_normals_index
Path to panel of normals index
