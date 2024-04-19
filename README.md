# NeoSarek

A Nextflow pipeline designed for use subsequent to the nf-core Sarek pipeline
that performs multi-sample / multi-regional mutation calling, copy number
calling, HLA-typing, Neoantigen calling, and phylogenetic reconstruction.

# Usage

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
