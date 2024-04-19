#!/opt/conda/envs/conipher/bin/Rscript

#library(CONIPHER)
library(tidyverse)

devtools::load_all('./CONIPHER')

args <- commandArgs(trailingOnly = TRUE)

patient = args[[1]]
conipher.prefix = args[[2]]
input_tsv = args[[3]]

#devtools::load_all('/conipherdir')
#patient = 'PEA038'
#conipher.prefix = 'PEA'
#input_tsv = 'test.tsv'

#export NUMBA_CACHE_DIR=/faildir/numba_cache
#mkdir /faildir/numba_cache

out_dir = paste0('./', patient, '_conipher_results/')

conipher_run(
    case_id = patient,
    prefix = conipher.prefix,
    out_dir = out_dir,
    input_tsv_loc = input_tsv,
    nProcs = 8
)
