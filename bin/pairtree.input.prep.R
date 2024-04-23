#!/nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/bin/Rscript

############################
#INPUTS 
############################

args <- commandArgs(trailingOnly = TRUE)
args2 = unlist(args)

patient = args2[[1]]
pairtree.data = read.delim(args2[[2]], sep = '\t')

#patient = 'PEA004'
#pairtree.data = read.delim('PEA004.combined.mutations.tsv', sep = '\t')

############################
#FUNCTIONS / LIBRARIES 
############################

library(jsonlite)
library(reshape2)
library(GenomicRanges)
library(dplyr)


############################
#PAIRTREE FILE PREP 
############################

pairtree.data$CHROM = as.numeric(pairtree.data$CHROM)

z = arrange(pairtree.data, CHROM, POS, sample_id)
samps1 = sort(unique(pairtree.data$sample))

nalt.vals = lapply(unique(z$key), function(i){
    x = lapply(samps1, function(ii){
        g = pairtree.data[pairtree.data$key == i & pairtree.data$sample_id == ii, ]
        if(nrow(g) == 0) return(0)
        return(as.numeric(g$nalt))
    }) %>% paste0(collapse = ',')
}) %>% unlist

total.reads.vals = lapply(unique(z$key), function(i){
    lapply(samps1, function(ii){
        g = pairtree.data[pairtree.data$key == i & pairtree.data$sample_id == ii, ]
        if(nrow(g) == 0) return(0)
        return(as.numeric(g$depth))
    }) %>% paste0(collapse = ',')
}) %>% unlist

prob.vals = lapply(unique(z$key), function(i){
    lapply(samps1, function(ii){
        g = pairtree.data[pairtree.data$key == i & pairtree.data$sample_id == ii, ]
        if(nrow(g) == 0) return(0.5) #pairtree prob can't be 0
        return(as.numeric(g$pairtree.prob))
    }) %>% paste0(collapse = ',')
}) %>% unlist

df1 = data.frame(
    id = paste0('s', 0:(length(unique(z$key)) - 1)),
    name = unique(z$key),
    var_reads = nalt.vals,
    total_reads = total.reads.vals,
    var_read_prob = prob.vals
)

write.table(df1, paste0(patient, '.pairtree.ssm'), sep = '\t', quote = F, row.names = F)

############################
#BUILD JSON FILE 
############################

p.split = split(pairtree.data, pairtree.data$key)

#check that each variant belongs to the same cluster across all samples
lapply(p.split, function(x){
    length(unique(x$cluster_id))
}) %>% unlist %>% unname %>% table

df2 = df1[c('id', 'name')]
colnames(df2) = c('ssm.id', 'key')
pairtree.data = left_join(pairtree.data, df2)

g = distinct(pairtree.data[c('ssm.id', 'cluster_id')]) %>%
    arrange(cluster_id, as.numeric(gsub('s', '', ssm.id)))

g.split = split(g, g$cluster_id)
g.split2 = lapply(g.split, function(x) x$ssm.id)

pairtree.json = list(
    samples = unique(pairtree.data$sample_id),
    clusters = g.split2,
    garbage = ""
)

pairtree.json[[2]] = unname(pairtree.json[[2]])
pairtree.json2 <- toJSON(pairtree.json)
pairtree.json2 = gsub('"garbage":\\[""\\]', '"garbage":[]', pairtree.json2)

write(pairtree.json2, file = paste0(patient, ".pairtree.input.json"))

