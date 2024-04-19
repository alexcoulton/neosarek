#!/nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/bin/Rscript

############################
#INPUTS 
############################

args <- commandArgs(trailingOnly = TRUE)
args2 = unlist(args)

patient = args2[[1]]
vcfs = args2[grepl('.vcf.gz$', args2)]
rds.files = args2[grepl('.rds$', args2)]
pyclone.file = args2[grepl('.tsv$', args2)]

############################
#FUNCTIONS / LIBRARIES 
############################

library(jsonlite)
library(reshape2)
library(GenomicRanges)
library(dplyr)

read.vcf = function(x){
    header = readLines(gzfile(x))
    header = gsub('^#', '', header[grepl('^#CHROM', header)])
    header = strsplit(header, '\t')[[1]]
    vcf = read.delim(gzfile(x), sep = '\t', comment.char = '#', header = F)
    colnames(vcf) = header
    vcf$key = paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, sep = '_')
    vcf
}

############################
#DEBUG 
############################

patient = 'PEA004'
vcfs = list.files('./', pattern = '.vcf.gz$')
rds.files = list.files('./', pattern = '.rds$')
pyclone.file = list.files('./', pattern = '.tsv$')

############################
#BUILD SSM FILE
############################

samples = gsub('(.*?)\\..*', '\\1', vcfs)

print(patient)
print(vcfs)
print(rds.files)
print(pyclone.file)

vcfs2 = lapply(vcfs, read.vcf)

#process samples one at a time within this patient
pairtree.data = lapply(samples, function(x){
    vcf = vcfs[grepl(paste0(x, '\\.forcecall'), vcfs)]
    rds = rds.files[grepl(paste0(x, '\\.facets'), rds.files)]

    vcf = read.vcf(vcf)
    rds = readRDS(rds)

    key.cols = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'key')

    vcf.melt = melt(vcf, id = key.cols)
    vcf.melt = vcf.melt[vcf.melt$variable == paste0(patient, '_', x), ] # filter out normal sample

    split.vals = strsplit(vcf.melt$value, ':')

    info1 = lapply(split.vals, function(x){
        gt = x[[1]]
        nref = strsplit(x[[2]], ',')[[1]][1]
        nalt = strsplit(x[[2]], ',')[[1]][2]
        if(is.na(nalt)) nalt = '.'
        depth = x[[4]]

        data.frame(
            gt,
            nref,
            nalt,
            depth
        )
    }) %>% bind_rows

    vcf.melt2 = bind_cols(vcf.melt, info1)

    vcf.melt3 = vcf.melt2[c(
        'key',
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'variable',
        'gt',
        'nref',
        'nalt',
        'depth'
    )]

    ############################
    #FACETS 
    ############################
    cn = rds[[2]]$cncf

    mut = vcf.melt3
    mut$start = mut$POS
    mut$end = mut$POS

    mut$CHROM = gsub('chr', '', mut$CHROM)
    mut$CHROM = gsub('X', '23', mut$CHROM)
    mut$CHROM = gsub('Y', '24', mut$CHROM)

    mut.ranges = GRanges(mut)
    cn.ranges = GRanges(cn)

    overlaps1 = as.data.frame(findOverlaps(mut.ranges, cn.ranges))

    mut2 = mut[overlaps1$queryHits, ]
    cn2 = cn[overlaps1$subjectHits, ]

    mut2$total.cn = cn2$tcn.em
    mut2$minor.cn = cn2$lcn.em

    mut2$purity = rds[[2]]$purity
    mut2$ploidy = rds[[2]]$ploidy

    pyclone.results = read.delim(pyclone.file, sep = '\t')
    py2 = pyclone.results[pyclone.results$sample_id == x, ]

    py2$key = py2$mutation_id

    mut2$key %in% py2$mutation_id

    mut3 = left_join(mut2, py2)

    mut3$vaf = as.numeric(mut3$nalt) / as.numeric(mut3$depth)
    mut3$ccf = mut3$cellular_prevalence

    mut3$mut.cn = 1
    mut3$mut.cn[mut3$nalt == 0] = 0

    #formula rearranged from Turajlic et al., 2018
    mut3$mut.cn = round((mut3$vaf * (((1 - mut3$purity) * 2) + (mut3$total.cn * mut3$purity))) / mut3$purity * mut3$ccf)

    mut3$mut.cn[mut3$mut.cn == 0 & mut3$nalt >= 5] = 1

    mut4 = mut3[!mut3$nalt < 5, ]
    mut4 = mut4[!is.na(mut4$key), ]
    mut4 = mut4[!is.na(mut4$ccf), ]
    mut4$major.cn = mut4$total.cn - mut4$minor.cn

    if(nrow(mut4[mut4$mut.cn > mut4$major.cn, ]) > 0){
        mut4$mut.cn[mut4$mut.cn > mut4$major.cn] =
            mut4$major.cn[mut4$mut.cn > mut4$major.cn]
    }

    mut4[mut4$mut.cn > 1, ][c('key', 'total.cn', 'minor.cn', 'mut.cn', 'vaf', 'ccf')]

    mut4$sample = x
    mut4$pairtree.prob = mut4$mut.cn / mut4$total.cn
    mut4
}) %>% bind_rows

write.table(pairtree.data, paste0('./', patient, '.combined.mutations.tsv'), row.names = F, quote = F, sep = '\t')

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

