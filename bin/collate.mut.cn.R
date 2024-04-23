#!/nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/bin/Rscript

#collates mutation data per patient (combines all samples into one file)
#estimates mutation CN
#performs some filtering of mutations (AD > 5)

#collates copy number data (segments, purity, ploidy) into one file per patient

############################
#INPUTS 
############################

args <- commandArgs(trailingOnly = TRUE)
args2 = unlist(args)

patient = args2[[1]]
vcfs = args2[grepl('.vcf.gz$', args2)]
rds.files = args2[grepl('.rds$', args2)]
pyclone.file = args2[grepl('.tsv$', args2)]

vcfs = vcfs[!grepl('anno.vcf.gz$', vcfs)]
annotations = args2[grepl('anno.vcf.gz$', args2)]

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

#patient = 'PEA004'
#vcfs = list.files('./', pattern = '.vcf.gz$')
#rds.files = list.files('./', pattern = '.rds$')
#pyclone.file = list.files('./', pattern = 'results.tsv$')
#annotations = list.files('./', pattern = 'anno.vcf.gz$')

############################
#COLLATE MUTATION DATA
############################

samples = gsub('(.*?)\\..*', '\\1', vcfs)

print(patient)
print(vcfs)
print(rds.files)
print(pyclone.file)

vcfs2 = lapply(vcfs, read.vcf)

#process samples one at a time within this patient
count1 = 1
mut.data = lapply(samples, function(x){
    print(count1)
    count1 <<- count1 + 1

    print('a1')

    vcf = vcfs[grepl(paste0(x, '\\.forcecall'), vcfs)]
    rds = rds.files[grepl(paste0(x, '\\.facets'), rds.files)]

    vcf = read.vcf(vcf)
    rds = readRDS(rds)

    key.cols = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'key')

    vcf.melt = melt(vcf, id = key.cols)
    vcf.melt = vcf.melt[vcf.melt$variable == paste0(patient, '_', x), ] # filter out normal sample

    split.vals = strsplit(vcf.melt$value, ':')

    print('a2')
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

    #vcf.melt3 = vcf.melt2[c(
        #'key',
        #'CHROM',
        #'POS',
        #'REF',
        #'ALT',
        #'INFO',
        #'variable',
        #'gt',
        #'nref',
        #'nalt',
        #'depth'
    #)]

    ############################
    #FACETS 
    ############################
    print('a3')
    cn = rds[[2]]$cncf

    mut = vcf.melt2
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

    print('a4')
    mut2$purity = rds[[2]]$purity
    mut2$ploidy = rds[[2]]$ploidy

    pyclone.results = read.delim(pyclone.file, sep = '\t')
    py2 = pyclone.results[pyclone.results$sample_id == x, ]

    py2$key = py2$mutation_id

    mut2$key %in% py2$mutation_id

    print('a5')
    mut3 = left_join(mut2, py2)

    mut3$vaf = as.numeric(mut3$nalt) / as.numeric(mut3$depth)
    mut3$ccf = mut3$cellular_prevalence

    mut3$mut.cn = 1
    mut3$mut.cn[mut3$nalt == 0] = 0

    #formula rearranged from Turajlic et al., 2018
    mut3$mut.cn = round((mut3$vaf * (((1 - mut3$purity) * 2) + (mut3$total.cn * mut3$purity))) / mut3$purity * mut3$ccf)

    mut3$mut.cn[mut3$mut.cn == 0 & mut3$nalt >= 5] = 1

    print('a6')
    mut4 = mut3[!mut3$nalt < 5, ]
    mut4 = mut4[!is.na(mut4$key), ]
    mut4 = mut4[!is.na(mut4$ccf), ]
    mut4$major.cn = mut4$total.cn - mut4$minor.cn

    if(nrow(mut4[mut4$mut.cn > mut4$major.cn, ]) > 0){
        mut4$mut.cn[mut4$mut.cn > mut4$major.cn] =
            mut4$major.cn[mut4$mut.cn > mut4$major.cn]
    }

    print('a7')
    #mut4[mut4$mut.cn > 1, ][c('key', 'total.cn', 'minor.cn', 'mut.cn', 'vaf', 'ccf')]

    mut4$sample = x
    mut4$pairtree.prob = mut4$mut.cn / mut4$total.cn
    mut4
}) %>% bind_rows

mut.data$rescue = T
mut.data$rescue[grep('DENOVO_CALLED', mut.data$INFO)] = F

############################
#ADD FUNCOTATOR ANNOTATIONS 
############################

anno2 = read.vcf(annotations)
anno3 = distinct(anno2[c('key', 'INFO')])

funco = unlist(lapply(strsplit(anno2$INFO, ';'), function(x) x[[5]]))
funco = gsub('\\]', '', gsub('FUNCOTATION=\\[', '', funco))
funco2 = strsplit(funco, '\\|')

funco.gene = unlist(lapply(funco2, function(x) x[[1]]))
funco.mut.type = unlist(lapply(funco2, function(x) x[[6]]))
funco.cat = unlist(lapply(funco2, function(x) x[[8]]))

anno.df = data.frame(
    key = anno3$key,
    funco.gene,
    funco.mut.type,
    funco.cat
)

mut.data2 = left_join(mut.data, anno.df, by = 'key')

############################
#WRITE DATA 
############################

write.table(mut.data2, paste0('./', patient, '.combined.mutations.tsv'), row.names = F, quote = F, sep = '\t')

############################
#PROCESS SEGMENTS 
############################

comb.seg = lapply(rds.files, function(x){
    x2 = readRDS(x)
    cn = x2[[2]]$cncf
    x.sample = gsub('.facets.rds$', '', x)
    cn$sample = x.sample
    cn$patient = patient
    cn$purity = x2[[2]]$purity
    cn$ploidy = x2[[2]]$ploidy
    cn
}) %>% bind_rows

write.table(
    comb.seg,
    paste0('./', patient, '.combined.seg.tsv'),
    row.names = F,
    quote = F,
    sep = '\t'
)
