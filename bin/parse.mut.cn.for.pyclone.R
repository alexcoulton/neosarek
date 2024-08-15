#!/nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/bin/Rscript
library(GenomicRanges)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

patient.name = args[[1]]
sample.name = args[[2]]
mutect2.vcf = args[[3]]
facets.output = args[[4]]
conipher.prefix = args[[5]]

#patient.name = 'PEA038'
#sample.name = 'SPA659A60'
#mutect2.vcf = './SPA659A60.forcecall.norm.vcf.gz'
#facets.output = 'SPA659A60.facets.rds'

#patient.name = 'PEA020'
#sample.name = 'SPA201A94'
#mutect2.vcf = './SPA201A94.forcecall.norm.isec.vcf.gz'
#facets.output = './SPA201A94.facets.rds'

read.vcf = function(vcf.path){
    x = read.delim(gzfile(vcf.path), sep = '\t', comment.char = '#', header = F)
    x2 = readLines(gzfile(vcf.path))
    header = gsub('^#', '', strsplit(x2[grep('^#CHROM', x2)], '\t')[[1]])
    colnames(x) = header
    x
}

mutect2.vcf = read.vcf(mutect2.vcf)
facets.output = readRDS(facets.output)

#facets.output[[1]]$jointseg
#facets.output[[1]]$out
#facets.output[[1]]$nX
#facets.output[[1]]$chromlevels
#facets.output[[1]]$dipLogR
#facets.output[[1]]$alBalLogR
#facets.output[[1]]$mafR.thresh
#facets.output[[1]]$flags

purity = facets.output[[2]]$purity
ploidy = facets.output[[2]]$ploidy

seg = facets.output[[2]]$cncf

mutations.granges = mutect2.vcf
mutations.granges = mutations.granges[c('CHROM', 'POS')]
mutations.granges$end = mutations.granges$POS
colnames(mutations.granges) = c('chr', 'start', 'end')

mutations.granges$chr = gsub('chr', '', mutations.granges$chr)
mutations.granges$chr = gsub('X', '23', mutations.granges$chr)

overlaps1 = as.data.frame(findOverlaps(
    GRanges(mutations.granges),
    GRanges(seg)
))

mutect2.vcf$major_cn = ""
mutect2.vcf$minor_cn = ""

#get major copy number by subtracting facets minor cn from total cn
mutect2.vcf$major_cn[overlaps1$queryHits] = seg$tcn.em[overlaps1$subjectHits] - seg$lcn.em[overlaps1$subjectHits]
#minor copy number is already present as a column in facets seg file
mutect2.vcf$minor_cn[overlaps1$queryHits] = seg$lcn.em[overlaps1$subjectHits]


# remove mutations where we do not have copy number info
## NAs
mutect2.vcf = mutect2.vcf[!(is.na(mutect2.vcf$minor_cn) | is.na(mutect2.vcf$major_cn)), ] 
## Blanks
mutect2.vcf = mutect2.vcf[!(mutect2.vcf$major_cn == "" & mutect2.vcf$minor_cn == ""), ]

#remove indels for pyclone
#think pyclone only deals with SNPs?
mutect2.vcf = mutect2.vcf[!(nchar(mutect2.vcf$REF) > 1 | nchar(mutect2.vcf$ALT) > 1), ]

mutation.id = with(mutect2.vcf, paste(CHROM, POS, REF, ALT, sep = '_'))

x = mutect2.vcf[[paste0(patient.name, '_', sample.name)]]

x2 = gsub('.*?:(.*?):.*', '\\1', x)
refn = as.numeric(gsub('(.*?),.*', '\\1', x2))
altn = as.numeric(gsub('.*?,(.*)', '\\1', x2))

depth = as.numeric(unlist(lapply(strsplit(x, ':'), function(y) y[[4]])))

############################
#CONIPHER INPUT 
############################

conipher.input = data.frame(
    CASE_ID = patient.name,
    SAMPLE = sample.name,
    CHR = mutect2.vcf$CHROM,
    POS = mutect2.vcf$POS,
    REF = mutect2.vcf$REF,
    ALT = mutect2.vcf$ALT,
    REF_COUNT = refn,
    VAR_COUNT = altn,
    DEPTH = depth,
    COPY_NUMBER_A = mutect2.vcf$major_cn,
    COPY_NUMBER_B = mutect2.vcf$minor_cn,
    ACF = purity,
    PLOIDY = ploidy
)

conipher.input$CHR = gsub('chr', '', conipher.input$CHR)
conipher.input$CHR = gsub('X', '23', conipher.input$CHR)
conipher.input$CHR = as.numeric(conipher.input$CHR)
conipher.input$SAMPLE = paste0(conipher.prefix, '_', conipher.input$SAMPLE)

write.table(
    conipher.input,
    paste0('./', patient.name, '.', sample.name, '.conipher.input.tsv'),
    quote = F,
    row.names = F,
    sep = '\t'
)

############################
#PYCLONE INPUT 
############################

#pyclonevi.input = data.frame(
    #mutation_id = mutation.id,
    #sample_id = sample.name,
    #ref_counts = refn,
    #alt_counts = altn,
    #major_cn = mutect2.vcf$major_cn,
    #minor_cn = mutect2.vcf$minor_cn,
    #normal_cn = 2,
    #tumour_content = purity
#)

pyclone.input = data.frame(
    mutation_id = mutation.id,
    ref_counts = refn,
    var_counts = altn,
    normal_cn = 2,
    minor_cn = mutect2.vcf$minor_cn,
    major_cn = mutect2.vcf$major_cn
)

pyclone.input = pyclone.input[pyclone.input$major_cn != 0, ]

which(duplicated(pyclone.input))
pyclone.input = pyclone.input[!duplicated(paste0(pyclone.input$mutation_id, '_', pyclone.input$sample_id)), ]

write.table(
    pyclone.input,
    paste0('./', patient.name, '.', sample.name, '.pyclone.input.tsv'),
    quote = F,
    row.names = F,
    sep = '\t'
)

writeLines(as.character(purity), paste0(patient.name, '.', sample.name, '.purity.txt'))
