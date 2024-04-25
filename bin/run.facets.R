#!/nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
tryCatch(args[[2]], error = function(e) 'No SNP pileup file provided')

sample = args[[1]]
snp.pileup.file = args[[2]]
rds.output = args[[3]]
plot.output = args[[4]]
diag.plot = args[[5]]
genome.version = args[[6]]

############################
#DEBUG 
############################

#sample = 'SPA817A1'
#snp.pileup.file = 'SPA817A1.facets.pileup.csv'
#rds.output = 'SPA817A1.facets.rds'
#plot.output = 'SPA817A1.facets.plot.pdf'
#diag.plot = 'SPA817A1.facets.diag.plot.pdf'

############################
#CONTINUE 
############################

library(facets)

snps1 = readSnpMatrix(snp.pileup.file)
snps2 = preProcSample(snps1, gbuild = genome.version)
snps3 = procSample(snps2)
fit = emcncf(snps3)

facets.solution = list(snps3, fit)

saveRDS(facets.solution, rds.output)

pdf(plot.output, 10, 10)
plotSample(x = snps3, emfit = fit)
dev.off()

pdf(diag.plot, 10, 10)
logRlogORspider(snps3$out, snps3$dipLogR)
dev.off()
