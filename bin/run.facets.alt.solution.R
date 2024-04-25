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

############################
#CALCULATE ALTERNATE SOLUTION BASED ON MOST BALANCED SEGMENT 
############################

x = facets.solution

x[[2]]$cncf$length.mb = 
    (x[[2]]$cncf$end - x[[2]]$cncf$start) / 1000000

#remove lower quartile of segments by length
z = x[[2]]$cncf[x[[2]]$cncf$length.mb > quantile(x[[2]]$cncf$length.mb, seq(0, 1, 0.25))[[2]], ]

z$abs.mafR = abs(z$mafR)

#get top 20% of most balanced segs
balanced.segs = z[z$abs.mafR <= quantile(z$abs.mafR, seq(0, 1, 0.1))[[3]], ]

balanced.segs.cnlr = sort(balanced.segs$cnlr.median)

############################
#POTENTIAL dipLogR VALUES 
############################

new.dip = balanced.segs.cnlr[which.min(abs(sort(balanced.segs.cnlr)))]
#median.diplogr = median(x[[2]]$cncf$cnlr.median)

snps1 = readSnpMatrix(snp.pileup.file)
snps2 = preProcSample(snps1)
snps3 = procSample(snps2, dipLogR = new.dip)
fit = emcncf(snps3)

facets.solution = list(snps3, fit)

saveRDS(facets.solution, rds.output)

pdf(plot.output, 10, 10)
plotSample(x = snps3, emfit = fit)
dev.off()

pdf(diag.plot, 10, 10)
logRlogORspider(snps3$out, snps3$dipLogR)
dev.off()
