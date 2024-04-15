#!/nemo/project/proj-tracerX/working/CMELA/alex/anaconda_nemo/bin/Rscript

############################
#FUNCTIONS / LIBRARIES 
############################

library(facets)
library(gdata)
library(dplyr)
library(ggplot2)

make.copy.number.plot.w.patient.heatmap = function(
    cn.df,
    ylim1 = c(0, 10),
    ybreaks = seq(0, 10, 2),
    colnames1 = c('sample', 'chr', 'start', 'end', 'cn', 'patient'),
    yintercept1 = 2
    ){
    colors = c(
        '0' = "#0c00ff",
        '1' = "#b38bff",
        '2' = "#ffffff",
        '3' = "#ffe4da",
        '4' = "#ffc8b6",
        '5' = "#ffac93",
        '6' = "#ff8f70",
        '7' = "#ff704f",
        '8' = "#ff4b2d",
        '9' = "#ff0000",
        '10' = "#ff0000",
        '11' ="#ff0000",
        '12' ="#ff0000",
        '13' ="#ff0000",
        '14' ="#ff0000"
    )
    #args:
        #cn.df: data.frame with columns sample, chr, start, end, cn and patient

    segmentation.plot = ggplot(cn.df, aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = cn)) +
        geom_rect() +
        #geom_text(data = sampdb2, aes(x = 90000000, y = 0.5, label = size), inherit.aes = F, size = 5) +
        facet_grid(cols = vars(chr), rows = vars(sample), space = 'free', scales = 'free') +
        rot.lab(size = 4) +
        scale_fill_manual(values = colors, drop = F) +
        theme(
            panel.spacing.y = unit(0, "lines"),
            panel.spacing.x = unit(0, "lines"),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            strip.text.y.right = element_text(angle = 0)
        )

    segmentation.plot
}

rot.lab = function(size, angle){
    if(missing(size)) size = 10
    if(missing(angle)) angle = 40
    theme(axis.text.x = element_text(angle = angle, hjust = 1, size = size))
}

############################
#CODE
############################

args <- commandArgs(trailingOnly = TRUE)

patient = args[[1]]
run_type = args[[2]]
rds.files = args[3:length(args)]

sample.names = gsub('(.*?)\\..*', '\\1', rds.files)

all.facets = lapply(rds.files, readRDS)

df1 = Map(function(x, y){
    z = y[[2]]$cncf
    z$patient = patient
    z$sample = x

    z2 = z[c('sample', 'chrom', 'start', 'end', 'tcn.em', 'patient')]
    z2
}, sample.names, all.facets)

df2 = bind_rows(df1)
df2$tcn.em = as.numeric(df2$tcn.em)
df2$tcn.em[df2$tcn.em > 14] = 14 #cap extremely high copy number levels for plotting
df2$tcn.em = factor(df2$tcn.em, levels = unique(sort(df2$tcn.em)))
colnames(df2) = c('sample', 'chr', 'start', 'end', 'cn', 'patient')

plot1 = make.copy.number.plot.w.patient.heatmap(df2)

ggsave(paste0(patient, '_', run_type, '_total.copy.number.plot.pdf'), plot = plot1, width = 10, height = 10)
