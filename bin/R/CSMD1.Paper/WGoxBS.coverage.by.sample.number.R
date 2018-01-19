#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
library(data.table)
library(ggplot2)
source("~/Pipelines/config/graphic.R")

dt.foo<-fread("/scratch/obsgynae/POPS/results/RnBeads/Boy.Girl/CpG/all/reports-2015-03-30_04_47_PM/quality_control_data/coverage_interrogated.csv")

p<-ggplot(dt.foo, aes(`Support (number of samples)`,`CpGs (million)`)) + 
	geom_point(aes(col=factor(`Minimal coverage`)),size=5) + 
	geom_line(aes(color=factor(`Minimal coverage`))) + 
	scale_x_reverse() + 
	scale_colour_manual(values=cbPalette,name="CpG Depth of\nCoverage Threshold") + 
	scale_y_continuous(breaks=seq(10,60,5)) + 
	theme_Publication()

pdf(file="~/results/CSMD1.2016.paper/Figures/Suppl/Suppl.Fig1.CpG.coverage.pdf", width=8, height=12, title="CpG depth of coverage")
print(p)
dev.off()

tiff(file="~/results/CSMD1.2016.paper/Figures/Suppl/Suppl.Fig1.CpG.coverage.tiff", width=8, height=12, units="in",res=300, compression = 'lzw', title="CpG depth of coverage")
print(p)
dev.off()
