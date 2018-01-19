#/usr/bin/Rscript --vanilla
rm(list=ls())
library(ggplot2)

pdf ( file=paste0(format(Sys.time(), "%Y-%m-%d_%I%p"), ".pdf") )

#methyl<-read.table("~/data/methyl_seq_mapping.csv", header=TRUE)
methyl<-read.table("~/Pipelines/data/methyl_seq_mapping.txt", header=TRUE)
qplot(Raw.Read, DoC, data=methyl, group=SLX, color=SLX) + scale_x_continuous("NO Raw Reads (10^6)") + scale_y_continuous("Depth of CpG Coverage")
qplot(Raw.Read, BoC_1x, data=methyl, group=SLX, color=SLX) + scale_x_continuous("NO Raw Reads (10^6)") + scale_y_continuous("% CpG Coverage at 1x")
qplot(Raw.Read, BoC_3x, data=methyl, group=SLX, color=SLX) + scale_x_continuous("NO Raw Reads (10^6)") + scale_y_continuous("% CpG Coverage at 3x")
qplot(Raw.Read, BoC_5x, data=methyl, group=SLX, color=SLX) + scale_x_continuous("NO Raw Reads (10^6)") + scale_y_continuous("% CpG Coverage at 5x")
qplot(Raw.Read, MET_1x, data=methyl, group=SLX, color=SLX) + scale_x_continuous("NO of Raw Reads (10^6)") + scale_y_continuous("% methylated CpG at 1x")
qplot(Raw.Read, MET_3x, data=methyl, group=SLX, color=SLX) + scale_x_continuous("NO of Raw Reads (10^6)") + scale_y_continuous("% methylated CpG at 3x")
qplot(Raw.Read, MET_5x, data=methyl, group=SLX, color=SLX) + scale_x_continuous("NO of Raw Reads (10^6)") + scale_y_continuous("% methylated CpG at 5x")

dev.off()
