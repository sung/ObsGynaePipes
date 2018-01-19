#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 18/Oct/2016
# Last modified 16/Nov/2016

source("~/Pipelines/bin/R/CSMD1.Paper/local.R")

##################
## Main Figures ##
##################
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Fig1.R") # beta-distribution
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Fig2.R") # Meth-diff boxplot
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Fig3.R") # Manhanttan plot
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Fig4.R") # CSMD1 DMR + 2 exons (5')
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Fig5.R") # Methylation and expression profile over the full-length CSMD1
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Fig6.R") # Exon-level pheatmap across multiple tissues

########################### 
## Supplementary Figures ##
########################### 
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Suppl.Fig.putative.promoter.DMR.R") # CSMD1 DMR from the putative promoter
source("~/Pipelines/bin/R/CSMD1.Paper/CSMD1.Suppl.Fig.Exon.Level.Exp.R") # P-value from exon-level expression & read coverage
