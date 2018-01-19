#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(methylKit)

top_result_dir<-"/scratch/ssg29/results"

# filtering raw CpG
min.count<-7
max.perc<-99.9
# differentially methylated region
myDiffvalue=25
myQvalue=0.01
# For DEG
myDegDiffvalue=5
myDegQvalue=0.01

# cores to use for the diffmet
np<-16
# promoter definitions
promoter.up=1500
promoter.down=500

# annotation files 
gene.file='~/data/Annotation/ENSG.hg19.bed'
cgi.file='~/data/Annotation/cpgi.hg19.bed'

# input definitions for 8 oxBS samples
my_prefix='8oxBS'
#SLX-8074.SLX-8080.v1
#oxBS: 
#A003 (68, SGA, D705_D505)
#A001 (66, AGA, D705_D506)
#BS
#A012 (68, SGA, D705_D505)
#A010 (66, AGA, D705_D506)
project1<-'SLX-8074.SLX-8080'
project_dir1<-paste0(top_result_dir,"/",project1,".v1") # /scratch/ssg29/results/SLX-8074.v1
bismark_call_dir1<-paste0(project_dir1,"/Bismark/MethylCall") # /scratch/ssg29/results/SLX-8074.v1/Bismark/MethylCall
bissnp_call_dir1<-paste0(project_dir1,"/BisSNP") # /scratch/ssg29/results/SLX-8074.v1/BisSNP

#SLX-8075.SLX-8077.SLX-8081.v1 (oxBS only)
#D705_D503 74 A007 SGA 
#D705_D504 72 A005 AGA
#D706_D501 77 A004 SGA
#D706_D502 75 A002 AGA
#D706_D504 65 A006 SGA
#D706_D505 73 A008 AGA
project2<-'SLX-8075.SLX-8077.SLX-8081'
project_dir2<-paste0(top_result_dir,"/",project2,".v1") 
bismark_call_dir2<-paste0(project_dir2,"/Bismark/MethylCall") 
bissnp_call_dir2<-paste0(project_dir2,"/BisSNP") 

# ~/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/MethylKit
# ~/scratch/results/Methyl-Seq/SGA.AGA/methylKit
methylkit_dir<-"~/scratch/results/Methyl-Seq/SGA.AGA/methylKit"

# edgeR result for 8 oxBS samples
toptags=paste0(top_result_dir,'/SLX-8549.v1/edgeR.oxBS/SLX-8549.toptags_pair_edgeR.csv')
deg_PValue <- 1.0e-02
