#!/usr/bin/Rscript --vanilla

# FLD
my.files=c(
`set1`="~/results/RnBeads/SGA.AGA.FLD.v1.8/CpG/all/reports-2015-04-22_04_20_PM/differential_methylation_data/diffMethTable_site_cmp1.csv.gz",
`set2`="~/results/RnBeads/SGA.AGA.FLD.min.50/CpG/all/reports-2015-04-23_05_08_PM/differential_methylation_data/diffMethTable_site_cmp1.csv.gz",
`set3`="~/results/RnBeads/SGA.AGA.FLD.min.50/CpG/all/reports-2015-04-23_06_42_PM/differential_methylation_data/diffMethTable_site_cmp1.csv.gz"
)

diff.sites<-lapply(my.files, read.csv)

diff.sites<-lapply(diff.sites, function(df) subset(df, select=-c(id))) # drop "id" column

diff.sites<-lapply(diff.sites, function(df) cbind(df, End=df$Start+1)) # add "End" column

diff.sites<-lappy(diff.sites, makeGRangesFromDataFrame) # from GenomicRanges
