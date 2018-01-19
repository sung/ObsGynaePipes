library(data.table)

#  
foo=list(`5kbp`="/scratch/obsgynae/POPS/results/RnBeads/Boy.Girl/CpG/all/AGA.Boy.Girl.tiling.top100.by.rank.csv",
         `gene`="/scratch/obsgynae/POPS/results/RnBeads/Boy.Girl/CpG/all/AGA.Boy.Girl.genes.top100.by.rank.csv",
         `promoter`="/scratch/obsgynae/POPS/results/RnBeads/Boy.Girl/CpG/all/AGA.Boy.Girl.promoters.top100.by.rank.csv",
         `cpgi`="/scratch/obsgynae/POPS/results/RnBeads/Boy.Girl/CpG/all/AGA.Boy.Girl.cpgislands.top100.by.rank.csv"
         )

dt.wgoxbs=lapply(foo, fread)
lapply(dt.wgoxbs, function(i) i[,V1:=NULL])

bar=list(`5kbp`="/scratch/obsgynae/POPS/results/RnBeads/SureSelect.Boy.Girl/CpG/all/reports-2016-02-04_04_21_PM/differential_methylation_data/diffMethTable_region_cmp1_tiling.csv.gz",
         `gene`="/scratch/obsgynae/POPS/results/RnBeads/SureSelect.Boy.Girl/CpG/all/reports-2016-02-04_04_21_PM/differential_methylation_data/diffMethTable_region_cmp1_genes.csv.gz",
         `promoter`="/scratch/obsgynae/POPS/results/RnBeads/SureSelect.Boy.Girl/CpG/all/reports-2016-02-04_04_21_PM/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv.gz",
         `cpgi`="/scratch/obsgynae/POPS/results/RnBeads/SureSelect.Boy.Girl/CpG/all/reports-2016-02-04_04_21_PM/differential_methylation_data/diffMethTable_region_cmp1_cpgislands.csv.gz")

dt.sureselect=lapply(bar, function(i) fread(paste("zcat ",i)))
lapply(dt.sureselect, function(i) i[,V1:=NULL])


dt.join=list()
for(i in names(foo)){
    dt.join[[i]]=merge(dt.wgoxbs[[i]], dt.sureselect[[i]], by=c("Chromosome","Start","End"), all.x=T)

}

lapply(names(dt.join), function(i) write.csv(
                                             dt.join[[i]][,.(Chromosome,Start,End,Strand,mean.mean.g1,mean.mean.g2,
                                                             mean.mean.diff.x,mean.mean.quot.log2.x,comb.p.val.x,
                                                             comb.p.adj.fdr.x,num.sites.x,combinedRank.x, 
                                                             mean.mean.F,mean.mean.M,mean.mean.diff.y,
                                                             mean.mean.quot.log2.y,comb.p.val.y,comb.p.adj.fdr.y,
                                                             num.sites.y,combinedRank.y)][order(combinedRank.x)],
                                             file=file.path("~/results/CSMD1.2016.paper/RData/SureSelect",paste0(i,".WGoxBS.top100.SureSelect.csv"))
                                             ))


load("~/data/Annotation/Ensembl/GRCh37/dt.ensg.RData") # dt.ensg, dt.enst, gr.ensg, gr.enst, gr.exon
dt.ensg[,chromosome_name:=paste0("chr",chromosome_name)] # 1 => chr1
#dt.ensg[,start_position:=start_position-1500]
#dt.ensg[,end_position:=end_position+1500]
setkeyv(dt.ensg, c("chromosome_name","start_position","end_position")) # regions of interests 

for(i in c("5kbp","cpgi")){
    dt.query=dt.wgoxbs[[i]]
    setkeyv(dt.query, c("Chromosome","Start","End")) # regions of interests 
    system.time(dt.overlap<-foverlaps(dt.query, dt.ensg, type="any", nomatch=NA)) # query=dt.wgoxbs[["5kbp"]]

    write.csv(dt.overlap[,.(hgnc_symbols=paste(hgnc_symbol, collapse=","), ensgs=paste(ensembl_gene_id,collapse=","), 
                combinedRank=mean(combinedRank)),
                    by=.(Chromosome,Start,End)][order(combinedRank)],
        file=file.path("~/results/CSMD1.2016.paper/RData/SureSelect",paste0(i,".WGoxBS.top100.SureSelect.gene.name.csv")))

}
