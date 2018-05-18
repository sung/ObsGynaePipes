deseq.RData="~/results/RNA-Seq/Plasma.2017/DESeq2/ALL.featureCount.PE75/deseq.ALL.RData"

dt.count<-data.table(`ensembl_gene_id`=rownames(counts(dds)),counts(dds))
dt.meanCount<-data.table(`ensembl_gene_id`=rownames(counts(dds)), `Plasma.all`=rowMeans(counts(dds)))

mySum=merge(dt.meanCount, dt.ensg)[,sum(Plasma.all)] # 19966212
merge(dt.meanCount, dt.ensg)[,.(`read`=round(sum(Plasma.all),1),`read.ratio`=round(sum(Plasma.all)/mySum*100,1),`No.gene`=.N),gene_biotype]

merge(dt.count, dt.ensg)
merge(dt.count, dt.ensg)[gene_biotype=="Mt_rRNA"]

write.csv(merge(dt.meanCount, dt.ensg)[,.(`read`=round(sum(Plasma.all),1),`read.ratio`=round(sum(Plasma.all)/mySum*100,1),`No.gene`=.N),gene_biotype], file=paste0("~/results/RNA-Seq/",myProject,"/Meta/meta.",myProject,".count.by.gene_biotype.csv"))
