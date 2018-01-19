# !/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
myCaller='Etc'
source ("~/Pipelines/config/DEG.R") # load config

# SE125 RData file
load("~/results/RNA-Seq/SGA.AGA.total.SE125/DESeq2/LOW.PAPPA/deseq.LOW.PAPPA.RData")
SE125<-list(FPKM=fpkm(ddsExonicGene), # isa 'matrix'
			FPM=fpm(ddsExonicGene), # isa 'matrix'
			rlogMat=assay(rld), # isa 'matrix'
			countMat=counts(dds), # isa 'matrix'
			res=res
		)
load("~/results/RNA-Seq/SGA.AGA.total.SE50/DESeq2/LOW.PAPPA/deseq.LOW.PAPPA.RData")
SE50<-list(FPKM=fpkm(ddsExonicGene), # isa 'matrix'
			FPM=fpm(ddsExonicGene), # isa 'matrix'
			rlogMat=assay(rld), # isa 'matrix'
			countMat=counts(dds), # isa 'matrix'
			res=res
		)


# Pearson's Correlation Co-efficient of FPKM between SE125 and SE50
sapply(colnames(SE125$FPKM), function(i) cor(SE125$FPKM[,i],SE50$FPKM[,i]))

# FPKM|FPM|Count of this sample
my.sample="71"
foo<-data.frame(fpkm.SE125=SE125$FPKM[,my.sample], fpkm.SE50=SE50$FPKM[,my.sample]) # FPKM
bar<-data.frame(fpm.SE125=SE125$FPM[,my.sample], fpm.SE50=SE50$FPM[,my.sample]) # FPM
foo<-cbind(foo, bar)
koo<-data.frame(count.SE125=SE125$countMat[,my.sample], count.SE50=SE50$countMat[,my.sample]) # raw-count
foo<-cbind(foo, koo[rownames(foo),])

#keep = rownames(bar[rowSums(bar >1) >= ncol(bar),]) # genes of FPM>1
keep=TRUE # all genes (no filter)

my.title=paste0("Pearson's r=",round(cor(foo)["fpkm.SE125","fpkm.SE50"],2), " (Sample: ",my.sample,")") # Pearson-corellation coefficieny
p1<-ggplot(foo, aes(log2(fpkm.SE125+1),log2(fpkm.SE50+1))) + 
	geom_point(aes(colour = fpm.SE50+fpm.SE125), alpha = 0.4, size=5) + 
	ggtitle(my.title) +
	scale_colour_gradient2(limit=c(0,2))
tiff(filename=paste0("~/results/RNA-Seq/SGA.AGA/SE125.SE50.FPKM.",my.sample,".tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
print(p1)
dev.off()

# Annotate Genes of FPKM.SE50>4 & FPKM.SE125<1
dummy<-log2(foo+1)
#select<-(dummy$fpkm.SE125<4 & dummy$fpkm.SE50>5) | (dummy$fpkm.SE125<1 & dummy$fpkm.SE50>4) | (dummy$fpkm.SE125>6 & dummy$fpkm.SE50<4)
select<-abs(dummy$fpkm.SE50-dummy$fpkm.SE125)>3
anno<-getBM(attributes = my.fields, filters = my.filter, values = rownames(dummy[select,]), mart = grch37)
anno<-merge(cbind(foo[select,], ensembl_gene_id=rownames(foo[select,])), anno)

write.csv(anno, file=paste0("~/results/RNA-Seq/SGA.AGA/SE125.SE50.FPKM.",my.sample,".outliers.csv"))

p2<-ggplot(foo, aes(log2(fpkm.SE125+1),log2(fpkm.SE50+1))) + 
	geom_point(colour = 'grey', alpha = 0.4, size=5) + 
	geom_point(data=anno, aes(log2(fpkm.SE125+1),log2(fpkm.SE50+1), colour=gene_biotype), alpha=0.4, size=5) 
	#geom_text(data=anno, aes(log2(fpkm.SE125+1),log2(fpkm.SE50+1), label=paste(hgnc_symbol,gene_biotype,sep=":")), hjust=0, vjust=0, size=5)
tiff(filename=paste0("~/results/RNA-Seq/SGA.AGA/SE125.SE50.FPKM.",my.sample,".outliers.tiff"),width=10,height=9,units="in",res=300, compression = 'lzw')
print(p2)
dev.off()


# Colour genes by their gene_biotype
my.group.by<-"gene_biotype"
anno <- getBM(attributes = my.fields, filters = my.filter, values = rownames(foo), mart = grch37)
anno.collapsed <- plyr::ddply(anno, my.filter, function(df) paste(unique(df[[my.group.by]]), collapse=';')) 
colnames(anno.collapsed) <- c(my.filter,my.group.by)

foo[,my.filter]<-rownames(foo) # add a new column id
foo<-merge(foo, anno.collapsed, by=my.filter,  all.x=TRUE) # left join to add column 'hgnc_symbol'
rownames(foo)<-foo[[my.filter]]
foo[[my.filter]]<-NULL; 

p2<-ggplot(foo, aes(log2(fpkm.SE125+1),log2(fpkm.SE50+1))) + 
	geom_point(aes(colour = factor(gene_biotype)), alpha = 0.3, size=4) +
	ggtitle(my.title)


# FPKM SE125 vs. FPKM SE50
p3<-ggplot(foo[keep,], aes(fpkm.SE125,fpkm.SE50)) + geom_point(alpha = 1/10, size=4) 


# rlog of this sample
kar<-data.frame(rlog.SE125=SE125$rlogMat[,my.sample], rlog.SE50=SE50$rlogMat[,my.sample])
p4<-ggplot(kar, aes(rlog.SE125,rlog.SE50)) + geom_point(alpha = 1/10, size=4) 

# p-value of this set
dummy<-data.frame(pval.SE125=SE125$res[,"pvalue"], p.valSE50=SE50$res[,"pvalue"])
p5<-ggplot(dummy, aes(pval.SE125,pval.SE50)) + geom_point(alpha = 1/10, size=4) 
