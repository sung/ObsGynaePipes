library(data.table)
##########
# Config #
##########
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM

my.tissue<-list(
	# 1st batch (N=15)
	`Adipose_Tissue`="Adipose - Visceral (Omentum)",
	`Adrenal_Gland`="Adrenal Gland",
	`Blood`="Whole Blood",
	`Blood_Vessel`="Artery - Aorta",
	`Brain`="Brain - Cortex",
	`Colon`="Colon - Sigmoid",
	`Esophagus`="Esophagus - Gastroesophageal Junction",
	`Heart`="Heart - Left Ventricle",
	`Liver`="Liver",
	`Lung`="Lung",
	`Pancreas`="Pancreas",
	`Pituitary`="Pituitary",
	`Small_Intestine`="Small Intestine - Terminal Ileum",
	`Spleen`="Spleen",
	`Thyroid`="Thyroid",
	# 2nd batch (N=4)
	`Breast`="Breast - Mammary Tissue",
	`Muscle`="Muscle - Skeletal",
	`Nerve`="Nerve - Tibial",
	`Skin`="Cells - Transformed fibroblasts",
	`Stomach`="Stomach",

	`Placenta`="Placenta"
)
dt.tissue<-data.table(SMTS=names(my.tissue), SMTSD=as.vector(unlist(my.tissue)))

##################
## Balaton 2015 ##
##################
if(!exists("dt.balaton")){
	dt.balaton<-fread("~/data/X-inactivation/Balaton.2015/Balaton.2015.csv", na.strings="") # GRCh37 base
	dt.balaton[,XCI:=ifelse(grepl("PAR",Y_homology),"PAR", ifelse(is.na(Balaton)|Balaton=="No call",'Unknown',ifelse(grepl("VE$",Balaton),"Variable", ifelse(grepl("E$",Balaton),"Escaped",ifelse(grepl("S$",Balaton),"Inactivated",Balaton)))))]
	dt.balaton[hgnc_symbol=="XIST",XCI:="Escaped"]
	dt.balaton[XCI=="Discordant", XCI:=ifelse(is.na(Carrel), XCI , ifelse(Carrel==1,"Escaped",ifelse(Carrel==0,"Inactivated","Variable")))] # Classify Discondordant case

}

###########################
## Tukiainen 2017 Nature ##
###########################
dt.tuki<-fread("~/data/X-inactivation/Tukiainen.2017/tukiainen_XCI.csv", na.strings="")
dt.tuki$ensembl_gene_id=substr(dt.tuki$ensembl_gene_id,1,15) # ENSG00000182378.8 => ENSG00000182378

makeMeta<-function(){
	dt.gtex.cnt<-fread("zcat /home/ssg29/data/Annotation/GTEx.v6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz", na.strings="")
	dt.gtex.subject<-fread("/home/ssg29/data/Annotation/GTEx.v6p/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt", na.strings="", select=1:2)
	dt.dummy<-fread("/home/ssg29/data/Annotation/GTEx.v6p/GTEx_Data_V6_Annotations_SampleAttributesDS.txt", na.strings="")

	# SMRIN: RIN
	# SMMAPRT: mapping ratio
	# SMEXNCRT: mapping ratio to exon
	# SMMNCPB: coverage per base
	dt.gtex.sample<-dt.dummy[!is.na(SMRIN) & SMMAPRT>0.9 & SMEXNCRT>0.8 & SMMNCPB>20, list(SAMPID,SMTS,SMTSD,SMRIN)]

	# http://stackoverflow.com/questions/18154556/split-text-string-in-a-data-table-columns
	# A-B-C-D => A-B
	if(TRUE){# data.table > v1.9.6+
		dt.gtex.sample[,c("foo","bar"):= tstrsplit(SAMPID, "-", fixed=TRUE, keep=c(1,2))]
		dt.gtex.sample[,SUBJID:=paste(foo,bar,sep="-")]
		dt.gtex.sample[,c("foo","bar"):=list(NULL,NULL)]
	}else{
		dt.gtex.sample$SUBJID= paste(as.character(lapply(strsplit(as.character(dt.gtex.sample$SAMPID), split="-"), "[", 1)),
								as.character(lapply(strsplit(as.character(dt.gtex.sample$SAMPID), split="-"), "[", 2)),
								sep="-")
	}

	dt.gtex.sample<-merge(dt.gtex.sample, dt.gtex.subject)
	dt.gtex.sample[,Sex:=ifelse(GENDER==1,"M","F")]
	dt.gtex.sample[,GENDER:=NULL]
	dt.gtex.sample[,SMTS:=gsub(" ","_",SMTS)]

	# 30 GTEx tissues so far
	#dt.gtex.sample[,.N,"SMTS,SMTSD,Sex"][order(SMTS,SMTSD,Sex)][,.N,"SMTS,SMTSD"][N==1] # one-sex tissues (N=7, Cervix_Uteri, Fallopian_Tube, Ovary, Prostate, Testis, Uterus, Vagina)
	single.sex.SMTSD= dt.gtex.sample[,.N,"SMTS,SMTSD,Sex"][order(SMTS,SMTSD,Sex)][,.N,"SMTS,SMTSD"][N==1,unique(SMTSD)]  # 8 sub-tissues

	# samples less than 20 per tissue per sex
	#dt.gtex.sample[,.N,"SMTS,SMTSD,Sex"][order(SMTS,SMTSD,Sex)][N<20] # Bladder, Kidney
	less.than.20.SMTSD<-dt.gtex.sample[,.N,"SMTS,SMTSD,Sex"][order(SMTS,SMTSD,Sex)][N<20,unique(SMTSD)] # 12 sub-tissues

	# final tissue and sub-tissues 
	#dt.gtex.sample[ !SMTSD %in% c(single.sex.SMTSD,less.than.20.SMTSD),.N,"SMTS,SMTSD,Sex"][order(SMTS,SMTSD,Sex)]
	dt.gtex.sample[ SMTSD %in% my.tissue,.N,"SMTS,SMTSD,Sex"][order(SMTS,SMTSD,Sex)]

	# Meta file
	write.csv(dt.gtex.sample[ SMTSD %in% my.tissue, list(SampleName=SAMPID, HtseqFile=paste0("/home/ssg29/data/Annotation/GTEx.v6p/CNT/",SAMPID,".txt"),Subject=SUBJID,Tissue=SMTS,Sex)], file="~/results/RNA-Seq/GTEx/Meta/meta.ALL.csv", row.names=F)
	#write.csv(dt.gtex.sample[ !SMTSD %in% c(single.sex.SMTSD,less.than.20.SMTSD), list(SampleName=SAMPID, HtseqFile=paste0("/home/ssg29/data/Annotation/GTEx.v6p/CNT/",SAMPID,".txt"),Subject=SUBJID,Tissue=SMTS,Sex)], file="~/results/RNA-Seq/GTEx/Meta/meta.ALL.csv", row.names=F)

	##################
	# gene-count file
	##################
	dt.gtex.cnt[,Name:=substr(Name,1,15)] # ENSG00000223972.4 => ENSG00000223972
	setnames(dt.gtex.cnt,"Name","ENSG")

	#lapply(dt.gtex.sample[ SMTSD %in% my.tissue,SAMPID], function(i) write.table(dt.gtex.cnt[,.(ENSG,get(i))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0("~/data/Annotation/GTEx.v6p/CNT/",i,".txt")) )
	parallel::mclapply(dt.gtex.sample[ SMTSD %in% my.tissue,SAMPID], function(i){write.table(dt.gtex.cnt[,.(ENSG,get(i))], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0("~/data/Annotation/GTEx.v6p/CNT/",i,".txt"))},mc.cores=8 )
} # end of makeMeta

# merge DEG from POPS and GTEx
mergeDEG.DESeq2.1.10<-function(flag='withChrY'){
	my.RData=ifelse(flag=="withChrY", "~/results/RNA-Seq/GTEx/DESeq2.1.10/res.all.tissues.RData", "~/results/RNA-Seq/GTEx/DESeq2.1.10.withoutChrY/res.all.tissues.RData")
	if(!file.exists(my.RData)){
		library(DESeq2) # to load DESeq2

		# patch for FG (25/Jul/2017)
		dl.dummy<-list()
		for(i in c(dt.tissue$SMTS,"Placenta")){
			if(i=="Placenta"){
				file<- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.10/AGA.withoutChrY/toptags_all_deseq.AGA.csv.gz"
			}else{
				file<- file.path("/home/ssg29/results/RNA-Seq/GTEx/DESeq2.1.10.withoutChrY",paste0(i, "/toptags_all_deseq.",i,".csv.gz"))
			}
			dt.dummy<-fread(paste("zcat ",file))
			dt.dummy$Tissue=i
			dl.dummy[[i]]<-dt.dummy[ensembl_gene_id %in% target.ensg[["PT.Novel.Escaped"]]]
		}

		dl.exp.gtex<-list() # initialise 
		for(i in dt.tissue$SMTS){
			cat(paste(i,"...\n"))
			# ~/results/RNA-Seq/GTEx/DESeq2.1.10.withoutChrY/Pancreas/deseq.Pancreas.RData
			# POPS (GRCh37)
			if(i=="Placenta"){
				deseq.RData=ifelse(flag=="withChrY","/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.10/AGA/deseq.AGA.RData","/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.10/AGA.withoutChrY/deseq.AGA.RData")
				#deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.10/AGA.withoutChrY/deseq.AGA.RData"
				#deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.10/AGA/deseq.AGA.RData"
			# GTEx 
			}else{
				deseq.RData=ifelse(flag=="withChrY",file.path("/home/ssg29/results/RNA-Seq/GTEx/DESeq2.1.10",paste0(i, "/deseq.",i,".RData")), file.path("/home/ssg29/results/RNA-Seq/GTEx/DESeq2.1.10.withoutChrY",paste0(i, "/deseq.",i,".RData")))
			}
			if(file.exists(deseq.RData)){
				load(deseq.RData) 
				cat("dds, res, rld?, ddsExonicGene? loaded\n")

				my.res<-res
				my.res$ensembl_gene_id=rownames(my.res)
				dt.exp.gtex=as.data.table(as.data.frame(my.res))
				dt.exp.gtex=merge(dt.exp.gtex, dt.ensg[,.(ensembl_gene_id,hgnc_symbol,gene_biotype,chromosome_name)], by="ensembl_gene_id")
				dt.exp.gtex$Tissue=i
				dl.exp.gtex[[i]]<-dt.exp.gtex
			}else{
				stop(paste(deseq.RData, "not found\n"))
			}
		} # end of SMTS
		#save(dl.exp.gtex, file="~/results/RNA-Seq/GTEx/DESeq2.1.10.withoutChrY/res.all.tissues.RData")
		save(dl.exp.gtex, file="~/results/RNA-Seq/GTEx/DESeq2.1.10/res.all.tissues.RData")
	}else{
		load(my.RData)
	}
	cat("dl.exp.gtex loaded\n")
	return (dl.exp.gtex)
} # end of 'mergeDEG'

# chrY included
mergeDEG.DESeq2.1.18.1<-function(){
    library(DESeq2)
	my.RData="~/results/RNA-Seq/GTEx/DESeq2.1.18.1/GTEx.PT.DESeq2.RData"
	if(!file.exists(my.RData)){
		dl.gtex.deseq<-list()
		dl.Fpkm<-list()
		for(i in dt.tissue$SMTS){
			if(i=="Placenta"){
				file<-"/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.18.1/AGA/filtered_toptags_deseq.all.gene.AGA.csv.gz"
				deseq.RData="/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.18.1/AGA/deseq.AGA.RData"
			}else{
				file<- file.path("/home/ssg29/results/RNA-Seq/GTEx/DESeq2.1.18.1",paste0(i, "/filtered_toptags_deseq.all.gene.",i,".csv.gz"))
				deseq.RData=file.path("/home/ssg29/results/RNA-Seq/GTEx/DESeq2.1.18.1",paste0(i,"/deseq.",i,".RData" ))
			}
			dt.dummy<-fread(paste("zcat ",file))
			dt.dummy$Tissue=i
			dl.gtex.deseq[[i]]<-dt.dummy
			if(file.exists(deseq.RData)){
				load(deseq.RData) 
				cat(paste("dds and res loaded for",i,"\n"))
                ddsFpkm <-fpkm(dds) # isa 'matrix'
	            dt.Fpkm<-data.table(`Tissue`=i, `ensembl_gene_id`=rownames(ddsFpkm), ddsFpkm)
			    dl.Fpkm[[i]]<-dt.Fpkm
			}else{
				stop(paste(deseq.RData, "not found\n"))
			}
		}
		save(dl.gtex.deseq, dl.Fpkm, file=my.RData)
	    cat("dl.gtex.deseq and dl.Fpkm saved\n")
    }else{
		load(my.RData)
    }    
	cat("dl.gtex.deseq and dl.Fpkm loaded\n")
} # end of 'mergeDEG'

# it loads 'dl.Fpkm' if not
getGTExFpkm<-function(my.target.ensg){
    if (length(my.target.ensg)==0) stop ("target ENSG<1")
    if(!exists("dl.Fpkm")){mergeDEG.DESeq2.1.18.1()}
    dt.fpkm.tissue<-rbindlist(
                        lapply(dl.Fpkm, function(i)
                                melt.data.table(i[ensembl_gene_id %in% my.target.ensg], 
                                                id.vars=c("Tissue","ensembl_gene_id"), 
                                                variable.name="SampleID", 
                                                value.name="FPKM"
                                                ) # wide (col-based)to long (row-based)
                    )
            )
    return(dt.fpkm.tissue)
}

#dl.exp.gtex=mergeDEG()
my.RData="~/results/chrX/RData/deg.XCI.all.tissues.RData"
if(file.exists(my.RData)){
	load(my.RData)
	cat("dl.exp.gtex, dt.chrX.gene, dt.exp.all.chrX, dt.deg.all.chrX, dt.balaton.all.deg loaded\n")

	if(FALSE){
	write.csv(merge(dt.deg.all.chrX, dt.tuki[,.(ensembl_gene_id,reported_XCI,sex_biased_gtex,XCI_across_tissues,XCI_single_cell,new_escape)],all.x=T)[order(Tissue,-log2FoldChange)], file=gzfile("~/results/chrX/CSV/female.biased.genes.placenta.and.19.GTEx.tissues.GRCh37.csv.gz"), row.names=F, quote=F)
	}

}else{
	dt.exp.all<-rbindlist(dl.exp.gtex)

	# No. of DEG by chr per tissue
	dcast(dt.exp.all[chromosome_name %in%  my.chr.order & padj<0.01,.N,"Tissue,chromosome_name"], Tissue ~ chromosome_name, value.var="N") # my.chr.order defined in 'config/graphic.R'
	dcast(dt.exp.all[chromosome_name %in%  my.chr.order & padj<0.01,.N,"Tissue,chromosome_name"], chromosome_name ~ Tissue, value.var="N")

	# ChrX only & no 'Breast' tissue
	dt.exp.all.chrX<-dt.exp.all[chromosome_name=="X" & Tissue!="Breast"]
	#dt.exp.all.chrX<-dt.exp.all[chromosome_name=="X"]
	# adjust p-value based on chrX only
	dt.exp.all.chrX[,`new.padj` := p.adjust(pvalue, method="BH"),by="Tissue"]

	#get gene-name
	#fields <- c("ensembl_gene_id", "hgnc_symbol", "start_position", "gene_biotype", "description")
	#dt.chrX.gene<-as.data.table(getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.exp.all.chrX[, unique(ensembl_gene_id)], mart = myMart)) # is a data.table
	dt.chrX.gene<-dt.ensg[ensembl_gene_id %in% dt.exp.all.chrX[, unique(ensembl_gene_id)]]

	# new.padj<0.01 
	# abs(log2FoldChange)>log2(1.1)  at least 10% difference 
	# baseMean>10 : at least 10 reads
	dt.deg.all.chrX<-merge(dt.exp.all.chrX[new.padj<0.01 & abs(log2FoldChange)>log2(1.1) & baseMean>10], dt.chrX.gene, by="ensembl_gene_id")

	# Placenta-specific DEG
	#dt.deg.all.chrX[Tissue=="Placenta" & !ensembl_gene_id %in% dt.deg.all.chrX[Tissue!="Placenta",unique(ensembl_gene_id)]]

	# No. of nw.padj chrX DEG by Tissue
	#data.table(dcast(dt.deg.all.chrX[,.N,.(Tissue,`F>M`=log2FoldChange>0)], Tissue~`F>M`))

	# Annotate XCI status defined by Balaton et al.
	#dt.balaton.all.deg<-merge(dt.deg.all.chrX, dt.balaton, by="hgnc_symbol", all.x=TRUE)[,list(Tissue,hgnc_symbol,ensembl_gene_id,gene_biotype,baseMean,log2FoldChange,padj,new.padj,Balaton,Carrel,other_literature,Y_homology,description)]
	dt.balaton.all.deg<-merge(dt.deg.all.chrX, dt.balaton, by="hgnc_symbol", all.x=TRUE)

	# No. of genes by XCI status per tissue
	#dcast(dt.balaton.all.deg[,.N,"Tissue,Expression?,XCI"], Tissue + `Expression?` ~ XCI)

	save(dl.exp.gtex, dt.chrX.gene, dt.exp.all.chrX, dt.deg.all.chrX, dt.balaton.all.deg, file=my.RData)
}

#############################################
## Used by "bin/R/chrX/tissue_comparison.R"
## it uses 'dt.chrX.gene'
## 'target.ensg'
#############################################
plotHeatMaplog2FC<-function(my.target, is.order.position=T, is.sort.col=T, is.sort.row=F, is.save=F, my.fontsize_row=12){
	my.dt.target<-dt.exp.all.chrX[ensembl_gene_id %in% target.ensg[[my.target]]]

	dt.foo<-dcast(my.dt.target, ensembl_gene_id~Tissue, value.var="log2FoldChange")
	if(is.order.position){
		dt.bar<-data.table(merge(dt.foo, dt.chrX.gene[,.(ensembl_gene_id, hgnc_symbol, start_position)]))[order(start_position)]
	}else{
		dt.bar<-data.table(merge(dt.foo, dt.chrX.gene[,.(ensembl_gene_id, hgnc_symbol, start_position)]))
	}
	dt.bar[,gene_name:=ifelse(hgnc_symbol=="",as.character(ensembl_gene_id),hgnc_symbol)] # replace empy hgnc_symbol with ensembl_gene_id

	my.end<-ncol(dt.bar)-3
	merged.log2fc<-as.matrix(dt.bar[,2:my.end])
	rownames(merged.log2fc)<-dt.bar$gene_name

	my.min.scale=floor(min(as.vector(merged.log2fc),na.rm=T))
	my.max.scale=ceiling(max(as.vector(merged.log2fc),na.rm=T))
	# breaks: a sequence of numbers that covers the range of values in mat and is one element longer than color vector
	if(my.min.scale < -1){
		if(my.max.scale>1){
			my.break=c(my.min.scale,seq(-1,1,by=0.1),my.max.scale)
		}else{
			my.break=c(my.min.scale,seq(-1,1,by=0.1))
		}
	}else{
		# assume min.scale < 1
		if(my.max.scale>1){
			my.break=c(seq(my.min.scale, 1, by=0.1), my.max.scale)
		}else{
			my.break=c(seq(my.min.scale, 1, by=0.1))
		}
	}
	#table(cut(as.vector(merged.log2fc), breaks=my.break, include.lowest =T, right=F))

	my.num.less.than.zero=sum(as.numeric(my.break<0))
	# colors: one less than breaks
	my.col.blue=colorRampPalette(colors = c("blue","#f2f2ff"))(my.num.less.than.zero)
	my.col.red=colorRampPalette(colors = c("#f7eaea","red"))(length(my.break)-my.num.less.than.zero-1)
	my.col.scale=c(my.col.blue, my.col.red)

	if(is.save){
		my.filename=ifelse(is.sort.row,paste0("~/results/chrX/Figures/pheatmap.",my.target,".gene.clustered.tiff"),paste0("~/results/chrX/Figures/pheatmap.",my.target,".tiff"))
		pheatmap::pheatmap(merged.log2fc, cluster_rows=is.sort.row, cluster_cols=is.sort.col, breaks=my.break, color =my.col.scale, width=10, height=15, fontsize_col=12, fontsize_row=my.fontsize_row, filename=my.filename) 
	}else{
		pheatmap::pheatmap(merged.log2fc, cluster_rows=is.sort.row, cluster_cols=is.sort.col, breaks=my.break, color =my.col.scale) 
	}
	#pheatmap::pheatmap(merged.log2fc, cluster_rows=F, cluster_cols=T, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}
