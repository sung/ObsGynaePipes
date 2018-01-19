library(GenomicFeatures)
library(biomaRt)
library(rtracklayer)
##############################################################
# use this file for the methylation analysis only (e.g. DMR) #
# for DEG, use config/DEG.R                                  #
# ENS_VER defined from the caller                            #
##############################################################
if(TR_PREFIX=="GRCh38"){
	myGenome="hg38" # for goseq (run supportedGenomes())
	#my.local.db<-"~/data/Annotation/GRCh38.p3.TxEnsembl.sqlite" # this will not be used anymore
	#mysql -uanonymous -P3306 -hensembldb.ensembl.org -Densembl_production_82 -e "select distinct(name),biotype_group from biotype where db_type like '%core%' and is_current=1 order by biotype_group,name;" >  ~/data/Annotation/ensembl.grch38.82.biotype.txt
	my.biotype.file<-file.path("~/data/Annotation", paste("ensembl",TR_PREFIX,ENS_VER,"biotype.txt",sep="."))
	#bioType
	#http://rest.ensembl.org/info/biotypes/homo_sapiens
	myBiotype=read.delim(my.biotype.file) # ed by bin/R/Cuffcompare/cuffcompare.R
	my.cgi.file<-'~/data/Annotation/CPGi/cpgi.hg38.bed' # only for the methylation analysis

	#biomart
    #m<-useMart("ensembl")
	#grch38 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset = "hsapiens_gene_ensembl")
    #grch38 = useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
	#myMart = grch38

	my.gene.gtf<-file.path("~/data/genome/Homo_sapiens/Ensembl", paste0(TR_PREFIX,"/Annotation/Genes/Homo_sapiens.",TR_PREFIX,".",ENS_VER,".gtf.gz"))
}else{
	myGenome="hg19" # for goseq (run supportedGenomes())
	my.local.db<-"~/data/Annotation/hg19ensGene.sqlite" #annotation db
	#my.local.db<-"~/data/Annotation/GRCh37.p13.TxEnsembl.sqlite"
	hg19ensGene<- loadDb(my.local.db) # from AnnotationDbi via GenomicFeatures
									  # makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")
									  # see 'bin/R/Etc/make.TxEnsembl.R' for details
	my.cgi.file<-'~/data/Annotation/CPGi/cpgi.hg19.bed' # only for the methylation analysis (hg19-basis)

	#biomart
	grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
	myMart = grch37

	my.gene.gtf<-"~/data/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.75.gtf.local"
}
#gr.ensg<-genes(loadDb(my.local.db)) # from GenomicFeatures
#grl.exons.by.gene <- exonsBy(loadDb(my.local.db),"gene") # isa "GRangeList"
#gr.enst<-transcripts(loadDb(my.local.db))

my.dt.ensg.file<-paste0("~/data/Annotation/Ensembl/",TR_PREFIX,".",ENS_VER,"/dt.ensg.RData")
if(file.exists(my.dt.ensg.file)){
	load(my.dt.ensg.file)
}else{
	library(data.table)
	my.ensg.colnames=c("gene_id","gene_name","gene_biotype")
	my.enst.colnames=c(my.ensg.colnames,"transcript_id")
	#gr.ensg<-rtracklayer::import.gff2(my.gene.gtf, format="gtf", feature.type="gene", colnames=my.ensg.colnames) # isa 'GRanges' # set strands as '*'
	gr.ensg<-rtracklayer::import.gff2(my.gene.gtf, format="gtf", feature.type="gene") # isa 'GRanges' 
	gr.enst<-rtracklayer::import.gff2(my.gene.gtf, format="gtf", feature.type="transcript") # isa 'GRanges' 
	gr.exon<-rtracklayer::import.gff2(my.gene.gtf, format="gtf", feature.type="exon") # isa 'GRanges' 

	dt.ensg<-data.table(data.frame(chromosome_name=seqnames(gr.ensg),
				start_position=start(gr.ensg),
				end_position=end(gr.ensg),
				strand=strand(gr.ensg),
				ensembl_gene_id=mcols(gr.ensg)$gene_id,
				hgnc_symbol=mcols(gr.ensg)$gene_name,
				gene_biotype=mcols(gr.ensg)$gene_biotype
				))

	dt.enst<-data.table(data.frame(chromosome_name=seqnames(gr.enst),
				start_position=start(gr.enst),
				end_position=end(gr.enst),
				strand=strand(gr.enst),
				ensembl_gene_id=mcols(gr.enst)$gene_id,
				ensembl_transcript_id=mcols(gr.enst)$transcript_id,
				hgnc_symbol=mcols(gr.enst)$gene_name
				#transcript_biotype=mcols(gr.enst)$transcript_biotype
				#tsl=mcols(gr.enst)$transcript_support_level
				))
	dt.exon<-data.table(data.frame(chromosome_name=seqnames(gr.exon),
				start_position=start(gr.exon),
				end_position=end(gr.exon),
				strand=strand(gr.exon),
				ensembl_gene_id=mcols(gr.exon)$gene_id,
				ensembl_transcript_id=mcols(gr.exon)$transcript_id,
				ensembl_exon_id=mcols(gr.exon)$exon_id,
				exon_number=mcols(gr.exon)$exon_number,
				hgnc_symbol=mcols(gr.exon)$gene_name
				))
	#fields <- c("ensembl_gene_id","ensembl_transcript_id")
    #dt.ensg.enst <- data.table(getBM(attributes =fields, filters = "ensembl_gene_id", values = dt.ensg$ensembl_gene_id, mart = myMart))

	save(dt.ensg, dt.enst, dt.exon, gr.ensg, gr.enst, gr.exon, file=my.dt.ensg.file)
	load(my.dt.ensg.file)
}
cat("gr.ensg and gt.ensg loaded\n")
