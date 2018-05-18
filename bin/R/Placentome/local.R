#library(data.table)

#############
## Samples ##
#############
outLiers=list(`Failed`=c("91"), `IGFBP1`=c("84C","88","80C"), `SGA.PET`=c("77","71","79","19P","02P","72P"))
dt.samples.all<-fread("/home/ssg29/results/RNA-Seq/Placentome/Meta/reconstruct.stringtie.csv", na.strings="") # n=325 for all samples
cat("Reading input file...\n")
dt.read.cnt<-fread(file.path("~/results/RNA-Seq/Placentome/Meta", paste0("FG.JD.SE125.",TR_PREFIX,".mapping.txt")), header=F)
setnames(dt.read.cnt,c("Library","BarCode","Category","Read"))

## FG ^ JD-BR
#merge(dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Source=="FG",.(SampleName,CRN,Source)], dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Source=="JD-BR",.(SampleName,CRN,Source)], by="CRN")
## JD ^ JD-BR
#merge(dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Source=="JD",.(SampleName,CRN,Source)], dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Source=="JD-BR",.(SampleName,CRN,Source)], by="CRN")
## FG ^ JD
#merge(dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Source=="FG",.(SampleName,CRN,Source)], dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Source=="JD",.(SampleName,CRN,Source)], by="CRN")

#########
## CTL ##
# n=152 #
#########
# FG & JD 
foo=merge(dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="FG",.(CRN,Source)], dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="JD",.(CRN,Source)],by="CRN")$CRN
# FG & BR
bar=merge(dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="FG",.(CRN,Source)], dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="JD-BR",.(CRN,Source)],by="CRN")$CRN
# JD & BR
tar=merge(dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="JD",.(CRN,Source)], dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="JD-BR",.(CRN,Source)],by="CRN")$CRN
# 152 unique samples
dt.ctl<-rbind(dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="JD-BR"], # breech samples only (n=24, 24 unique)
			dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="FG" & ! CRN %in% bar], # FG samples excluding BR (n=57, 55 unique)
			dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==0 & Source=="JD" & ! CRN %in% c(foo,tar)] # JD samples excluding FG & BR (n=92, 73 unique)
			)
#write.table(dt.ctl[,StringTieFile],file="~/results/RNA-Seq/Placentome/Meta/CTL.stringtie.GRCh38.82.gtf.list", quote=F,row.names=F,col.names=F)

#########
## SGA ##
# n=52  #
#########
#dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==1 & Source=="FG"] # n=55
dt.sga<-dt.samples.all[!SampleName %in% unname(unlist(outLiers)) & Condition==1 & Source=="FG"] # n=52
#write.table(dt.samples.all[!SampleName %in% unname(unlist(outLiers)) & Condition==1 & Source=="FG",StringTieFile], file="~/results/RNA-Seq/Placentome/Meta/SGA.stringtie.GRCh38.82.gtf.list", quote=F,row.names=F,col.names=F)

#########
## PET ##
# n=91  #
#########
#dt.samples.all[!SampleName %in% outLiers[["IGFBP1"]] & Condition==1 & Source=="JD"] # n=94
dt.pet<-dt.samples.all[!SampleName %in% unname(unlist(outLiers)) & Condition==1 & Source=="JD"] # n=91
#write.table(dt.samples.all[!SampleName %in% unname(unlist(outLiers)) & Condition==1 & Source=="JD",StringTieFile], file="~/results/RNA-Seq/Placentome/Meta/PET.stringtie.GRCh38.82.gtf.list", quote=F,row.names=F,col.names=F)

dt.pops<-rbind(dt.ctl[,`Cohort`:="CTL"], dt.sga[,`Cohort`:="SGA"], dt.pet[,`Cohort`:="PET"])
dt.pops$Cohort=factor(dt.pops$Cohort,levels=c("CTL","SGA","PET"))

#####################
## Read Count Data ##
## see ~/results/RNA-Seq/Placentome/Meta/README
#####################
dt.pops.read.cnt<-merge(dt.pops, dt.read.cnt, by=c("Library","BarCode"))
#dt.sample.cohort=unique(dt.pops.read.cnt[,.(Library,BarCode,Cohort)])[,list(`Sample`=paste(Library,BarCode,sep="."),Cohort)] # use dt.pops above


prep.cuff<-function(){
	# ~/results/RNA-Seq/Placentome/Meta/CTL.stringtie.GRCh38.82.gtf.list
	input.gtf.list=file.path(top.dir, paste0("Meta/",myCohort,".stringtie.",TR_PREFIX,".",ENS_VER,".gtf.list"))  
	my.tracking.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".",ENS_VER,".cuffcompare.tracking")) # per transcript (TCON_)
	my.gtf.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".",ENS_VER,".cuffcompare.combined.gtf")) 
	#my.tracking.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".",ENS_VER,".cuffcompare.",myTarget,".tracking")) # per transcript (TCON_)
	#my.gtf.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".",ENS_VER,".cuffcompare.",myTarget,".gtf")) 

	##########################
	## Transfrags Definition #
	##########################
	cat("creating Transfrag...\n")
	my.cmd<-paste("awk 'BEGIN{FS=\"\t\";OFS=\",\"}{hgnc=\"NA\";enst=\"NA\";if($3!=\"-\"){split($3,ref,\"|\"); hgnc=ref[1]; enst=ref[2]}print $1,$2,hgnc,enst,$4}' ", my.tracking.file)
	dummy<-system(my.cmd, intern=T)
	cuff.transfrag<-t(simplify2array(strsplit(dummy, ","))) # row: TCON, columns: TCON,XLOC,HGNC,ENST,CODE 
	cuff.transfrag[cuff.transfrag[,3]=="NA",3]<-NA # convert character "NA" to NA
	cuff.transfrag[cuff.transfrag[,4]=="NA",4]<-NA # convert character "NA" to NA
	rownames(cuff.transfrag)<-cuff.transfrag[,1]
	colnames(cuff.transfrag)<-c("tcon","xloc","hgnc_symbol","ensembl_transcript_id","class_code")   # isa 'matrix'

	########################
	## Ensembl Annotation ##
	########################
	cat("creating anno...\n")
	fields <- c("ensembl_transcript_id", "transcript_biotype")
	anno <- getBM(attributes =fields, filters = fields[1], values = unique(cuff.transfrag[,"ensembl_transcript_id"]), mart = myMart)
	rownames(anno)<-anno$ensembl_transcript_id
																									# 'hgnc_symbol' not really! (read the ensembl genebuild)
	##########
	## FPKM ##
	##########
	cat("creating Fpkm...\n")
	my.cmd<-paste("awk -f ~/Pipelines/bin/parse.cuffcompare.fpkm.awk ", my.tracking.file)
	dummy<-system(my.cmd, intern=T)
	dummy<-t(simplify2array(strsplit(dummy, ","))) # row: TCON, columns: samples

	cuff.fpkm<-matrix(as.numeric(dummy[,2:ncol(dummy)]), ncol=ncol(dummy)-1) # isa 'matrix', NAs introduced by coercion
	rownames(cuff.fpkm)<-dummy[,1]
	colnames(cuff.fpkm)<-unname(sapply(read.table(input.gtf.list, header=FALSE, stringsAsFactors=FALSE)[,"V1"],function(i) paste(unlist(strsplit(unlist(strsplit(i, "/"))[8], ".", fixed=T))[1], unlist(strsplit(i, "/"))[7], sep=".") ))
	# isa matrix

	#######################
	## Depth of Coverage ##
	#################e#####
	cat("creating Coverage...\n")
	my.cmd<-paste("awk -f ~/Pipelines/bin/parse.cuffcompare.cov.awk ", my.tracking.file)
	dummy<-system(my.cmd, intern=T)
	dummy<-t(simplify2array(strsplit(dummy, ","))) # row: TCON, columns: samples

	cuff.cov<-matrix(as.numeric(dummy[,2:ncol(dummy)]), ncol=ncol(dummy)-1) # isa 'matrix', NAs introduced by coercion
	rownames(cuff.cov)<-dummy[,1]
	colnames(cuff.cov)<-unname(sapply(read.table(input.gtf.list, header=FALSE, stringsAsFactors=FALSE)[,"V1"],function(i) paste(unlist(strsplit(unlist(strsplit(i, "/"))[8], ".", fixed=T))[1], unlist(strsplit(i, "/"))[7], sep=".") ))
	# isa matrix

	##########################
	## Length of Transcript ##
	##########################
	cat("creating Length...\n")
	my.cmd<-paste("awk -f ~/Pipelines/bin/parse.cuffcompare.len.awk ", my.tracking.file)
	dummy<-system(my.cmd, intern=T)
	dummy<-t(simplify2array(strsplit(dummy, ",")))

	cuff.len<-matrix(as.numeric(dummy[,2:ncol(dummy)]), ncol=ncol(dummy)-1) # isa 'matrix', NAs introduced by coercion
	rownames(cuff.len)<-dummy[,1]
	colnames(cuff.len)<-unname(sapply(read.table(input.gtf.list, header=FALSE, stringsAsFactors=FALSE)[,"V1"],function(i) paste(unlist(strsplit(unlist(strsplit(i, "/"))[8], ".", fixed=T))[1], unlist(strsplit(i, "/"))[7], sep=".") ))
	# isa matrix

	######################
	## GFF of TransFrag ##
	######################
	cat("loading GFF...\n")
	cuff.gff=rtracklayer::import.gff2(my.gtf.file)
	#tcons.gr.list<-split(cuff.gff, mcols(cuff.gff)[["transcript_id"]])
	#tcons.exon.cnt=sapply(tcons.gr.list, length) # takes too long

	##########################
	# Exon cnt per TransFrag #
	##########################
	cat("creating tcons.exon.cnt...\n")
	my.cmd<-paste("awk -f ~/Pipelines/bin/parse.cuffcompare.exon.cnt.awk ", my.gtf.file)
	dummy<-system(my.cmd, intern=T)
	dummy<-t(simplify2array(strsplit(dummy, " ")))
	tcons.exon.cnt<-data.frame(row.names=dummy[,1], exon.cnt=as.numeric(dummy[,2]))

	######################
	## Update Transfrags #
	######################
	cat("updating Transfrags...\n")
	cuff.transfrag<-as.data.frame(cuff.transfrag, stringsAsFactors=FALSE)
	## apply(x, 1) : for each row (transfrag)
	## quantile(x) assume at least 4 elements
	cuff.transfrag[["fpkm.iqr"]]=apply(cuff.fpkm, 1, function(i) ifelse(sum(!is.na(i))>=3,mean(i[i>=quantile(i,0.25,na.rm=T) & i<=quantile(i,0.75,na.rm=T)],na.rm=T),mean(i,na.rm=T)))
	cuff.transfrag[["cov.iqr"]]=apply(cuff.cov, 1, function(i) ifelse(sum(!is.na(i))>=3,mean(i[i>=quantile(i,0.25,na.rm=T) & i<=quantile(i,0.75,na.rm=T)],na.rm=T),mean(i,na.rm=T)))
	cuff.transfrag[["len.iqr"]]=apply(cuff.len, 1, function(i) ifelse(sum(!is.na(i))>=3,mean(i[i>=quantile(i,0.25,na.rm=T) & i<=quantile(i,0.75,na.rm=T)],na.rm=T),mean(i,na.rm=T)))
	cuff.transfrag[["exon.cnt"]]=tcons.exon.cnt[rownames(cuff.transfrag),]
	cuff.transfrag[["evi.ratio"]]=rowSums( !is.na( cuff.fpkm ) ) / ncol(cuff.fpkm)
	cuff.transfrag[["transcript_biotype"]]=ifelse(is.na(cuff.transfrag[["ensembl_transcript_id"]]), NA, anno[cuff.transfrag[["ensembl_transcript_id"]],"transcript_biotype"])
	#cuff.transfrag[["transcript_biotype"]]=anno[cuff.transfrag[["ensembl_transcript_id"]],"transcript_biotype"]

	################
	## Save RData ##
	################
	cat("Saving RData...\n")
	assign(myProject, list(`Transfrag`=as.data.table(cuff.transfrag),`Fpkm`=cuff.fpkm,`Coverage`=cuff.cov,`Length`=cuff.len, `Gff`=cuff.gff))
	save(list=myProject, file=cuffRData)
	cat("A list named ", myProject, "saved to ",cuffRData,"\n")
	#return(myProject)
}

do.go.cuff<-function(){
	cat("doing GO analysis...\n")
	library(goseq) # for goseq
	library(GO.db) # for ,GOBPPARENTS
	library(GOstats) # for GOGraph

	cat("fetching ensembl_gene_id...\n")
	#select<-dt.frag$fpkm.iqr >= minFpkm & dt.frag$class_code=="=" & dt.frag$transcript_biotype=="protein_coding"
	dummy<-sapply(c(">=95%",">=5%"), function(i) as.character(unique(dt.frag[fpkm.iqr >= minFpkm & class_code=="=" & transcript_biotype=="protein_coding" & evi.ratio >= sampleFreq[i],ensembl_transcript_id]))) # isa 'list'
	dummy<-lapply(dummy, function(i) unique(getBM(attributes ="ensembl_gene_id", filters ="ensembl_transcript_id", values = i, mart = myMart)[[1]])) #isa 'list'
	myGenes<-as.integer(dummy[[">=5%"]] %in% dummy[[">=95%"]])
	names(myGenes)<-dummy[[">=5%"]]

	if(myGenome=="hg19"){
		pwf=goseq::nullp(myGenes,myGenome,"ensGene") # isa data.frame  (check available ID via supportedGeneIDs())
		GO.wall=goseq::goseq(pwf,myGenome,"ensGene") # isa data.frame
		#GO.samp=goseq::goseq(pwf,myGenome,"ensGene",test.cats=c("GO:MF"),method="Sampling",repcnt=1000)
	}else{
		cat("calulating gene length...\n")
		grl.exons.by.gene <- exonsBy(loadDb(my.local.db),"gene") # isa "GRangeList"
		gene_lengths=sapply(names(myGenes),function(x) if(x %in% names(grl.exons.by.gene)){sum(width(reduce(grl.exons.by.gene[[x]])))}else{NA})
		pwf=goseq::nullp(myGenes,bias.data=gene_lengths)

		#cat("fetching go_id...\n")
		#go_map=sapply(names(myGenes),function(x) getBM(attributes ="go_id", filters ="ensembl_gene_id", values = x, mart = myMart)[[1]] ) # isa 'list'
		#GO.wall=goseq::goseq(pwf,myGenome,"ensGene", gene2cat=go_map) # isa data.frame
		GO.wall=goseq::goseq(pwf,myGenome,"ensGene") # isa data.frame
	}
	enriched.GO.wall=list('over'=GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.001,], # isa data.frame
						'under'=GO.wall[p.adjust(GO.wall$under_represented_pvalue, method="BH")<.001,] # isa data.frame
					)

	# Find the level from the root (depth-first or shortest) given a GO term
	for(j in names(enriched.GO.wall)){
		enriched.GO.wall[[j]]$depth<-as.numeric(sapply(enriched.GO.wall[[j]]$category, function(i) {
					if(Ontology(i)=="BP"){
						graph::DFS(GOGraph(i,GOBPPARENTS), "all")[i]
					}else if(Ontology(i)=="MF"){
						graph::DFS(GOGraph(i,GOMFPARENTS), "all")[i]
					}else{
						graph::DFS(GOGraph(i,GOCCPARENTS), "all")[i]
					}
				} # end of function(i)
			)# end of sapply
		) # end of as.numeric
	}
	#DFS(GOGraph("GO:0044764",GOBPPARENTS), "all")["GO:0044764"]
	#plot(GOGraph("GO:0044764",GOBPPARENTS))

	# list to data.frame
	cat("writing GO analysis results...\n")
	write.csv(enriched.GO.wall[["over"]][enriched.GO.wall[["over"]]$depth-1>=5,], file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort, ".goseq.coding.95.vs.5.over.csv")))
	write.csv(enriched.GO.wall[["under"]][enriched.GO.wall[["under"]]$depth-1>=5,], file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort, ".goseq.coding.95.vs.5.under.csv")))

	############################
	## GO terms to Word Cloud ##
	############################
	cat("making WordCloud...\n")
	library(tm)
	library(wordcloud)
	library(SnowballC)

	terms<-as.character(unlist(sapply(enriched.GO.wall[["over"]][enriched.GO.wall[["over"]]$depth-1>=5,"term"], function(i) strsplit(i, " ")))) # isa character vector
	terms<-tm::Corpus(VectorSource(terms)) # isa 'VCorpus'

	terms<- tm_map(terms, content_transformer(tolower))
	terms<- tm_map(terms, removePunctuation)
	terms<- tm_map(terms, removeNumbers)
	terms<- tm_map(terms, removeWords, stopwords("english"))
	terms<- tm_map(terms, stripWhitespace)
	#terms<- tm_map(terms, stemDocument) #install.packages("SnowballC")
	wordcloud(terms, scale=c(5,0.5), max.words=100, random.order=FALSE, rot.per=0.35, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2")) # brewer.pal from RColorBrewer
}

# Suppl information from: https://academic.oup.com/bib/article/18/2/205/2562739
# Function require a vector with expression of one gene in different tissues.
# If expression for one tissue is not known, gene specificity for this gene is NA
# Minimum 2 tissues
fTau <- function(x){
    if(all(!is.na(x))){
        if(min(x, na.rm=TRUE) >= 0){
            if(max(x)!=0){
                x <- (1-(x/max(x)))
                res <- sum(x, na.rm=TRUE)
                res <- res/(length(x)-1)
            }else{
                res <- 0
            }
 		}else{
            res <- NA
            #print("Expression values have to be positive!")
 		} 
 	}else{
        res <- NA
        #print("No data for this gene avalable.")
    } 
    return(res)
}

