#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(DESeq2)
library(genefilter) # for rowVar
library(gplots) # for heatmap.2
library(BiocParallel)

source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
source("~/Pipelines/config/graphic.R")
NP=8 # No. of worker process for DESeq
register(MulticoreParam(NP))

#################
## PDN Samples ##
#################
#myProject="PDN"
#sampleType="YS" # "Placenta": Human 1st trimester Placenta,  "YS" # Human Yolk Sac (not for DEG, but for clustering by bin/R/TopExp/find.commonly.expressed.genes.R)

##############################
## DEGs from Boys and Girls ##
##############################
#myProject=paste("Boy.Girl.FG.JD",TR_PREFIX,sep=".") # Boy.Girl
#myProject=paste("Boy.Girl.exon",TR_PREFIX,sep=".")
#myProject=paste("Boy.Girl.FG.JD.mature.miRNA",TR_PREFIX,sep=".")
#myProject=paste("Boy.Girl.FG.JD.pre.miRNA",TR_PREFIX,sep=".")
#myProject=paste("Boy.Girl.FG.JD.piRNA",TR_PREFIX,sep=".")
#sampleType='ALL' # all samples from JD and FG (including chrY for Salmon GRCh38.82)
#sampleType='AGA' # healthy samples from JD and FG (including chrY for Salmon GRCh38.82)
#sampleType='SGA'
#sampleType='PET'
#sampleType='CASE'
#sampleType='BR' # chrY genes included for this type to compare with plasma samples
#sampleType='BR.Salmon' # chrY genes included for this type to compare with plasma samples

###################################
## CSMD1 exon count from RoadMap ##
###################################
#myProject="RoadMap.exon.GRCh37"
#sampleType="Schultz"
#sampleType="Consortium"
#sampleType="FetalEmbryo"
#sampleType="ALL"

##############################
## DEGs from PET vs. Normal ##
##############################
#myProject=paste("PET.RPT",TR_PREFIX,sep=".") # Re-constructed Placenta Transcrptome
#myProject=paste("PET.RPT.TX",TR_PREFIX,sep=".") # Re-constructed Placenta Transcrptome
#myProject=paste("PET",TR_PREFIX,"salmon",sep=".")  # DESeq ==1.10 and 1.18.1
#myProject=paste("PET",TR_PREFIX,"featureCount",sep=".")  # DESeq2 ==1.18.1
#sampleType='ALL'				   # All samples (i.e n=110 of FG samples)
#myProject="PET.piRNA"	

############################
# DEGs between SGA and AGA #
############################
#myProject=paste("SGA.AGA.SE125",TR_PREFIX,"salmon",sep=".")		# totalRNA-Seq (SGA.AGA.total): GRCh37|GRCh38
#myProject=paste("SGA.AGA.SE125",TR_PREFIX,"HTSeq",sep=".")		# totalRNA-Seq (SGA.AGA.total): GRCh37|GRCh38
#myProject=paste("SGA.AGA.SE125",TR_PREFIX,"featureCount",sep=".")		# totalRNA-Seq (SGA.AGA.total): GRCh37|GRCh38
#myProject=paste("SGA.AGA.miRNA",TR_PREFIX,sep=".")		# miRNA from the small RNA-Seq (v1: GRCh37, v2: GRCh38)
#myProject=paste("SGA.AGA.piRNA",TR_PREFIX,sep=".")		# piRNA from the small RNA-Seq (v1: GRCh37, v2: GRCh38)
#myProject=paste("SGA.AGA.RPT",TR_PREFIX,sep=".") # Re-constructed Placenta Transcrptome
#myProject=paste("SGA.AGA.RPT.TX",TR_PREFIX,sep=".") # Re-constructed Placenta Transcrptome (tcons_xx based)

#sampleType='ALL'				   # All samples (i.e n=110 of FG samples)
#sampleType='ALL.PHEN'			   # 11 SGA vs. 11 matched AGA where all measurements (pappa, gvel, doppler, pet are avaiable)
#sampleType='LOW.PAPPA'			   # 8 low PAPP-A SGA vs. 8 matched AGA of which 4 PAPPA are unknown <NA>
#sampleType='LOW.PAPPA.UNIQ'	   # 3 unique low PAPP-A SGA vs. 3 matched AGA
#sampleType='LOW.PAPPA.EXC'	       # 7 exclusive low PAPPA-A SGA vs. 7 matched AGA
#sampleType='NORMAL.PAPPA'		   # 14 normal PAPP-A SGA vs. 14 matched AGA
#sampleType='PAPPA.SGA'			   # 8 low PAPP-A SGA vs. 14 unmatched AGA
#sampleType='LOW.GVEL'			   # 13 Low GVel SGA vs. 13 matched AGA 
#sampleType='LOW.GVEL.UNIQ'		   # 7 unique Low GVel SGA vs. 7 matched AGA 
#sampleType='LOW.GVEL.EXC'		   # 7 exclusive Low GVel 
#sampleType='NORMAL.GVEL'		   # 40 Normal GVel SGA vs. 40 matched AGA 
#sampleType='HIGH.UAD'			   # 12 high (abnomral) Uterine Artery Doppler (UAD) SGA vs. 12 matched AGA of which 1 high UAD 
#sampleType='HIGH.UAD.UNIQ'		   # 2 unique high (abnomral) Uterine Artery Doppler (UAD) SGA vs. 2 matched AGA
#sampleType='HIGH.UAD.EXC'		   # 6 exclusive HIGH.UAD
#sampleType='NORMAL.UAD'		   # 40 nomral UAD SGA vs. 40 matched AGA
#sampleType='HIGH.UBD'			   # 11 high (abnormal) Umbilical Doppler (UBD) vs. 11 matched AGA of which 1 high UBD (SLX-9169 D711_D502  P53)
#sampleType='HIGH.UBD.UNIQ'		   # 5 unique high (abnormal) Umbilical Doppler (UBD) vs. 5 matched AGA
#sampleType='HIGH.UBD.EXC'		   # 8 exclusive HIGH.UBD
#sampleType='NORMAL.UBD'		   # 42 normal UBD SGA vs. 42 matched AGA
#sampleType='HIGH.DUS'			   # high Doppler Ultra Sound (DUS, UBD+UAD) which is abnomral 
#sampleType='NORMAL.DUS'		   # normal DUS 
#sampleType='HT'				   # 11 HT SGA and 11 matched AGA (HT from 7 PET + 3 GH + 1 ambigous)
#sampleType='HT.UNIQ'			   # 3 unique HT SGA and 3 matched AGA (HT from 7 PET + 3 GH + 1 ambigous)
#sampleType='HT.EXC'			   # 8 exclusive HT
#sampleType='NHT'				   # 45 NHT SGA vs. 45 matched AGA (excluding 11 pairs from PET + Hypertension + ambigous)
#sampleType='PET'				   # 6 PET SGA vs. 6 matched AGA
#sampleType='PET.SGA'			   # 6 PET SGA vs. 45 unmatched NHT SGA
#sampleType='OXBS'				   # MJ 4 low PAPP-A SGA vs. 4 matched AGA 
#sampleType='Boys'				   # FG 44 boys only 
#sampleType='Girls'				   # FG 56 girs only 
#sampleType='LOW.PAPPA.Boys'	# 4 low PAPP-A Boy SGA vs. 4 matched Boy AGA
#sampleType='LOW.PAPPA.Girls'	# 4 low PAPP-A Girl SGA vs. 4 matched Girl AGA

						   # by bin/R/TopExp/find.commonly.expressed.genes.R)
#sampleType='AGA'          # healty babies only (not for DEG, but for clustering)
#sampleType='SGA'          # small babies only  (not for DEG, but for clustering)

####################
## Irving Samples ##
## GRCh38         ##
####################
#myProject='Retosiban'
#myProject='Retosiban.miRNA'
#myProject='Retosiban.piRNA'
#myProject='Retosiban.GRCh38.2nd'
#sampleType='R8' #ALL, Stretch, R8, R6, Dosage
#myProject='Trophoblast.IA' # GRCh38.90 salmon
#sampleType='F.siRNA' #M.DFMO F.DFMO M.siRNA F.siRNA

####################################
## FG Plasma Samples              ##
## SLX-7630: SE50, 4 samples only ##
## SLX-9342: SE50, PE75, PE150    ##
## SLX-9345: SE50, PE75           ##
## All based on GRCh38            ##
## based on DESeq2 < 1.16         ##
####################################
#myProject="Plasma.2017"    # Plasma.2017 or Plasma.2016
#sampleType='ALL'    # PE75+PE150 Plasma.2017 (F-M)
#sampleType='ALL.Salmon'    # Plasma.2017 PE75 GRCh38.82 
#sampleType='GA.36'          # Plasma.2017 GA 12, 20, 28, 36 wk
#sampleType='PLS.PT'       # Plasma.2016 12WK vs. Placenta
#sampleType='PLS12WK.PT'   # Plasma.2016 12WK vs. Placenta
#sampleType='PLS36WK.PT'   # Plasma.2016 36WK vs. Placenta
#sampleType='PLS36WK.12WK' # Plasma.2016

###########################################
## Illumina Plasma Breach Samples         #
## N=83 samples (83 BAM files received )  # 
## from 22 patients (11 males 11 females) #
## BAM files mapped against hg19          #
## featureCount based on GRCh37.82        #
## Salmon based on GRCh38.82              #
## based on DESeq2 1.18                   #
###########################################
#myProject="Plasma.ILM.2017" #
#sampleType='ALL'      # (F-M)
#sampleType='ALL.Salmon'    # Plasma.2017 PE75 GRCh38.82 
#sampleType='GA.36'          # Plasma.2017 GA 12, 20, 28, 36 wk

###########################################
## Illumina Placenta Breach Samples       #
## N=19 samples (19 BAM files received )  # 
## from 22 patients (11 males 11 females) #
## but BR18 missing from library (female) #
## BR09 missing from BAM (male)           #
## BR14 missing from BAM (female)         #
## BR17,BR08 NOT SENT (n=2)               #
## BAM files mapped against hg19          #
## based on DESeq2 1.18                   #
###########################################
myProject=paste("Boy.Girl.ILM",TR_PREFIX,sep=".")
#sampleType="BR" # F-M
sampleType="BR.Salmon" # F-M

#############################
## SML Cholestatis samples ##
## based on DESeq2 1.18.1  ##
#############################
#myProject="Cholestasis.2017"
#sampleType='ALL' 

##############################
## GTEx V6                  ##
## GENCODE v19 (Ensembl 75) ##
## GRCh37                   ##
## bin/R/GTEx/local.R       ##
##############################
#myProject="GTEx"
#sampleType='Adipose_Tissue'      # Boy vs. Girl 
#sampleType='Adrenal_Gland'      # Boy vs. Girl 
#sampleType='Blood'      # Boy vs. Girl 
#sampleType='Blood_Vessel'      # Boy vs. Girl 
#sampleType='Brain'      # Boy vs. Girl 
#sampleType='Colon'      # Boy vs. Girl 
#sampleType='Esophagus'      # Boy vs. Girl
#sampleType='Heart'      # Boy vs. Girl 
#sampleType='Liver'      # Boy vs. Girl 
#sampleType='Lung'      # Boy vs. Girl 
#sampleType='Pancreas'      # Boy vs. Girl 
#sampleType='Pituitary'      # Boy vs. Girl 
#sampleType='Small_Intestine'      # Boy vs. Girl 
#sampleType='Spleen'      # Boy vs. Girl 
#sampleType='Thyroid'      # Boy vs. Girl 
#sampleType='Breast'      # Boy vs. Girl 
#sampleType='Muscle'      # Boy vs. Girl 
#sampleType='Nerve'      # Boy vs. Girl 
#sampleType='Skin'      # Boy vs. Girl 
#sampleType='Stomach'      # Boy vs. Girl 

###############
# Directories #
###############
resultDir <-'~/results'
#MetresultDir <-paste0(resultDir,'/Methyl-Seq/',myProject)
RNAresultDir <-paste0(resultDir,'/RNA-Seq/',myProject)
deg.metaDir  <-paste0(RNAresultDir,'/Meta')
#dex.metaDir  <-paste0(RNAresultDir,'/Meta.DEXSeq')
#RDataDir <-paste0(RNAresultDir,'/RData'); if(!file.exists(RDataDir)){dir.create(RDataDir)}
clusterHome <- paste0(RNAresultDir,'/Cluster'); if(!file.exists(clusterHome)){dir.create(clusterHome)}
cluster.dir=paste0(clusterHome,'/',sampleType); if(!file.exists(cluster.dir)){dir.create(cluster.dir)}
if(myProject!="PDN"){
	#edgeRHome <- paste0(RNAresultDir,'/edgeR'); if(!file.exists(edgeRHome)){dir.create(edgeRHome)}
	#edgeR.dir=paste0(edgeRHome,'/',sampleType); if(!file.exists(edgeR.dir)){dir.create(edgeR.dir)}
	#DegHome <- paste0(RNAresultDir,'/DEG'); if(!file.exists(DegHome)){dir.create(DegHome)}
	#deg.dir=paste0(DegHome,'/',sampleType); if(!file.exists(deg.dir)){dir.create(deg.dir)}
	#deseqHome <- paste0(RNAresultDir,'/DESeq2'); if(!file.exists(deseqHome)){dir.create(deseqHome)}
	deseqHome <- paste0(RNAresultDir,'/', paste("DESeq2",packageVersion("DESeq2"),sep=".")); if(!file.exists(deseqHome)){dir.create(deseqHome)}
	deseq.dir=paste0(deseqHome,'/',sampleType); if(!file.exists(deseq.dir)){dir.create(deseq.dir)}
}else{
	deseq.dir=RDataDir
}

if(myCaller=="DEXSeq"){
	metaDir=dex.metaDir
	flattenedfile="/home/ssg29/data/Annotation/GRCh37.genes.dexseq.gff"
	dexseqHome <- paste0(RNAresultDir,'/DEXSeq'); if(!file.exists(dexseqHome)){dir.create(dexseqHome)}
	dexseq.dir=paste0(dexseqHome,'/',sampleType); if(!file.exists(dexseq.dir)){dir.create(dexseq.dir)}
}else{
	metaDir=deg.metaDir
}

#################
# Annotation DB #
#################
is.BM=FALSE
# miRNA from the Small-RNA-Seq
if(grepl("miRNA",myProject)){ 
	my.bed=ifelse(TR_PREFIX=='GRCh38',"~/data/Annotation/miRbase.21/hsa.gff3","~/data/Annotation/miRbase.20/hsa.gff3")
	my.filter="mirbase_id"
	my.id=my.filter #my.id="hgnc_symbol" # gene name
	#my.fields <- c("chromosome_name", my.filter, "ensembl_gene_id", "hgnc_symbol","description", "gene_biotype") # hsa-mir-106a (my.filter) has 6 hgnc_symbol (my.id)
	my.target<-rtracklayer::import.gff3(my.bed) # isa 'GRanges'
	my.target.list<-split(my.target, mcols(my.target)[["Name"]]) # isa 'GRangeList'
# piRNA from the Small-RNA-Seq
}else if(grepl("piRNA",myProject)){ 
	my.bed=paste0("~/data/Annotation/piRBase.10/piR_",TR_PREFIX,"_v1.0.bed.gz") #~/data/Annotation/piRBase.10/piR_GRCh37_v1.0.bed.gz
	my.filter="pirbase_id"
	my.id="pirbase_id"
	my.target<-rtracklayer::import.bed(my.bed) # isa 'GRanges'
	my.target.list<-split(my.target, mcols(my.target)[["name"]]) # isa 'GRangeList'
# Boy.Girl.exon.GRCh38 | RoadMap.exon.GRCh37 
}else if(grepl("exon",myProject)){ 
	# exon-level expression analysis for a targetted gene(s)
	# exon definition
	my.bed=ifelse(TR_PREFIX=='GRCh38',"/home/ssg29/data/genome/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.82.CSMD1.exon.union.sorted.gtf",
				  					"/home/ssg29/data/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.75.CSMD1.exon.union.sorted.gtf")
	my.filter="ensembl_exon_id"
	my.id="hgnc_symbol" # gene name
	my.fields <- c("chromosome_name", my.filter, "hgnc_symbol","description", "gene_biotype")
	# exon definition
	my.target<-rtracklayer::import(my.bed)
	my.target.list<-split(my.target, mcols(my.target)[["exon_id"]]) # isa 'GRangeList'
# re-constructed trascriptome based on our own RNA-Seq
}else if(grepl("RPT",myProject)){ 
	my.bed="~/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.gtf"
	my.filter=ifelse(grepl("TX",myProject),"transcript_id","gene_name") # gene_name for 'hgnc_symbol'
	my.id=my.filter
	my.target<-rtracklayer::import(my.bed) # isa 'GRanges'
	my.target.list<-split(my.target, mcols(my.target)[[my.filter]]) # isa 'GRangeList'
	#my.target.list<-split(my.target, mcols(my.target)[["nearest_ref"]]) # isa 'GRangeList'
# total RNA-Seq (quantified either featureCount or Salmon)
}else{ 
	my.filter="ensembl_gene_id"
	my.id="hgnc_symbol"
	my.fields <- c("chromosome_name", my.filter, "hgnc_symbol","description", "gene_biotype")
	my.target.list<-split(gr.exon, mcols(gr.exon)[["gene_id"]]) # gr.exon from config/Annotation.R
}

############
# cuff-off #
############
myFDR<-0.4
myPValue<-1.0e-02 #  .01
#mylogFC<-0.4 # 1.319508 FC
#mylogFC<-0.33 #1.257013 FC (only for DESeq2)
#mylogFC<-0.3 # 1.231144 FC
#mylogFC<-0.2 # 1.148698 FC
#mylogFC<-0.18 # 1.132884 FC
#mylogFC<-0.1 # 1.071773 FC
minRead<-10 # (only for DESeq2)

#############
## Samples ##
#############
selectSamples<-function(my.type, my.all){
	# FG samples
	if(grepl("^ALL",my.type)){
		my.samples<-my.all
	}else if(my.type=='LOW.PAPPA'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$PAPPA) & my.all$PAPPA==1,]$Pair,] # Low PAPPA SGA & matched AGA
	}else if(my.type=='NORMAL.PAPPA'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$PAPPA) & my.all$PAPPA==0,]$Pair,] # Normal PAPPA SGA & matched AGA
	}else if(my.type=='LOW.GVEL'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$GVel) & my.all$GVel==1,]$Pair,] # Low GVel SGA & matched AGA
	}else if(my.type=='NORMAL.GVEL'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$GVel) & my.all$GVel==0,]$Pair,] # Normal GVel SGA & matched AGA
	}else if(my.type=='HIGH.UAD'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$UAD) & my.all$UAD==1,]$Pair,] # High (Abnormal) UAD SGA & matched AGA
	}else if(my.type=='NORMAL.UAD'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$UAD) & my.all$UAD==0,]$Pair,] # Normal UAD SGA & matched AGA
	}else if(my.type=='HIGH.UBD'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$UBD) & my.all$UBD==1,]$Pair,] # High (Abnormal) UBD SGA & matched AGA
	}else if(my.type=='NORMAL.UBD'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$UBD) & my.all$UBD==0,]$Pair,] # Normal UBD SGA & matched AGA
	}else if(my.type=='PAPPA.SGA'){
		my.samples<-my.all[my.all$Condition==1 & !is.na(my.all$PAPPA),] # Low PAPPA SGA & unmatched normal SGA
	}else if(my.type=='HT'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & my.all$HT==1,]$Pair,] # Non-Hypertensive SGA & matched AGA
	}else if(my.type=='NHT'){
		my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & my.all$HT==0,]$Pair,] # Non-Hypertensive SGA & matched AGA
	}else if(my.type=='PET'){ 
		if(grepl("^SGA.AGA",myProject)){ # FG Samples only
			my.samples<-my.all[my.all$Pair %in% my.all[my.all$Condition==1 & !is.na(my.all$PET) & my.all$PET==1,]$Pair,] # PET SGA & matched AGA
		}else if(grepl("^Boy.Girl",myProject)){ # FG & JD
			my.samples<-my.all[my.all$Condition==1 & my.all$Source=="JD",]
		}else{
			stop(paste(myProject, "not supported for PET\n"))
		}
	}else if(my.type=='AGA'){ # healthy babies only (no breach)
		my.samples<-my.all[my.all$Condition==0 & my.all$Source!="JD-BR",]
	}else if(my.type=='SGA'){ # un-healthy (i.e small) babies only from FG
		my.samples<-my.all[my.all$Condition==1 & my.all$Source=="FG",]
	}else if(my.type=='CASE'){ # either SGA or PET
		my.samples<-my.all[my.all$Condition==1,]
	}else if(my.type=='BR'){ # Breech samples from JD 
		my.samples<-my.all[my.all$Source=="JD-BR",]
	}else if(my.type=='BR.Salmon'){ # Breech samples from JD 
		my.samples<-my.all[my.all$Source=="JD-BR",]
	}else if(my.type=='Boys'){ # Boys only (n=44)
		my.samples<-my.all[my.all$Pair %in% names(which(sapply(split(my.all[my.all$Sex=="M",], my.all[my.all$Sex=="M",]$Pair), nrow)==2)),] # Boys
	}else if(my.type=='Girls'){ # Girls only (n=56)
		my.samples<-my.all[my.all$Pair %in% names(which(sapply(split(my.all[my.all$Sex=="F",], my.all[my.all$Sex=="F",]$Pair), nrow)==2)),] # Girls
	}else if(my.type=='LOW.PAPPA.Boys'){
		my.samples<-my.all[my.all$Pair %in% names(which(sapply(split(my.all[my.all$Sex=="M",], my.all[my.all$Sex=="M",]$Pair), nrow)==2)),] # Boys
		my.samples<-my.samples[my.samples$Pair %in% my.samples[my.samples$Condition==1 & !is.na(my.samples$PAPPA) & my.samples$PAPPA==1,]$Pair,] # Low PAPPA
	}else if(my.type=='LOW.PAPPA.Girls'){
		my.samples<-my.all[my.all$Pair %in% names(which(sapply(split(my.all[my.all$Sex=="F",], my.all[my.all$Sex=="F",]$Pair), nrow)==2)),] # Girls 
		my.samples<-my.samples[my.samples$Pair %in% my.samples[my.samples$Condition==1 & !is.na(my.samples$PAPPA) & my.samples$PAPPA==1,]$Pair,] # Low PAPPA
	#Boy.Girl
	}else if(sampleType=='MJ.AGA'|sampleType=="FG.AGA"){ 
		my.samples<-my.all[my.all$Condition==0,] # AGA only
	}else if(sampleType=='MJ.SGA'|sampleType=="FG.SGA"){ 
		my.samples<-my.all[my.all$Condition==1,] # SGA only
	#IA (Retosiban) samples
	}else if(my.type=='Stretch'){
		my.samples<-my.all[my.all$Condition %in% c("L", "H"),]
		my.samples$Condition<-relevel(my.samples$Condition,"L")
	}else if(my.type=='R8'){
		my.samples<-my.all[my.all$Condition %in% c("H", "R8"),]
		my.samples$Condition<-relevel(my.samples$Condition,"H")
	}else if(my.type=='R6'){
		my.samples<-my.all[my.all$Condition %in% c("H", "R6"),]
		my.samples$Condition<-relevel(my.samples$Condition,"H")
	}else if(my.type=='Dosage'){
		my.samples<-my.all[my.all$Condition %in% c("H", "R8", "R6"),]
		my.samples$Condition<-relevel(my.samples$Condition,"H")
	#IA Trophoblast samples
	}else if(my.type=='M.DFMO'){
		my.samples<-my.all[my.all$Sex=="M" & my.all$Condition %in% c("Vehicle", "DFMO"),]
		my.samples$Condition<-relevel(my.samples$Condition,"Vehicle")
	}else if(my.type=='F.DFMO'){
		my.samples<-my.all[my.all$Sex=="F" & my.all$Condition %in% c("Vehicle", "DFMO"),]
		my.samples$Condition<-relevel(my.samples$Condition,"Vehicle")
	}else if(my.type=='M.siRNA'){
		my.samples<-my.all[my.all$Sex=="M" & my.all$Condition %in% c("Control", "siRNA"),]
		my.samples$Condition<-relevel(my.samples$Condition,"Control")
	}else if(my.type=='F.siRNA'){
		my.samples<-my.all[my.all$Sex=="F" & my.all$Condition %in% c("Control", "siRNA"),]
		my.samples$Condition<-relevel(my.samples$Condition,"Control")
	# Plasma.2017 samples
    }else if(my.type=="GA.12"|my.type=="GA.20"|my.type=="GA.28"|my.type=="GA.36"){
        my.samples<-my.all[paste0("GA.",my.all$GA)==my.type,]
	# Plasma.2016 samples
	}else if(my.type=='PLS.PT'){
		my.samples<-my.all[my.all$Sex=="M" & my.all$Condition==0,] # all healthy (condition=0) males (sex=M) and 12wk GA for placenta
		my.samples$Source<-relevel(my.samples$Source,"Placenta")
	}else if(my.type=='PLS12WK.PT'){
		my.samples<-my.all[my.all$Sex=="M" & my.all$Condition==0 & my.all$GA %in% c("40", "12"),] # all healthy (condition=0) males (sex=M) and 12wk GA for placenta
		my.samples$Source<-relevel(my.samples$Source,"Placenta")
	}else if(my.type=='PLS36WK.PT'){
		my.samples<-my.all[my.all$Sex=="M" & my.all$Condition==0 & my.all$GA %in% c("40", "36"),] # all healthy (condition=0) males (sex=M) and 12wk GA for placenta
		my.samples$Source<-relevel(my.samples$Source,"Placenta")
	}else if(my.type=='PLS36WK.12WK'){
		my.samples<-my.all[my.all$Sex=="M" & my.all$Condition==0 & my.all$Source=="Plasma",] # all heathy male plasma samples only 
		my.samples$Source<-relevel(my.samples$Source,"12")
	# RoadMap
	}else if(sampleType=="Schultz"|sampleType=="Consortium"|sampleType=="Fetal"){ 
		my.samples<-my.all[my.all$Source==sampleType,]
	# RoadMap fetal + Manchester Embryo
	}else if(sampleType=="FetalEmbryo"){ 
		my.samples<-my.all[my.all$Origin=="Fetal-Somatic"|my.all$Origin=="Embryo-Somatic",]
	}else{
		stop("No such my.type to subset")
	}
	return(my.samples)
}# end of function 'selectSamples'

##################
## SAMPLES INFO ##
##################
#SGA.AGA, SGA.AGA.miRNA, SGA.AGA.piRNA 
if(grepl("^SGA.AGA",myProject)){
	if(sampleType=='OXBS'){ # 
		samples <- read.csv(paste0(metaDir,"/meta.",sampleType,".csv"))
		samples$Condition<-as.factor(samples$Condition)
	}else{
		samples.all = read.csv(paste0(metaDir,"/meta.ALL.csv"), stringsAsFactors=FALSE)
		samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # column 3-End as factor
		rownames(samples.all)<-samples.all$SampleName

		if(file.exists(paste0(metaDir,"/meta.Pheno.csv"))){
			samples.phen= read.csv(paste0(metaDir,"/meta.Pheno.csv"))
			samples.phen[c("MD","Smoke","FS","Smoker","PT")]=lapply(samples.phen[c("MD","Smoke","FS","Smoker","PT")], as.factor) # apply factor for the categorical values
			rownames(samples.phen)<-samples.phen$SampleName
		}
		#unaffected = read.csv(paste0(metaDir,"/unaffected.samples.txt"), stringsAsFactors=FALSE) # for mRNA-Seq only
		# select samples *unique* to the sub-type
		if(grepl(".UNIQ$", sampleType)){
			# 5 SGA sub-types
			sampleSubTypes=c("LOW.PAPPA", "LOW.GVEL", "HIGH.UAD", "HIGH.UBD", "HT") # 5 SGA sub-types
			subtype.list<-lapply(sampleSubTypes, selectSamples, samples.all) # list of 5 SGA sub-type samples
			names(subtype.list)<-sampleSubTypes # name each slice
			subtype.list<-lapply(subtype.list, function(i) i<-i[i$Condition==1,c("SampleName")]) # get SGA 'SampleName' only

			parentType<-unlist(strsplit(sampleType, ".UNIQ")) # e.g. LOW.PAPPA.UNIQ => LOW.PAPPA
			uniqSampleName<-Reduce(setdiff, subtype.list[c(parentType, setdiff(sampleSubTypes, parentType))]) # SampleName unique to parentType 
			samples<-selectSamples(parentType, samples.all) # samples for parentType
			samples<-samples[samples$Pair %in% samples[samples$SampleName %in% uniqSampleName,]$Pair,]
		}else if(grepl(".EXC$", sampleType)){
			# 5 SGA sub-types
			sampleSubTypes=c("LOW.PAPPA", "LOW.GVEL", "HIGH.UAD", "HIGH.UBD", "HT") # 5 SGA sub-types
			subtype.list<-lapply(sampleSubTypes, selectSamples, samples.all) # list of 5 SGA sub-type samples
			names(subtype.list)<-sampleSubTypes # name each slice
			subtype.list<-lapply(subtype.list, function(i) i<-i[i$Condition==1,c("SampleName")]) # get SGA 'SampleName' only

			parentType<-unlist(strsplit(sampleType, ".EXC")) # e.g. LOW.PAPPA.EXC=> LOW.PAPPA
			uniqSampleName<-Reduce(setdiff, subtype.list[c(parentType, setdiff(sampleSubTypes, parentType))]) # SampleName unique to parentType 

			if(sampleType=="LOW.PAPPA.EXC"){
				# LOW.PAPPA.UNIQ + intersect(LOW.PAPPA, HT)
				# 7 = 3 + 4 
				uniqSampleName=c(uniqSampleName, intersect(subtype.list$LOW.PAPPA, subtype.list$HT))
			}else if(sampleType=="LOW.GVEL.EXC"){
				# 7 = 7
				uniqSampleName=uniqSampleName
			}else if(sampleType=="HIGH.UAD.EXC"){
				# HIGH.UAD.UNIQ + intersect(HIGH.UAD + LOW.GVEL) + intersect(HIGH.UAD + HIGH.UBD)
				# 6 = 2 + 3 + 1
				uniqSampleName=c(uniqSampleName, intersect(subtype.list$HIGH.UAD, subtype.list$LOW.GVEL), intersect(subtype.list$HIGH.UAD, subtype.list$HIGH.UBD))
			}else if(sampleType=="HIGH.UBD.EXC"){
				# HIGH.UBD.UNIQ + intersect(HIGH.UBD, LOW.GVEL)
				# 8 = 5 + 3
				uniqSampleName=c(uniqSampleName, intersect(subtype.list$HIGH.UBD, subtype.list$LOW.GVEL))
			}else if(sampleType=="HT.EXC"){
				# HT.UNIQ + intersect(HT, HIGH.UAD) - LOW.PAPPA
				# 8 = 3 + 6 - 1
				uniqSampleName=setdiff( c(uniqSampleName, intersect(subtype.list$HT, subtype.list$HIGH.UAD)), subtype.list$LOW.PAPPA)
			}else{
				stop("No such exclusive sampleType to subset")
			}

			samples<-selectSamples(parentType, samples.all) # samples for parentType
			samples<-samples[samples$Pair %in% samples[samples$SampleName %in% uniqSampleName,]$Pair,]
		}else{
            ## SAMPLE EXCLUSION AS PER JUSTYNA AND FRANCESCA (ORIGINALLY FROM THE LANCET PAPER) 30/APR/2018
            ## SampleName: 88 decidual contaminated
            my.blacklist<-samples.all$SampleName=="88" | samples.all$PET==1 | samples.all$HT==1 | samples.all$GH==1
            my.bad.pair<-as.character(samples.all[my.blacklist,]$Pair)
            # sex-mismatch
            my.bad.pair<-c(my.bad.pair,"P13","P5")
            samples<-samples.all[!samples.all$Pair %in% my.bad.pair,] # n=80 (40 SGA & 40 matched controls)
            # sex-mismatch
            #dt.samples<-data.table(samples.all[!samples.all$Pair %in% samples.all[my.blacklist,]$Pair,])
            #dt.samples[,.N,"Pair,Sex"]
            #dcast.data.table(dt.samples[,.N,"Pair,Sex"], Pair~Sex)[!is.na(F*M),Pair] 
		}
    }
	##############
	## Encoding ##
	##############
	if(FALSE){
		samples$Code<-'A'
		samples$Code[samples$Condition==1]='S'
		# PAPPA
		samples$Code[is.na(samples$PAPPA)]=paste0(samples$Code[is.na(samples$PAPPA)],'*')
		samples$Code[!is.na(samples$PAPPA) & samples$PAPPA==1]=paste0(samples$Code[!is.na(samples$PAPPA) & samples$PAPPA==1],'p')
		samples$Code[!is.na(samples$PAPPA) & samples$PAPPA==0]=paste0(samples$Code[!is.na(samples$PAPPA) & samples$PAPPA==0],'P')
		# HT 
		samples$Code[samples$HT==1]=paste0(samples$Code[samples$HT==1],'h')
		samples$Code[samples$HT==0]=paste0(samples$Code[samples$HT==0],'H')
		# GVel
		samples$Code[is.na(samples$GVel)]=paste0(samples$Code[is.na(samples$GVel)],'*')
		samples$Code[!is.na(samples$GVel) & samples$GVel==1]=paste0(samples$Code[!is.na(samples$GVel) & samples$GVel==1],'g')
		samples$Code[!is.na(samples$GVel) & samples$GVel==0]=paste0(samples$Code[!is.na(samples$GVel) & samples$GVel==0],'G')
		# ubd
		samples$Code[is.na(samples$UBD)]=paste0(samples$Code[is.na(samples$UBD)],'*')
		samples$Code[!is.na(samples$UBD) & samples$UBD==1]=paste0(samples$Code[!is.na(samples$UBD) & samples$UBD==1],'b')
		samples$Code[!is.na(samples$UBD) & samples$UBD==0]=paste0(samples$Code[!is.na(samples$UBD) & samples$UBD==0],'B')
		# UAD
		samples$Code[is.na(samples$UAD)]=paste0(samples$Code[is.na(samples$UAD)],'*')
		samples$Code[!is.na(samples$UAD) & samples$UAD==1]=paste0(samples$Code[!is.na(samples$UAD) & samples$UAD==1],'u')
		samples$Code[!is.na(samples$UAD) & samples$UAD==0]=paste0(samples$Code[!is.na(samples$UAD) & samples$UAD==0],'U')
	}
#Boy.Girl, Boy.Girl.miRNA, Boy.Girl.piRNA
}else if(grepl("^Boy.Girl",myProject)){
	# black list to exclude (contaminated)
	#outLiers=list(`IGFBP1`=c("06C","80C","84C","88"),`SNORD3C`=c("02C","87C"), `Breech`=c("93C","94C"))
	#outLiers=list(`IGFBP1`=c("80C","84C","88"), `Breech`=c("93C","94C")) 
    outLiers=list(`Failed`=c("91"),`IGFBP1`=c("80C","84C","88")) # as of 17/Apr/2018
	if(grepl("iRNA",myProject)){outLiers[["Failed"]]=c("20P","20C")}

	#samples.all = read.csv(paste0(metaDir,"/meta.ALL.Salmon.csv"), stringsAsFactors=FALSE)
	samples.all = read.csv(paste0(metaDir,"/meta.ALL.csv"), stringsAsFactors=FALSE) # n=321 for ALL samples (FG=114, JD=188, JD-BR=19)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # 3rd column to the end: as factor 
	rownames(samples.all)<-samples.all$SampleName
	samples<-selectSamples(sampleType, samples.all)
	# add 'SampleId' column if not there (for 'collapseReplciates')
	if(!any(colnames(samples.all) %in% c("SampleId"))){
		samples$SampleId<-samples$SampleName
	}
	cat("Removing samples from the black list (outliers)...\n")
	samples<-samples[!samples$SampleName %in% unname(unlist(outLiers)),]
}else if(grepl("^PET",myProject)){
    # 12-pairs excluded either batch-effect or decidual contamination (30/APR/2018)
	outLiers=c("PET08","PET16","PET20","PET60","PET64","PET75","PET76","PET77","PET78","PET79","PET80","PET84")

	samples.all = read.csv(paste0(metaDir,"/meta.",sampleType,".csv"),stringsAsFactors=FALSE)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # 3rd column to the end: as factor 
	rownames(samples.all)<-samples.all$SampleName
	samples<-selectSamples(sampleType, samples.all)

	cat("Removing samples from the black list (outliers)...\n")
	samples<-samples[!samples$Pair%in% outLiers,]
}else if(myProject=="PDN"){
	samples.all <- read.csv(paste0(metaDir,"/meta.",sampleType,".csv"))
	samples<-samples.all
}else if(grepl("^Retosiban",myProject)){
	samples.all <- read.csv(paste0(metaDir,"/meta.",myProject,".csv"), stringsAsFactors=FALSE)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # column 3-End as factor 
	rownames(samples.all)<-samples.all$SampleName
	samples<-selectSamples(sampleType, samples.all)
}else if(myProject=="Trophoblast.IA"){
	samples.all = read.csv(paste0(metaDir,"/meta.ALL.csv"), stringsAsFactors=FALSE)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # column 3-End as factor 
	rownames(samples.all)<-samples.all$SampleName
	samples<-selectSamples(sampleType, samples.all)
}else if(grepl("^Plasma",myProject)){
	#samples.all = read.csv(paste0(metaDir,"/meta.",sampleType,".csv"),stringsAsFactors=FALSE)
	samples.all = read.csv(paste0(metaDir,"/meta.ALL.Salmon.csv"), stringsAsFactors=FALSE)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # column 3-End as factor 
	rownames(samples.all)<-samples.all$SampleName
	#samples<-samples.all
    samples<-selectSamples(sampleType, samples.all)
}else if(grepl("^RoadMap",myProject)){
	samples.all <- read.csv(paste0(metaDir,"/meta.ALL.csv"), stringsAsFactors=FALSE)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # column 3-End as factor 
	rownames(samples.all)<-samples.all$SampleName
	samples<-selectSamples(sampleType, samples.all)
}else if(grepl("^GTEx",myProject)){
	samples.all <- read.csv(paste0(metaDir,"/meta.ALL.csv"), stringsAsFactors=FALSE)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # column 3-End as factor 
	rownames(samples.all)<-samples.all$SampleName
	samples<-samples.all[samples.all$Tissue==sampleType,]
}else{
	samples.all <- read.csv(paste0(metaDir,"/meta.ALL.csv"), stringsAsFactors=FALSE)
	samples.all[c(3:ncol(samples.all))]=lapply(samples.all[c(3:ncol(samples.all))], as.factor) # column 3-End as factor 
	rownames(samples.all)<-samples.all$SampleName
	samples<-samples.all
	#stop("No such project\n")
}# end of Sample Info

if("Sex" %in% colnames(samples)){
	cat("releveling sex (F-M)...\n")
	samples$Sex<-relevel(samples$Sex,"M") # F-M (case: female, control: male)
}

# adjust levels
cat("Adjusting levels...\n")
for(i in colnames(samples[3:ncol(samples)]))(samples[i]<-droplevels(samples[i]))

############
## Design ##
############
#SGA.AGA, SGA.AGA.piRNA, SGA.AGA.miRNA
if(grepl("^SGA.AGA",myProject) | grepl("^PET",myProject) | grepl("^Retosiban",myProject) | grepl("^Cholestasis",myProject) | grepl("^Trophoblast.IA",myProject)){
	my.contrast="Condition"
	if(sampleType=='PAPPA.SGA'){
		design <- model.matrix(~PAPPA, samples) # isa 'matrix'
		deseq.design <- formula(~ PAPPA) 		# simple design.
		dexseq.design<- ~ sample + exon + PAPPA:exon 
	}else if(sampleType=='PET.SGA'){
		design <- model.matrix(~HT, samples) 	# isa 'matrix'
		deseq.design <- formula(~ HT) 			# simple design.
		dexseq.design<- ~ sample + exon + HT:exon 
	}else if(sampleType=='AGA' | sampleType=='SGA'){ 	# (un)healty babies only
		design <- model.matrix(~SampleName, samples) 		# isa 'matrix'
		deseq.design <- formula(~ SampleName) 			# individual difference (isa 'formula')
		dexseq.design<- ~ sample + exon + SampleName:exon 
	}else{
		design <- model.matrix(~ Pair + Condition, samples) # isa 'matrix' with paired design
		deseq.design <- formula(~ Pair + Condition) 		# paired design. isa 'formula'
		dexseq.design<- ~ sample + exon + Pair:exon + Condition:exon 
		formulaReducedModel = ~ sample + exon + Pair:exon
	}
#Boy.Girl, Boy.Girl.miRNA, Boy.Girl.piRNA
}else if(grepl("^Boy.Girl",myProject) | grepl("^GTEx",myProject)){
	my.contrast="Sex"
	if(sampleType=='BR'){
	# paired-design for breech samples (BR)
		design <- model.matrix(~ Pair + Sex, samples) # isa 'matrix' with paired design
		deseq.design <- formula(~ Pair + Sex) 		# paired design. isa 'formula'
		dexseq.design<- ~ sample + exon + Pair:exon + Sex:exon 
    }else if(sampleType=='BR.Salmon'){
	# paired-design for breech samples (BR.Salmon)
		design <- model.matrix(~ Pair + Sex, samples) # isa 'matrix' with paired design
		deseq.design <- formula(~ Pair + Sex) 		# paired design. isa 'formula'
		dexseq.design<- ~ sample + exon + Pair:exon + Sex:exon 
	}else{
	# simple design
		design <- model.matrix(~Sex, samples) 
		deseq.design <- formula(~ Sex) 
		dexseq.design<- ~ sample + exon + Sex:exon 
	}
}else if(myProject=="PDN"){
	my.contrast="Condition"
	design <- model.matrix(~SampleName, samples) 	# isa 'matrix'
	deseq.design <- formula(~ SampleName) 			# individual difference (isa 'matrix')
	dexseq.design<- ~ sample + exon + SampleName:exon 
}else if(grepl("^Plasma",myProject)){
	# Plasma.2016: PLS.PT, PLS12WK.PT or PLS36WK.PT
	if(grepl(".PT$",sampleType)){ 
		my.contrast="Source"
		design <- model.matrix(~ Source, samples) # isa 'matrix'
		deseq.design <- formula(~ Source) 		  # isa 'formula'
		dexseq.design<- ~ sample + exon + Source:exon 
	# Plasma.2016: PLS36WK.12WK
	}else if(grepl(".WK$",sampleType)){ 
		my.contrast="GA"
		design <- model.matrix(~ GA, samples) # isa 'matrix'
		deseq.design <- formula(~ GA) 		  # isa 'formula'
		dexseq.design<- ~ sample + exon + GA:exon 
	# ALL, ALL.Salmon, GA.12 
    }else{
		my.contrast="Sex"
		# simple design
		design <- model.matrix(~Sex, samples) 
		deseq.design <- formula(~ Sex) 
		dexseq.design<- ~ sample + exon + Sex:exon 
    }
}else if(grepl("^RoadMap",myProject)){
	my.contrast="Tissue"
	design <- model.matrix(~ Tissue, samples) # isa 'matrix'
	deseq.design <- formula(~ Tissue) 		  # isa 'formula'
}else{
	stop("No such project\n")
}

assignFDRcol <-function(my.top.deg){
	cat("assigning FDR colour\n")
	my.top.deg$col<-vector(length=nrow(my.top.deg))
	for(i in 1:nrow(my.top.deg)){
		if(myCaller=="DESeq2"){
			myFDR=my.top.deg[i,]$padj
		}else if(myCaller=="edgeR"){
			myFDR=my.top.deg[i,]$FDR
		}else{
			stop("myCaller not defined. Now Stopped!\n")
		}  

		if(myFDR <= 0.1){my.top.deg[i,]$col=myCol[1]}
		if(myFDR > 0.1 & myFDR <= 0.2 ){my.top.deg[i,]$col=myCol[2]}
		if(myFDR > 0.2 & myFDR <= 0.3 ){my.top.deg[i,]$col=myCol[3]}
		if(myFDR > 0.3 & myFDR <= 0.4 ){my.top.deg[i,]$col=myCol[4]}
		if(myFDR > 0.4 & myFDR <= 0.6 ){my.top.deg[i,]$col=myCol[5]}
		if(myFDR > 0.6 ){my.top.deg[i,]$col=myCol[6]}
	}
	return(my.top.deg)
}

########################################
### Plot Expression Level of Top DEGs ##
########################################
#my.entry: ensg_id
plotExp<-function(my.entry,my.metric="FPM",my.contrast="Condition",my.box=TRUE,my.save=FALSE){
	my.gene <- my.entry 

	foo <- data.frame(`Count`=counts(dds)[my.entry,], `FPKM`=ddsFpkm[my.entry,], `FPM`=ddsFpm[my.entry,],colData(dds))
	foo$SampleName <- rownames(foo)
	foo$Condition <- as.factor(ifelse(foo$Condition==1, "Case", "Control"))

	#from DESeq2
	my.pval<-round(my.res[my.entry,"pvalue"],2)
	my.padj<-round(my.res[my.entry,"padj"],2)
	my.fc<-round(2^abs(my.res[my.entry,"log2FoldChange"]),2)
	# p-value from wilcox test
	if(grepl("^Boy.Girl",myProject)){
		my.new.pval<-wilcox.test(foo[foo$Sex=="M","FPKM"], foo[foo$Sex=="F","FPKM"])$p.value
	}else{
		my.new.pval<-my.pval
	}
	#my.title=paste0("Expression level of ",my.gene," from ",sampleType, " (p-value: ", my.new.pval, ", p-value(DESeq2):",my.pval ,", q-value(DESeq2):",my.padj,", FC:",my.fc,")")
	#p1: per-group boxplot
	my.labels=ifelse(grepl("^Boy.Girl",myProject),"Gender",my.contrast)
	p1<-ggplot(foo, aes_string(x=my.contrast, y=my.metric)) +
		geom_boxplot() + 
		#geom_jitter(width=0.2) +
		geom_dotplot(aes_string(fill=my.contrast), binaxis='y', stackdir='center', dotsize=rel(.8)) +
		scale_fill_manual(values=my.col[[my.contrast]]) +
		scale_x_discrete(labels=names(my.col[[my.labels]])) +
		theme_Publication() +
		theme(legend.position="none")

	#p2: per-sample barchart
    if("Pair" %in% colnames(colData(dds))){
        p2<-ggplot(foo, aes_string(x="Pair", y=my.metric)) +
            geom_bar(stat="identity",aes_string(fill=my.contrast), position="dodge") + 
            ggtitle(paste(my.gene,"(padj=",my.padj,")")) + 
            scale_fill_manual(values=my.col[[my.contrast]]) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }else{
        p2<-ggplot(foo, aes_string(x="SampleName", y=my.metric, fill=my.contrast)) +
            geom_bar(stat="identity") + 
            ggtitle(paste(my.gene,"(padj=",my.padj,")")) + 
            scale_fill_manual(values=my.col[[my.contrast]]) +
            theme_Publication()
    }

	if(my.save){
		if(my.box){
			my.filename <- file.path(deseq.dir, paste(my.gene, "boxplot.exp.level",my.metric,sep="."))
            p.this<-p1
		}else{
			my.filename <- file.path(deseq.dir, paste(my.gene, "exp.level",my.metric,sep="."))
            p.this<-p2
		}
        #if(capabilities()[["tiff"]]){
        #    tiff(filename=paste0(my.filename,".tiff"),width=14,height=9,units="in",res=300, compression = 'lzw')
        #}else{
        #    jpeg(filename=paste0(my.filename,".jpeg"),width=14,height=9,units="in",res=300)
        #}
        print(p.this)
		#dev.off()
	}else{
		if(my.box){print(p1)}else{print(p2)}
	}
}# endof plotExp

#called by 'bin/R/MJ.RNA/get.expression.profile.R'
plotExpMulti<-function(my.entry,my.metric="Count"){
	if(my.metric=="Count"){
		foo<-reshape2::melt(counts(dds)[my.entry,])
	}else if(my.metric=="FPM"){
		foo<-reshape2::melt(ddsFpm[my.entry,])
	}else{
		foo<-reshape2::melt(ddsFpkm[my.entry,])
	}
	foo$Sample<-rownames(foo); dummy.sample<-"Sample"
	colnames(foo)<-c(my.metric,dummy.sample)
	foo<-cbind(foo, colData(dds))

	# by fetal sex
	p1 <-ggplot(foo,aes_string(x=dummy.sample,y=my.metric)) + 
				geom_bar(aes(fill=Sex),stat='identity') + 
				theme_Publication() +
				ggtitle(my.entry)

	# by condition
	p2 <-ggplot(foo,aes_string(x=dummy.sample,y=my.metric)) + 
				geom_bar(aes(fill=Condition),stat='identity') + 
				theme_Publication() +
				ggtitle(my.entry)

	# by collection time? 
	#p3 <-ggplot(foo,aes_string(x=dummy.sample,y=my.metric)) + 
	#			geom_bar(aes(fill=Ptype),stat='identity') + 
	#			ggtitle(my.entry)

	# both fetal sex and condition
	p4 <-ggplot(foo,aes_string(x=dummy.sample,y=my.metric)) + 
				geom_bar(stat='identity') + 
				facet_grid(Sex ~ Condition) +
				ggtitle(my.entry) +
				theme(strip.text = element_text(face="bold"),
					axis.text.x = element_text(angle=45),
					#strip.background = element_rect(colour="#B2B26B", fill="#FFFF99")
					strip.background = element_rect(colour="#664C00", fill="#CC9900")
				)
	multiplot(p1, p2, p4, cols=3) # ~/lib/multiplot.R
}#end of plotFPM

#################
# Volcanio Plot #
#################
getTopDeseq<-function(my.top.deg, my.name){
    if(my.name=="all.gene"){
        p <- ggplot(my.top.deg, aes(log2FoldChange,-log10(padj))) + 
            geom_point(data=my.top.deg[padj>=0.05],alpha=0.5,size=1) + 
            geom_point(data=my.top.deg[padj<0.05 & log2FoldChange>=0],alpha=0.5,size=1,col="red") + 
            geom_point(data=my.top.deg[padj<0.05 & log2FoldChange<0],alpha=0.5,size=1,col="blue") + 
            geom_hline(aes(yintercept=-log10(0.05)),col="grey", linetype="dashed") +
            ggtitle(paste("All ", nrow(my.top.deg)," genes")) +
            theme_Publication()
	}else if(grepl("^top",my.name)){
        p <- ggplot(my.top.deg, aes(log2FoldChange,-log10(padj))) + 
            geom_point(alpha=0.8,size=1.5) + 
            geom_hline(yintercept=-log10(0.01),col="grey", linetype="dashed") +
            geom_vline(xintercept=c(-log2(1.5),log2(1.5)),col="grey", linetype="dashed") +
            geom_text(data=my.top.deg[log2FoldChange>=0],aes_string(label=my.id),col="red",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
            geom_text(data=my.top.deg[log2FoldChange<0],aes_string(label=my.id),col="blue",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
            ggtitle(paste("Top ", nrow(my.top.deg)," DEGs")) +
            theme_Publication()
    }else{
        p <- ggplot(my.top.deg, aes(log2FoldChange,-log10(padj))) + 
            geom_point(alpha=0.8,size=1.5) + 
            geom_hline(yintercept=-log10(my.pval),col="grey", linetype="dashed") +
            geom_vline(xintercept=c(-log2(1.5),log2(1.5)),col="grey", linetype="dashed") +
            geom_text(data=my.top.deg[log2FoldChange>=0],aes_string(label=my.id),col="red",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
            geom_text(data=my.top.deg[log2FoldChange<0],aes_string(label=my.id),col="blue",hjust=0,nudge_x=0.05, nudge_y=0.02,size=4.5) + 
            ggtitle(paste(nrow(my.top.deg),"DEGs of padj<",my.pval)) +
            theme_Publication()
    }
    print(p)
	
	################################ 
	# Write top deg into .csv file #
	################################ 
	write.csv(my.top.deg, file=gzfile(paste0(deseq.dir,"/filtered_toptags_deseq.",my.name,".",sampleType,".csv.gz")), row.names=F, quote=F) # isa 'Data.Frame'
}#end of getTopDeseq 

getTopDeseq_OLD_NO_USE <-function(my.top.deg, my.name){
	my.top.deg<-as.data.frame(my.top.deg)
	my.top.deg[[my.filter]]<-rownames(my.top.deg)

	minDeseqlogFC=round(min(abs(my.top.deg$log2FoldChange)),3)
	
    # miRNA, piRNA, exon, RPT
    if(exists("my.target")){
		#sapply(my.target.list[names(my.target.list) %in% rownames(my.top.deg),], function(i) as.character(seqnames(i)))
		bar<-unique(with(as.data.frame(my.target), data.frame(seqnames,Name))); colnames(bar)<-c("chromosome_name",my.filter)
		deseq.top.deg.anno=merge(my.top.deg, bar, all.x=TRUE)
    # ensembl_gene_id based analysis
    }else{
        deseq.top.deg.anno<-merge(my.top.deg, dt.ensg[,.(ensembl_gene_id,chromosome_name,hgnc_symbol,gene_biotype)], all.x=TRUE)
    }
	deseq.top.deg.anno<-assignFDRcol(deseq.top.deg.anno) # assignFDRcol defined within the config

	# Append FPM and FPKM
	if(nrow(deseq.top.deg.anno)>1){
		foo<-data.frame(
			meanFpm=rowMeans(ddsFpm[my.top.deg[my.top.deg[[my.filter]] %in% rownames(ddsFpm),my.filter],]),
			meanFpkm=rowMeans(ddsFpkm[my.top.deg[my.top.deg[[my.filter]] %in% rownames(ddsFpkm),my.filter],])
		)
		foo[[my.filter]]<-my.top.deg[my.top.deg[[my.filter]] %in% rownames(ddsFpm),my.filter]
		deseq.top.deg.anno<-merge(deseq.top.deg.anno, foo, all.x=TRUE)
	}else{
		deseq.top.deg.anno[["meanFpm"]]=mean(ddsFpm[deseq.top.deg.anno[deseq.top.deg.anno[[my.filter]] %in% rownames(ddsFpm),my.filter],])
		deseq.top.deg.anno[["meanFpkm"]]=mean(ddsFpkm[deseq.top.deg.anno[deseq.top.deg.anno[[my.filter]] %in% rownames(ddsFpkm),my.filter],])
	}

	############
	## ggplot ##
	############
	dummy<-unique(deseq.top.deg.anno[,c(my.id,"baseMean","log2FoldChange","padj","col")])
	my.main=paste0(nrow(my.top.deg)," top DEGs for ",sampleType," (logFC>=",minDeseqlogFC,") from ", my.name)
	p<-ggplot(dummy, aes(log10(baseMean),log2FoldChange)) + 
        geom_point() + 
        geom_hline(yintercept=c(-minDeseqlogFC, minDeseqlogFC), col="blue", linetype=2) + 
        ggtitle(my.main) + 
        geom_text(aes_string(label=my.id, colour="col")) + 
        scale_colour_manual(name='Adjusted P-value', values = unname(myCol[myCol %in% dummy$col]), labels=names(myCol[myCol %in% dummy$col]))
	print(q+theme_Publication())

	################################ 
	# Write top deg into .csv file #
	################################ 
	write.csv(deseq.top.deg.anno[order(deseq.top.deg.anno$pvalue),], file=gzfile(paste0(deseq.dir,"/filtered_toptags_deseq.",my.name,".",sampleType,".csv.gz")), row.names=F, quote=F) # isa 'Data.Frame'
}#end of getTopDeseq 

getTopEdgeR <-function(my.top.deg,my.name){
	minlogFC=round(min(abs(my.top.deg$logFC)),3)

	################
	## Annotation ##
	################
	top.deg.ens=getBM(attributes = my.fields, filters = my.filter, values = rownames(my.top.deg), mart = myMart)
	my.top.deg[,my.filter]<-rownames(my.top.deg)

	# join two data.frame by a common column
	top.deg.anno=as.data.frame(merge(my.top.deg, top.deg.ens, by=my.filter, all.x=TRUE))
	# FDR distribution
	#hist(top.deg.anno$FDR, breaks=10,xlab='FDR')

	####################
	## Plot Gene Name ##
	####################
	# assign color
	# myCol defined within the config
	top.deg.anno<-assignFDRcol(top.deg.anno) # assignFDRcol defined within the config
	minlogFC=round(min(abs(my.top.deg$logFC)),3)
	plot(y=top.deg.anno$logFC, x=top.deg.anno$logCPM, main=paste0('edgeR (',sampleType,'): Top ', length(rownames(top.deg.anno)),' DEGs from ',my.name), ylab='logFC', xlab='logCPM', cex=0.2)
	text(y=top.deg.anno$logFC, x=top.deg.anno$logCPM, labels=top.deg.anno$hgnc_symbol, col=top.deg.anno$col, cex=0.75)
	abline(h=c(-minlogFC, minlogFC), col="blue")
	legend("topright", legend=names(myCol), fill=myCol)

	write.csv(top.deg.anno[order(top.deg.anno$logFC,decreasing=TRUE),], file=gzfile(paste0(edgeR.dir,'/filtered_toptags_edgeR.',my.name,".",sampleType,'.csv.gz')))
}

# merge edgeR & DESeq2
mergeDEG<-function(my.type, my.name){
	cat(paste0("merging ", my.type, my.name,"...\n"))
	#~/results/RNA-Seq/SGA.AGA/edgeR/LOW.PAPPA/filtered_toptags_edgeR.FDR.4.LOW.PAPPA.csv
	#~/results/RNA-Seq/SGA.AGA/DESeq2/LOW.PAPPA/filtered_toptags_deseq.padj.4.LOW.PAPPA.csv
	if(my.type=="FDR"){
		edgeR.result1 = read.csv(paste0(edgeR.dir,"/filtered_toptags_edgeR.FDR.",my.name,".",sampleType,".csv"), stringsAsFactors=FALSE)
		deseq.result1 = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.padj.",my.name,".",sampleType,".csv"), stringsAsFactors=FALSE)
	}else if(my.type=="top"){
		edgeR.result1 = read.csv(paste0(edgeR.dir,"/filtered_toptags_edgeR.top",my.name,".",sampleType,".csv"), stringsAsFactors=FALSE)
		deseq.result1 = read.csv(paste0(deseq.dir,"/filtered_toptags_deseq.top",my.name,".",sampleType,".csv"), stringsAsFactors=FALSE)
	}else if(my.type=="all"){
		edgeR.result1 = read.csv(paste0(edgeR.dir,"/toptags_all_edgeR.",sampleType,".csv"), stringsAsFactors=FALSE)
		deseq.result1 = read.csv(paste0(deseq.dir,"/toptags_all_deseq.",sampleType,".csv"), stringsAsFactors=FALSE)

		names(edgeR.result1)[names(edgeR.result1) == 'X'] <- my.filter
		names(deseq.result1)[names(deseq.result1) == 'X'] <- my.filter
	}else{
		stop("No such option\n")
	}

	union.top.deg1=merge(edgeR.result1, deseq.result1, by=my.filter, all=T) # isa 'data.frame'   # union between two
	if(my.type=="all"){
		common.top.deg1=union.top.deg1[!is.na(union.top.deg1$logFC) & !is.na(union.top.deg1$log2FoldChange), ] # common DEG (intersection)
	}else{
		common.top.deg1=union.top.deg1[!is.na(union.top.deg1$X.x) & !is.na(union.top.deg1$X.y), ] # common DEG (intersection)
	}
	write.csv(union.top.deg1, file=gzfile(paste0(deg.dir,"/unified.top.deg.",my.type,".",my.name,".",sampleType,".csv.gz")))
	write.csv(common.top.deg1, file=gzfile(paste0(deg.dir,"/common.top.deg.",my.type,".",my.name,".",sampleType,".csv.gz")))
}

initDDS<-function(collapsed=FALSE){
	cat("DESeqDataSetFromHTSeqCount...\n")
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = samples, design= deseq.design ) # isa 'DESeqDataSet'
	if(collapsed){
		cat("collapsing ddsHTSeq by SampleName...\n")
		ddsCollapsed <- collapseReplicates(ddsHTSeq, groupby=ddsHTSeq$SampleName, run=ddsHTSeq$LibName) # isa 'DESeqDataSet'
																										# ddsHTSeq$SampleName will be rownames of colData(dds)
																										# ddsHTSeq$LibName will be join by ','
		dds <- DESeq(ddsCollapsed) # isa 'DESeqDataSet'
	}else{
		cat("init dds...\n")
		dds <- DESeq(ddsHTSeq) # isa 'DESeqDataSet'
	}
	cat("making rld...\n")
	rld <- rlog(dds) # isa 'SummarizedExperiment'

	####################################
	# set exon size per gene for FPKM  #
	####################################
	keep <- rownames(dds) %in% names(my.target.list) # remove ensg of no exons info
	newCounts <- counts(dds)[keep,]
	ddsMatrix <- DESeqDataSetFromMatrix(newCounts, colData(dds)[1:ncol(colData(dds))-1], deseq.design) # isa 'DESeqDataSet'
	ddsExonicGene <- DESeq(ddsMatrix) # isa 'DESeqDataSet'
	# 1. by rowData-way 
	rowData(ddsExonicGene) <- my.target.list 
	# 2. by mcols[["basepairs"]] (method 1 is much faster) (https://www.biostars.org/p/83901/)
	#exonic.gene.sizes <- lapply(my.target.list,function(x){sum(width(reduce(x)))}) # isa 'list'
	#mcols(ddsExonicGene)[["basepairs"]]<-unlist(exonic.gene.sizes[rownames(ddsExonicGene)]) # isa column 'vector'
	#mcols(ddsExonicGene)$basepairs<-exonic.gene.sizes[rownames(dds),] # isa column 'vector'

	cat("saving DESeq2 RData\n")
	if(collapsed){
		save(ddsHTSeq,dds,rld,ddsExonicGene, file=deseq.RData) #  save 'dds'
	}else{
		save(dds,rld,ddsExonicGene, file=deseq.RData) #  save 'dds'
	}
}


# rownames(my.table) should be 
# ensg for total RNA  
# mirbase_id for miRNA
# pirbase_id for piRNA
getTopExp<-function(my.table, iqr=FALSE, upto=60){
	if(iqr){
		select <- head( order( apply(my.table, 1, function(i) mean(i[i>=quantile(i,0.25,na.rm=T) & i<=quantile(i,0.75,na.rm=T)],na.rm=T)) , decreasing=TRUE), n=upto)
	}else{
		select <- head( order(rowMeans(my.table), decreasing=TRUE ), n=upto ) # top n FPKM genes
	}
	my.table <- as.data.frame(my.table[select,])
    if(my.filter=="ensembl_gene_id"){
		anno <- dt.ensg[ensembl_gene_id %in% rownames(my.table),.(ensembl_gene_id,hgnc_symbol,gene_biotype)]		
		#anno <- getBM(attributes = my.fields, filters = my.filter, values = rownames(my.table), mart = myMart)
		# >head(anno)
		#   chromosome_name ensembl_gene_id hgnc_symbol                                                          description   gene_biotype
		#1               17 ENSG00000161970       RPL26                 ribosomal protein L26 [Source:HGNC Symbol;Acc:10327] protein_coding
		#2               11 ENSG00000213934        HBG1                    hemoglobin, gamma A [Source:HGNC Symbol;Acc:4831] protein_coding
		#3               10 ENSG00000222414    RNU2-59P  RNA, U2 small nuclear 59, pseudogene [Source:HGNC Symbol;Acc:48552]          snRNA
		#4                3 ENSG00000223247    RNU2-64P  RNA, U2 small nuclear 64, pseudogene [Source:HGNC Symbol;Acc:48557]          snRNA
		#5                6 ENSG00000234426                                                                                         lincRNA
		anno.collapsed <- plyr::ddply(anno, my.filter, function(df) paste(df[[my.id]], collapse=';'))  # isa 'data.frame'
		colnames(anno.collapsed) <- c(my.filter,my.id)
		# to collapse (group by ensembl_gene_id) ensembl gene of multiple gene names(then concat hgnc_symbol)
		# > head(anno.collapsed)
		#  ensembl_gene_id    hgnc_symbol 
		#1 ENSG00000002726   AOC1
		#2 ENSG00000003989 SLC7A2
		my.table[,my.filter]<-rownames(my.table) # add a new column id
		my.table<-merge(my.table, anno.collapsed, by=my.filter,  all.x=TRUE) # left join to add column 'hgnc_symbol'
		rownames(my.table) <- ifelse(my.table[[my.id]]=="" | is.na(my.table[[my.id]]), my.table[[my.filter]], my.table[[my.id]]) # replace with ensembl id if no gene name
		my.table[[my.filter]]<-NULL; my.table[[my.id]]<-NULL # drop columns
		# re-order ensg 
		if(iqr){
			my.table<-my.table[order( apply(my.table, 1, function(i) mean(i[i>=quantile(i)[2] & i<=quantile(i)[4]]) ), decreasing=TRUE),] # re-order by the sum of FPKM of IQR
		}else{
			my.table<-my.table[order(rowMeans(my.table), decreasing=TRUE),] # re-order by the mean of FPKM 
		}
		return(list(fpkm=my.table, select=select[1:upto], gene=anno))
	}else{
		return(list(fpkm=my.table, select=select[1:upto]))
	}
}

# 5-set Venn Diagram
# called by bin/R/DEG/VennDiagram.R
draw.venn<-function(x,file.name){
	venn.diagram(
		x,
		filename = file.name,
		col = "black",
		fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		alpha = 0.50,
		cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
		1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
		cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		cat.cex = 1.1,
		cat.fontface = "bold",
		margin = 0.07,
		rotation.degree = -32
	)
}

#https://support.bioconductor.org/p/64693/
collapseReplicatesBugFix <- function(object, groupby, run, renameCols=TRUE) {
	if (!is.factor(groupby)) groupby <- factor(groupby)
		groupby <- droplevels(groupby)
		stopifnot(length(groupby) == ncol(object))
		sp <- split(seq(along=groupby), groupby)
		countdata <- sapply(sp, function(i) rowSums(assay(object)[,i,drop=FALSE]))
		mode(countdata) <- "integer"
		colsToKeep <- sapply(sp, `[`, 1)
		collapsed <- object[,colsToKeep]
		assay(collapsed) <- countdata
		if (!missing(run)) {
			stopifnot(length(groupby) == length(run))
			colData(collapsed)$runsCollapsed <- sapply(sp, function(i) paste(run[i],collapse=","))
		}
		if (renameCols) {
			colnames(collapsed) <- levels(groupby)
		}
		stopifnot(sum(as.numeric(assay(object))) == sum(as.numeric(assay(collapsed))))
		collapsed
}

