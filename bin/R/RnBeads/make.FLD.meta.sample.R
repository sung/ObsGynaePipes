#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

my.fld.ver="FLD.v4" # FLD.v2, FLD.v3, WGBS (target-only)
this.sets<-c(`The Eight Original Samples`="the.eight",`Low-PAPPA-SGA vs Normal-PAPPA-AGA`="sga.aga.validation",
			`Low-PAPPA-AGA vs Normal-PAPPA-AGA`="aga.low.vs.normal.pappa",`PET vs Controls`="pet.control",
			`Cord SGA vs Cord AGA`="cord.sga.aga") #, `Doppler analysis group`="doppler.group") 

pdf_file_name<-paste0("~/results/RnBeads/",my.fld.ver,".sample.cov.plot")
pdf(file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

#fld.file<-"~/results/RnBeads/Meta/FLD.v2/FLD.v2.sample.info.csv" # 141 FLD v1 samples 
fld.file<-paste0("~/results/RnBeads/Meta/",my.fld.ver,"/",my.fld.ver,".sample.info.csv") # 141 FLD v2 samples 
if(!exists(fld.file)){
	if(my.fld.ver=="WGBS"){
		fld.file<-"~/results/RnBeads/Meta/WGBS.nkx1.2.fld.ontarget.csv" # 8 WGBS samples ontarget of FLD NKX1.2 
		this.sets<-c(`The Eight Original Samples`="the.eight")
	}else if(grepl("FLD",my.fld.ver)){
		#echo "Barcode BedFile" > ~/Pipelines/data/MJ.FLD/FLD.v2/FLD.BedFile.txt
		#for j in $(for i in `ls ~/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done | uniq);do printf "$j "; echo ~/results/$SLX.Homo_sapiens.v1/BisSNP/$j/$SLX.$j.cpg.filtered.CG.bed; done >> ~/Pipelines/data/MJ.FLD/FLD.v2/FLD.BedFile.txt
		#fld.bed.file<-file.path(paste0("~/Pipelines/data/MJ.FLD",my.fld.ver,"/",my.fld.ver,".chr10.BedFile.txt"))
		fld.bed.file<-file.path(paste0("~/Pipelines/data/MJ.FLD/",my.fld.ver,"/",my.fld.ver,".BedFile.txt"))
		fld.bed<-read.table(fld.bed.file,header=T)

		fld.info.file<-file.path(paste0("~/Pipelines/data/MJ.FLD/",my.fld.ver,"/mj.FLD.RnBeads.sample.sheet.txt"))
		fld.info<-read.table(fld.info.file, header=T)
		if(my.fld.ver=="FLD.v2"){
			colnames(fld.info)<-c("SampleId","Barcode","CRN","Pair","PairGroup","Condition","DopplerGroup","DopplerCategory","PAPPA")
		}else if(my.fld.ver=="FLD.v3" | my.fld.ver=="FLD.v4"){
			colnames(fld.info)<-c("SampleId","Barcode.oxBS","Barcode.BS","CRN","Pair","PairGroup","Condition","DopplerGroup","DopplerCategory","PAPPA")
		}else{
			stop("this assay not supported")
		}
		fld.info$PAPPA<-as.numeric(!fld.info$PAPPA) # Michelle's 0/1 PAPPA reversed!

		if(my.fld.ver=="FLD.v2"){
			fld.sex.file<-file.path(paste0("~/Pipelines/data/MJ.FLD/",my.fld.ver,"/FLD .sex.csv"))
			fld.sex<-read.csv(fld.sex.file,header=T)
			fld.info<-merge(fld.info,fld.sex, all.x=TRUE) # attach 'FetalSex'
			fld.info$SampleName<-paste0(ifelse(is.na(fld.info$Condition),as.character(fld.info$DopplerCategory),as.character(fld.info$Condition)), ifelse(is.na(fld.info$Pair),"",fld.info$Pair), fld.info$FetalSex) # create 'SampleName'
			fld.info$Pair<-ifelse(is.na(fld.info$Pair),NA,paste0("Pair",fld.info$Pair))

			fld.info<-merge(fld.info,fld.bed, all.x=TRUE) # attach 'BedFile'
			write.csv(fld.info,file=fld.file, row.names=F)
		}else if(my.fld.ver=="FLD.v3" | my.fld.ver=="FLD.v4"){
			fld.pheno.file<-file.path(paste0("~/Pipelines/data/MJ.FLD/",my.fld.ver,"/variables.txt"))
			fld.pheno<-read.delim(fld.pheno.file, header=T)

			fld.info.merged<-merge(fld.info,fld.pheno, all.x=T, by=c("SampleId","Barcode.oxBS","CRN")) # by=c("SampleId","CRN","Barcode.oxBS","Barcode.BS")
			fld.info.merged$Barcode.BS<-fld.info.merged$Barcode.BS.x
			fld.info.merged$Barcode.BS.x<-NULL; fld.info.merged$Barcode.BS.y<-NULL

			fld.info.merged$FetalSex<-factor(toupper(fld.info.merged$FetalSex))
			fld.info.merged$SampleName<-paste0(ifelse(is.na(fld.info.merged$Condition),as.character(fld.info.merged$DopplerCategory),as.character(fld.info.merged$Condition)), ifelse(is.na(fld.info.merged$Pair),"",fld.info.merged$Pair), fld.info.merged$FetalSex) # create 'SampleName'
			fld.info.merged$Pair<-ifelse(is.na(fld.info.merged$Pair),NA,paste0("P",fld.info.merged$Pair))

			dummy1<-merge(subset(fld.info.merged, select=-Barcode.oxBS), fld.bed, by.x="Barcode.BS", by.y="Barcode", all.x=T); dummy1$assay<-"BS"
			dummy2<-merge(subset(fld.info.merged, select=-Barcode.BS), fld.bed, by.x="Barcode.oxBS", by.y="Barcode", all.x=T); dummy2$assay<-"oxBS"
			colnames(dummy1)<-c("Barcode",colnames(dummy1)[2:ncol(dummy1)])
			colnames(dummy2)<-c("Barcode",colnames(dummy2)[2:ncol(dummy2)])
			fld.meta.sample<-rbind(dummy1,dummy2)
			write.csv(fld.meta.sample, file=fld.file, row.names=F)
		}else{
			stop("this assay not supported")
		}
	}else{
		stop("this assay not supported")
	}
}

## Define Functions
check.paired<-function(this.samples){
	this.samples$Pair<-droplevels(this.samples$Pair) # relevel Pair
	paired<-sapply(split(this.samples, this.samples$Pair), function(i) nrow(i)==2)
	select<-this.samples$Pair %in% names(paired[paired])
	this.samples<-this.samples[select,]
	this.samples[order(this.samples$PairGroup, this.samples$Pair, this.samples$Condition),]
}
filter.by.min.cpg<-function(this.samples, this.set, this.min.cpg){
	select<-!is.na(this.samples$numCpG.d10) & this.samples$numCpG.d10>=min.cpg
	this.samples<-this.samples[select,]
	# test whether it's paired or not (SGA/AGA)
	this.samples<-check.paired(this.samples)

	my.csv<-paste0("~/results/RnBeads/Meta/",my.fld.ver,"/",my.fld.ver,".",this.set,".min.",min.cpg,".cpg.csv")
	that.samples<-with(this.samples, data.frame(BedFile, SampleId, SampleName, Pair, Condition, PairGroup, FetalSex))
	write.csv(that.samples, file=my.csv, row.names=FALSE, quote=FALSE)
	cat(my.csv," has written\n")
	return(this.samples)
}

fld.meta.sample<-read.csv(fld.file)
####################
## get CpG number ##
####################
library(plyr)
library(ggplot2)
fld.meta.sample$numCpG<-unlist( lapply(as.character(fld.meta.sample$BedFile), function(i) if(file.exists(i)){if(length(readLines(i))>1){nrow(read.table(i, skip=1))}else{0}}else{NA} ) )
fld.meta.sample$numCpG.d10<-unlist(lapply(as.character(fld.meta.sample$BedFile), function(i) if(file.exists(i)){if(length(readLines(i))>=2){nrow(read.table(i, skip=1)[read.table(i, skip=1)[,"V5"]>=10,])}else{0}}else{NA}))

# SGA Low-PAPPA in any cases
for(i in seq(1:length(this.sets))){
	min.cpg=110 # at least this number of CpGs
	# q1=original 4 SGA cases and AGA control pairs
	if(this.sets[i]=="the.eight"){
		select<-!is.na(fld.meta.sample$PairGroup) & fld.meta.sample$PairGroup=="q1" & fld.meta.sample$assay=="oxBS"
		#select<-!is.na(fld.meta.sample$Condition) & (fld.meta.sample$Condition=="SGA" | fld.meta.sample$Condition=="AGA") & fld.meta.sample$PairGroup=='q1'
	# q2=additional 15 SGA cases and AGA control pairs
	}else if(this.sets[i]=="sga.aga.validation"){
		select<-!is.na(fld.meta.sample$PairGroup) & fld.meta.sample$PairGroup=="q2" & fld.meta.sample$assay=="oxBS"
		#select<-!is.na(fld.meta.sample$Condition) & ((fld.meta.sample$Condition=="SGA" & fld.meta.sample$PAPPA==0) | (fld.meta.sample$Condition=="AGA" & fld.meta.sample$PAPPA==1)) & fld.meta.sample$PairGroup!='q6'
	# q3=10 AGA with low PAPP-A cases and AGA with normal PAPP-A controls pairs
	}else if(this.sets[i]=="aga.low.vs.normal.pappa"){
		select<-!is.na(fld.meta.sample$PairGroup) & fld.meta.sample$PairGroup=="q3" & fld.meta.sample$assay=="oxBS"
	# q4=8 severe pre-eclampsia and normal control pairs
	}else if(this.sets[i]=="pet.control"){
		select<-!is.na(fld.meta.sample$PairGroup) & fld.meta.sample$PairGroup=="q4" & fld.meta.sample$assay=="oxBS"
	# q6=14 cord DNA SGA cases and cord DNA AGA control pairs
	}else if(this.sets[i]=="cord.sga.aga"){
		select<-!is.na(fld.meta.sample$PairGroup) & fld.meta.sample$PairGroup=="q6"
	}else if(this.sets[i]=="doppler.group"){
		select<-!is.na(fld.meta.sample$DopplerGroup) & fld.meta.sample$DopplerGroup=="q5"
	}else{
		stop(paste0(this.sets[i], " not supported. Stoppped!"))
	}

	filter.fld.meta.sample<-fld.meta.sample[select,] # apply filter
	if(this.sets[i]!="cord.sga.aga"){filter.fld.meta.sample<-check.paired(filter.fld.meta.sample)} # test whether it's paired or not (SGA/AGA)
	filter.fld.meta.sample<-filter.by.min.cpg(filter.fld.meta.sample,this.sets[i],min.cpg)

	##
	## CpG Coverage by Depth and the number of supporting Samples
	merged.bed<-lapply(as.character(filter.fld.meta.sample$BedFile), function(i){x<-read.delim(i, skip=1, header=F); x$Barcode<-unlist(strsplit(i,'/'))[7]; return(x)})
	merged.bed<-data.table::rbindlist(merged.bed)
	merged.bed$ID<-paste(merged.bed$V1,merged.bed$V3,sep=".")

	my.table<-data.frame()
	for(min.depth in seq(1:10)){
		print(paste0("min.depth=",min.depth,"\n"))
		dummy<-plyr::ddply(merged.bed[merged.bed$V5>=min.depth,], c("ID"), function(df) length(unique(df$Barcode))) # group by ID (chr.position), V5: depth of coverage
		for(sample.cnt in rev(seq(1:nrow(filter.fld.meta.sample)))){
			my.table<-rbind(my.table, c(min.depth,sample.cnt,nrow(dummy[dummy$V1>=sample.cnt,])))
		}	
	}
	colnames(my.table)<-c("min.cov","num.supporting.sample","num.cpg")
	my.table$min.cov<-as.factor(my.table$min.cov)

	# coverage supported by the number of samples
	my.csv<-paste0("~/results/RnBeads/Meta/",my.fld.ver,"/",my.fld.ver,".",this.sets[i],".min.",min.cpg,".cpg.cov.stat.csv")
	write.csv(my.table, file=my.csv, row.names=FALSE, quote=FALSE)
	cat(my.csv," has written\n")

	hmcol <- colorRampPalette(c("black","royalblue"))(length(levels(my.table$min.cov)))

	print(
		ggplot(my.table, aes(x=reorder(my.table$num.supporting.sample, num.cpg), y=num.cpg)) + 
		geom_point() + 
		geom_line(aes(colour=min.cov, group=min.cov)) + 
		scale_colour_manual(values=hmcol) + 
		xlab("Number of Supporting Samples") + 
		ylab("Number of CpG") + 
		ggtitle(paste0(names(this.sets[i]), " (n=", nrow(filter.fld.meta.sample),")"))
	)
}

dev.off()
cat("All is done", "\n")
