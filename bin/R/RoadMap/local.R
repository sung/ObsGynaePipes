library(data.table)
################################
## RoadMap                    ##
## Schultz et al. Nature 2015 ##
## Tissue comparision         ##
################################

#for i in `ls ~/data/RoadMap/BS-Seq/Schultz2015/*.gz`; do echo $i | awk "BEGIN{FS=\"/\";OFS=\",\"}{split(\$8,a,\".tsv.gz\"); split(a[1],b,\" <- \"); print \$i,b[2],\"individual\"b[3],\"chr\"b[4]}"; done > ~/data/RoadMap/BS-Seq/Schultz2015.tissue.meta.csv
#cnt=0; for i in `ls ~/data/RoadMap/BS-Seq/Schultz2015/*.gz`; do let cnt++; if [ $cnt -eq 1 ];then printf "File,Tissue,Ind,Chr\n"; fi; echo $i | awk "BEGIN{FS=\"/\";OFS=\",\"}{split(\$8,a,\".tsv.gz\"); split(a[1],b,\"_\"); print \$i,b[2],\"individual\"b[3],\"chr\"b[4]}"; done  >  ~/data/RoadMap/BS-Seq/Schultz2015.tissue.meta.csv
schultz.meta<-read.csv("~/data/RoadMap/BS-Seq/Schultz2015.tissue.meta.csv", stringsAsFactor=FALSE)

# select tissue having both individual2 and individual3 only 
select<-schultz.meta$Ind=='individual2' | schultz.meta$Ind=='individual3'
dummy<-schultz.meta[select,]
tissue.filter<-sapply(split(dummy, dummy$Tissue), function(i) length(unique(i$Ind)))
avail.tissues<-names(tissue.filter[tissue.filter==2]) # 9 tissues ("AD" "AO" "EG" "FT" "GA" "PA" "PO" "SB" "SX") where both STL2(female) and STL3(male) are available
													  # PA: not imported yet
select<-dummy$Tissue %in% avail.tissues

#schultz.list<-lapply(split(schultz.meta, schultz.meta$Ind), function(i) split(i, i$Tissue)) # isa 'list'
schultz.list<-lapply(split(dummy[select,], dummy[select,]$Ind), function(i) split(i, i$Tissue)) # isa 'list'

#for i in `ls ~/data/RoadMap/BS-Seq/Schultz2015/*.gz`; do echo $i | cut -d/ -f 8 | cut -d. -f1 | cut -d'_' -f2; done | sort | uniq 
tissue.list=list(
	AD='Adrenal',
	AO='Aorta',
	BL='Bladder',
	EG='Oesophagus',
	FT='Fat',
	GA='Gastric',
	LG='Lung',
	LI='Liver',
	LV='Left Ventricle',
	OV='Ovary',
	PA='Pancreas',
	PO='Psoas',
	RA='Right Atrium',
	RV='Right Ventricle',
	SB='Small Bowel',
	SG='Sigmoid Colon',
	SX='Spleen',
	TH='Thymus',
	PT='Placenta'
	)
##########
# Config #
##########
#my.classes=c("character","integer","character","character","integer","integer","integer")
min.doc <- 10 # depth of coverage
num.cores <-4  # Set Multi-cores

#options(scipen=999) # diable scientific notation
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM

# called by process.roadmap.tissues.male.female.R
prep.roadmap<-function(my.tissue,is.filter=FALSE, my.cpg.type){
	cat(paste0("Processing ", my.tissue, "...\n"))
	# data will be saved in this RData file
	my.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste0(my.tissue,".ind2.ind3.RData"))

	#####################################################################################
	# 1. non-parallel version to avlid "long vectors not supported yet: memory.c:3323"  #
	# IT WORKS BUT SLOW
	#     user   system  elapsed 
	# 2590.541  350.311 2955.665 
	#####################################################################################
	#system.time(female<-lapply(schultz.list$individual2[[my.tissue]]$File, read.delim, header=F, comment.char="", colClasses=my.classes)) # isa 'list'
	#system.time(male<-lapply(schultz.list$individual3[[my.tissue]]$File, read.delim, header=F, comment.char="", colClasses=my.classes)) # isa 'list'

	##################################################################################################################
	# 2. parallel version: read the methylation file                                                                ##
	# (NB. mc.preschedule=F) http://stackoverflow.com/questions/23231183/mclapply-long-vectors-not-supported-yet    ##
	# FAILED: "long vectors not supported yet: memory.c:3323"                                                       ##
	##################################################################################################################
	# http://www.biostat.jhsph.edu/~rpeng/docs/R-large-tables.html
	#system.time(female<-parallel::mclapply(schultz.list$individual2[[my.tissue]]$File, read.delim, header=F, comment.char="", colClasses=my.classes, mc.cores=num.cores, mc.preschedule=F)) # isa 'list'
	#system.time(male<-parallel::mclapply(schultz.list$individual3[[my.tissue]]$File, read.delim, header=F, comment.char="", colClasses=my.classes, mc.cores=num.cores, mc.preschedule=F)) # isa 'list'

	##########################################################################################
	## 3. foreach way
	# http://www.r-bloggers.com/import-all-text-files-in-a-folder-with-parallel-execution/  ##
	# FAILED: long vectors not supported yet: memory.c:3323
	##########################################################################################
	#library(doParallel)
	#library(foreach)
	#registerDoParallel(cores = num.cores)
	#system.time( female <- foreach(i = schultz.list$individual2[[my.tissue]]$File) %dopar% read.delim(i, header=F, comment.char="", colClasses=my.classes) )
	#system.time( male <- foreach(i = schultz.list$individual3[[my.tissue]]$File) %dopar% read.delim(i, header=F, comment.char="", colClasses=my.classes) )

	# filter CpG of by coverage
	# i: isa 'data.frame'
	#X       2699425 +       CCC     0       1       0
	#X       2699426 +       CCG     0       1       0
	#X       2699427 +       CGG     2       2       0
	#i$V7==1: CpG only
	#i$V6>=10: at least 10x coverage
	#female<-parallel::mclapply(female, function(i){select<-i$V7==1 & i$V6>=min.doc; return(i[select,])}, mc.cores=num.cores) # isa 'list'
	#male<-parallel::mclapply(male, function(i){select<-i$V7==1 & i$V6>=min.doc; return(i[select,])}, mc.cores=num.cores) # isa 'list'

	#female.cpg<-Reduce(function(...) rbind(...), female) # combine all chr
	#male.cpg<-Reduce(function(...) rbind(...), male) # combine all chr

	#female.cpg<-do.call(rbind, female) # isa 'data.frame'
	#male.cpg<-do.call(rbind, male) # isa 'data.frame'

	# merge by: Chr, Position, Strand
	#merged.cpg<-merge(female.cpg, male.cpg, by=c("V1","V2","V3"))

	#> head(merged.cpg)
	#  V1        V2 V3 V4.x V5.x V6.x V7.x V4.y V5.y V6.y V7.y
	#1  1  10000005  +  CGT   14   14    1  CGT   17   19    1
	#2  1  10000006  -  CGG   19   19    1  CGG   26   28    1
	#3  1 100000180  +  CGA   15   25    1  CGA   46   91    1

	############################################################
	## 4. data.table way (parallel)                            #
	# https://github.com/ActiveAnalytics/Parallel-R-File-Read  #
	# FAILED: long vectors not supported yet: memory.c:3323
	############################################################
	#library(data.table)
	#system.time(female<-parallel::mclapply(paste('zcat', schultz.list$individual2[[my.tissue]]$File), fread, mc.cores=num.cores)) # isa 'list' of data.table
	#system.time(male<-parallel::mclapply(paste('zcat', schultz.list$individual3[[my.tissue]]$File), fread, mc.cores=num.cores)) # isa 'list' of data.table

	###############################################################
	## 5. data.table way (un-parallel)                            #
	# https://github.com/ActiveAnalytics/Parallel-R-File-Read     #
	# BEST SPEED!
	#   user  system elapsed
	#534.251  37.097 578.091
	###############################################################
	print(system.time(female<-lapply(paste('zcat', schultz.list$individual2[[my.tissue]]$File), fread))) # isa 'list' of data.table (by chr)
	print(system.time(male<-lapply(paste('zcat', schultz.list$individual3[[my.tissue]]$File), fread))) # isa 'list' of data.table (by chr)

	# fread & filer
	#system.time(female<-lapply(paste('zcat', schultz.list$individual2[[my.tissue]]$File), function(i){DT=fread(i); return(DT[V7==1 & V6>=min.doc,,]) })) # isa 'list' of data.table
	#system.time(male<-lapply(paste('zcat', schultz.list$individual3[[my.tissue]]$File), function(i){DT=fread(i); return(DT[V7==1 & V6>=min.doc,,]) })) # isa 'list' of data.table

	# name list by chr
	names(female)<-schultz.list$individual2[[my.tissue]]$Chr
	names(male)<-schultz.list$individual3[[my.tissue]]$Chr
	assign(my.tissue, list(`STL2`=female,`STL3`=male))

	###########
	# Save it #
	###########
	save(list=my.tissue, file=my.RData)
	cat(paste0(my.tissue, " saved\n"))

	##########################################
	# Further filter by depth and only CpG  ##
	##########################################
	#V1: chr
	#V2: position
	#V3: str
	#V4: triplet-cpg-context (e.g. CTT)
	#V5: number of C (methylated)
	#V6: number of C+T (methylated + unmethylated)
	#V7: flag for more-than-noise

	#10      101081  -       CCA     0       10      0
	#10      101082  +       CGG     6       7       1
	#10      101083  -       CGC     7       10      1
	#10      101084  -       CCG     0       11      0
	#10      101086  +       CGA     10      11      1
	#10      101087  -       CGT     10      12      1
	if(is.filter){
		# data will be saved into this file
		my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,my.cpg.type,"dt.merged.RData",sep="."))
		cat(paste0("Filtering ", my.tissue, "...\n"))

		# for each chr
		if(my.cpg.type=="CG"){
			system.time(dt.merged.list<-lapply(names(female), function(i){
													merge(
														female[[i]][substr(V4, 1,2)=="CG" & V6>=min.doc], 
														male[[i]][substr(V4, 1,2)=="CG" & V6>=min.doc], 
														by=c("V1","V2","V3","V4")
													)
												}
				) # isa 'list' of data.table (by chr)
			)
		}else if(my.cpg.type=="CHH"){
			system.time(dt.merged.list<-lapply(names(female), function(i){
													merge(
														female[[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)!="G" & V6>=min.doc], 
														male[[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)!="G" & V6>=min.doc], 
														by=c("V1","V2","V3","V4")
													)
												}
				) # isa 'list' of data.table (by chr)
			)
		}else if(my.cpg.type=="CHG"){
			system.time(dt.merged.list<-lapply(names(female), function(i){
													merge(
														female[[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)=="G" & V6>=min.doc], 
														male[[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)=="G" & V6>=min.doc], 
														by=c("V1","V2","V3","V4")
													)
												}
				) # isa 'list' of data.table (by chr)
			)
		}else{
			stop(paste0(my.cpg.type, " Not supported\n"))
		}

		names(dt.merged.list)<-names(female)
		assign(paste0(my.tissue,".dt.merged"), dt.merged.list)
		###########
		# Save it #
		###########
		save(list=paste0(my.tissue,".dt.merged"), file=my.new.RData)
		cat(paste0(my.tissue, ".dt.merged saved\n"))
		rm(dt.merged.list)
	}
	rm(female,male) # free the memory
	cat(paste0("Done for ", my.tissue, "\n"))
}#end of function 'prep.roadmap'

# prep.roadmap must have been done previously before this
# this function does the post-processing if 'is.filter'==1 for prep.roadmap
# called by process.roadmap.tissues.dt.merged.R
dt.merge.roadmap<-function(my.tissue, my.cpg.type){
	my.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste0(my.tissue,".ind2.ind3.RData"))
	if(!file.exists(my.RData)){
		stop(paste0(my.RData," not found\n"))
	}else{
		cat(paste0("Loading ", my.RData, "...\n"))

		## Load 
		print(system.time(load(my.RData))) # load list of data.table
		cat(paste0(my.tissue, " loaded\n"))
		this.tissue<-get(my.tissue) # isa list of data.table (STL2, STL3)
		eval(parse(text=paste('rm(', my.tissue, ')'))) # remove the original object from memory

		####################
		#V1: chr
		#V2: position
		#V3: str
		#V4: triplet-cpg-context (e.g. CTT)
		#V5: number of C (methylated)
		#V6: number of C+T (methylated + unmethylated)
		#V7: flag for more-than-noise
		####################
		# data will be saved into this file
		if(my.cpg.type=="all"){
			# 1. CG
			my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,"CG.dt.merged.RData",sep="."))
			cat(paste0("Filtering ", my.tissue, " for CG...\n"))
			# for each chr
			print(system.time(dt.merged.list<-lapply(names(this.tissue[["STL2"]]), function(i){
													merge(
														this.tissue[["STL2"]][[i]][substr(V4, 1,2)=="CG"], 
														this.tissue[["STL3"]][[i]][substr(V4, 1,2)=="CG"], 
														by=c("V1","V2","V3","V4")
													)
												}
				) # isa 'list' of data.table (by chr)
			))
			names(dt.merged.list)<-names(this.tissue[["STL2"]])
			assign(paste0(my.tissue,".CG.dt.merged"), dt.merged.list)
			save(list=paste0(my.tissue,".CG.dt.merged"), file=my.new.RData)
			rm(dt.merged.list)
			eval(parse(text=paste0('rm(', my.tissue, '.CG.dt.merged)'))) # remove the original object from memory

			# 2. CHH
			my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,"CHH.dt.merged.RData",sep="."))
			cat(paste0("Filtering ", my.tissue, " for CHH...\n"))
			# for each chr
			print(system.time(dt.merged.list<-lapply(names(this.tissue[["STL2"]]), function(i){
													merge(
														this.tissue[["STL2"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)!="G" & V6>=min.doc], 
														this.tissue[["STL3"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)!="G" & V6>=min.doc], 
														by=c("V1","V2","V3","V4")
													)
												}
				) # isa 'list' of data.table (by chr)
			))
			names(dt.merged.list)<-names(this.tissue[["STL2"]])
			assign(paste0(my.tissue,".CHH.dt.merged"), dt.merged.list)
			save(list=paste0(my.tissue,".CHH.dt.merged"), file=my.new.RData)
			rm(dt.merged.list)
			eval(parse(text=paste0('rm(', my.tissue, '.CHH.dt.merged)'))) # remove the original object from memory

			# 3. CHG
			my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,"CHG.dt.merged.RData",sep="."))
			cat(paste0("Filtering ", my.tissue, " for CHG...\n"))
			# for each chr
			print(system.time(dt.merged.list<-lapply(names(this.tissue[["STL2"]]), function(i){
													merge(
														this.tissue[["STL2"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)=="G" & V6>=min.doc], 
														this.tissue[["STL3"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)=="G" & V6>=min.doc], 
														by=c("V1","V2","V3","V4")
													)
												}
				) # isa 'list' of data.table (by chr)
			))
			names(dt.merged.list)<-names(this.tissue[["STL2"]])
			assign(paste0(my.tissue,".CHG.dt.merged"), dt.merged.list)
			save(list=paste0(my.tissue,".CHG.dt.merged"), file=my.new.RData)
			rm(dt.merged.list)
			eval(parse(text=paste0('rm(', my.tissue, '.CHG.dt.merged)'))) # remove the original object from memory
		}else{
			my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,my.cpg.type,"dt.merged.RData",sep="."))
			cat(paste0("Filtering ", my.tissue, "...\n"))
			# for each chr
			if(my.cpg.type=="CG"){
				print(system.time(dt.merged.list<-lapply(names(this.tissue[["STL2"]]), function(i){
														merge(
															this.tissue[["STL2"]][[i]][substr(V4, 1,2)=="CG"], 
															this.tissue[["STL3"]][[i]][substr(V4, 1,2)=="CG"], 
															by=c("V1","V2","V3","V4")
														)
													}
					) # isa 'list' of data.table (by chr)
				))
			}else if(my.cpg.type=="CHH"){
				print(system.time(dt.merged.list<-lapply(names(this.tissue[["STL2"]]), function(i){
														merge(
															this.tissue[["STL2"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)!="G" & V6>=min.doc], 
															this.tissue[["STL3"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)!="G" & V6>=min.doc], 
															by=c("V1","V2","V3","V4")
														)
													}
					) # isa 'list' of data.table (by chr)
				))
			}else if(my.cpg.type=="CHG"){
				print(system.time(dt.merged.list<-lapply(names(this.tissue[["STL2"]]), function(i){
														merge(
															this.tissue[["STL2"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)=="G" & V6>=min.doc], 
															this.tissue[["STL3"]][[i]][substr(V4, 1,2)!="CG" & substr(V4, 3,3)=="G" & V6>=min.doc], 
															by=c("V1","V2","V3","V4")
														)
													}
					) # isa 'list' of data.table (by chr)
				))
			}else{
				stop(paste0(my.cpg.type, " Not supported\n"))
			}
			names(dt.merged.list)<-names(this.tissue[["STL2"]])
			assign(paste0(my.tissue,my.cpg.type,".dt.merged"), dt.merged.list) # SX.CG.dt.merged
			save(list=paste0(my.tissue,my.cpg.type,".dt.merged"), file=my.new.RData)
			rm(dt.merged.list)
			eval(parse(text=paste0('rm(', my.tissue,'.',my.cpg.type,'.dt.merged)'))) # remove the original object from memory

		}
		rm(this.tissue)
	}
}# end of 'dt.merge.roadmap'

# called by 'bin/R/RoadMap/local.R'
# called by 'bin/R/RoadMap/find.overlap.promoter.cpgi.R'
# called by 'bin/R/RoadMap/get.meth.level.roadmap.gender.by.chr.tissue.context.R'
# called by 'bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.chr.tissue.context.R'
# called by 'bin/R/RoadMap/plot.roadmap.gender.meth.distribution.by.chr.tissue.context.R'
# called by 'bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.XCI.Carrel.science.2005.R'
get.region<-function(i){
	if(i=="Genes"){
		#RnBeads::rnb.get.annotation(type="genes", assembly = "hg19") # isa 'GRangesList' splitted by chr
		subject<-Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
	}else if(i=="miR_Genes"){
		subject<-rtracklayer::import.gff("~/data/Annotation/miRbase.20/hsa.gff3")
	#+100b, -100b
	}else if(i=="miR_Genes_01.01"){
		subject<-rtracklayer::import.gff("~/data/Annotation/miRbase.20/hsa.gff3")
		ranges(subject)=IRanges(start=start(subject)-100,end=end(subject)+100,names=names(subject))
	#+1kb, -1kb
	}else if(i=="miR_Genes_10.10"){
		subject<-rtracklayer::import.gff("~/data/Annotation/miRbase.20/hsa.gff3")
		ranges(subject)=IRanges(start=start(subject)-1000,end=end(subject)+1000,names=names(subject))
	#+1kb, -1kb
	}else if(i=="Genes_10.10"){
		subject=Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
		ranges(subject)=IRanges(start=start(subject)-1000,end=end(subject)+1000,names=names(subject))
	#+2kb, -2kb
	}else if(i=="Genes_20.20"){
		subject=Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
		ranges(subject)=IRanges(start=start(subject)-2000,end=end(subject)+2000,names=names(subject))
		#below takes too long
		#subject<-unlist(GRangesList(lapply(gr.genes, function(i) promoters(i, upstream=2000, downstream=end(i)-start(i)+1+2000)))) #includes 1000 +TSS and 1000 from +TES 
	#+5kb, -5kb
	}else if(i=="Genes_50.50"){
		subject=Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
		ranges(subject)=IRanges(start=start(subject)-5000,end=end(subject)+5000,names=names(subject))
	# -1.5kb, +.5kb
	}else if(i=="Promo_15.05"){
		subject<-Reduce(c,RnBeads::rnb.get.annotation(type="promoters", assembly = "hg19"))
	# -1.5kb, +1.5kb
	}else if(i=="Promo_15.15"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1500, downstream=1500)
	# -1.5kb, +1.0kb
	}else if(i=="Promo_15.10"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1500, downstream=1000)
	# -1.5kb, +0kb
	}else if(i=="Promo_15.00"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1500, downstream=0)
	# -1.0kb, +1.5kb
	}else if(i=="Promo_10.15"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=1500)
	# -1.0kb, +1.0kb
	}else if(i=="Promo_10.10"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=1000)
	# -1.0kb, +.5kb
	}else if(i=="Promo_10.05"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=500) 
	# -1.0kb, +.3kb
	}else if(i=="Promo_10.03"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=300) # Promo_10.03
	# -1.0kb, +.0kb
	}else if(i=="Promo_10.00"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=0)
	# -2.5kb, +.0kb
	}else if(i=="Promo_25.00"){
		subject<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=2500, downstream=0)
	}else if(i=="miR_Promo"){
		subject<-rtracklayer::import.bed("~/data/Annotation/PROmiRNA/PROmiRNA.grch37.bed", genome="hg19")
	# -100b, +100b
	}else if(i=="miR_Promo_01.01"){
		subject<-rtracklayer::import.bed("~/data/Annotation/PROmiRNA/PROmiRNA.grch37.bed", genome="hg19")
		ranges(subject)=IRanges(start=start(subject)-100,end=end(subject)+100,names=names(subject))
	}else if(i=="Tiling5K"){
		subject<-Reduce(c,RnBeads::rnb.get.annotation(type="tiling", assembly = "hg19"))
	}else if(i=="Tiling500"){
		#method1 via RnBeads
		#RnBeads::rnb.load.annotation("~/data/Annotation/RnBeads/hg19.genome.tiling500.rnb.RData", i)
		#print(system.time(subject<-unlist(RnBeads::rnb.get.annotation(type=i, assembly = "hg19"))))
		#print(system.time(subject<-Reduce(c,RnBeads::rnb.get.annotation(type=i, assembly = "hg19")))) #same as above, but take much longer

		#method2 via GenomicRanges
		print(system.time(subject<-GenomicRanges::tileGenome(seqlengths(gr.ensg)[1:24], tilewidth=500, cut.last.tile.in.chrom=TRUE))) # the size of first tile always 'tilewidth'
		#print(system.time(subject<-unlist(GenomicRanges::tileGenome(seqlengths(gr.ensg)[1:24], tilewidth=500)))) ' the size of first tile not 'tilewidth'
	}else if(i=="CPGi"){
		subject<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))
	}else if(i=="CPGi-with-promo"){
		gr.first<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))
		gr.second<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=1000) # Promo_10.10
		#gr.second<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=300) # Promo_10.03
		# intersection betweeen the two regions
		subject<-intersect(gr.first, gr.second, ignore.strand=TRUE)
	}else if(i=="CPGi-no-promo"){
		gr.first<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))
		gr.second<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=1000) # Promo_10.10
		#gr.second<-GenomicRanges::promoters(Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19")), upstream=1000, downstream=300) # Promo_10.03
		# setdiff betweeen the two regions
		subject<-setdiff(gr.first, gr.second, ignore.strand=TRUE)
	}else if(i=="CPGi-with-gene"){
		gr.first<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))
		gr.second<-Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
		# intersection betweeen the two regions
		subject<-intersect(gr.first, gr.second, ignore.strand=TRUE)
	}else if(i=="CPGi-no-gene"){
		gr.first<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))
		gr.second<-Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
		# setdiff betweeen the two regions
		subject<-setdiff(gr.first, gr.second, ignore.strand=TRUE)
	}else if(i=="CPGi-shores"){
		my.grange=Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))
		#https://github.com/al2na/methylKit/blob/master/R/annotate.R
		my.shores=c( IRanges::flank(my.grange,2000),IRanges::flank(my.grange,2000,FALSE) )
		subject=IRanges::reduce(IRanges::setdiff(my.shores, my.grange)) # ,erge overlapping shores remove CpG coordinates from all shores, das ist so cool!!
	}else if(i=="CPGi_Genes_5K"){
		# Genes_50.50 (5K up/down)
		gr.first=Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
		ranges(gr.first)=IRanges(start=start(gr.first)-5000,end=end(gr.first)+5000,names=names(gr.first))
		gr.first$ensg<-names(gr.first) # add ensg

		# CPGi
		gr.cpgi<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))

		# intersection betweeen the two regions
		gr.second<-intersect(gr.first, gr.cpgi, ignore.strand=TRUE)

		# get the ensg
		gr.overlap<-mergeByOverlaps(gr.first, gr.second) # isa 'DataFrame' (any overlap size>0)

		# attach gene info to the 'intersection'
		subject<-makeGRangesFromDataFrame(
					with(gr.overlap, data.frame(Chromosome=seqnames(gr.overlap$gr.second), Start=start(gr.overlap$gr.second), End=end(gr.overlap$gr.second), Strand=strand(gr.overlap$gr.second), ensg=ensg, Strand.gene=strand(gr.overlap$gr.first))),
					keep.extra.columns=TRUE
					)
	}else if(i=="CPGi_Genes_2K"){
		# Genes_20.20 (2K up/down)
		gr.first=Reduce(c,RnBeads::rnb.get.annotation(type="genes", assembly = "hg19"))
		ranges(gr.first)=IRanges(start=start(gr.first)-2000,end=end(gr.first)+2000,names=names(gr.first))
		gr.first$ensg<-names(gr.first) # add ensg

		# CPGi
		gr.cpgi<-Reduce(c,RnBeads::rnb.get.annotation(type="cpgislands", assembly = "hg19"))

		# intersection betweeen the two regions
		gr.second<-intersect(gr.first, gr.cpgi, ignore.strand=TRUE)

		# get the ensg
		gr.overlap<-mergeByOverlaps(gr.first, gr.second) # isa 'DataFrame' (any overlap size>0)

		# attach gene info to the 'intersection'
		subject<-makeGRangesFromDataFrame(
					with(gr.overlap, data.frame(Chromosome=seqnames(gr.overlap$gr.second), Start=start(gr.overlap$gr.second), End=end(gr.overlap$gr.second), Strand=strand(gr.overlap$gr.second), ensg=ensg, Strand.gene=strand(gr.overlap$gr.first))),
					keep.extra.columns=TRUE
					)
	}else if(i=="Enhancer"){
		subject<-rtracklayer::import.bed("~/data/Annotation/FANTOM5/permissive_enhancers.bed", genome="hg19")
	}else if(i=="Enhancer-ts"){
		subject<-Reduce(c,rtracklayer::import.bed("~/data/Annotation/FANTOM5/enhancer_tissue_specific.bed", genome="hg19"))
	}else if(i=="ref.genes"){
		subject<-gr.ensg # from Annotation.R
	}else if(i=="ref.promoter1k"){
		subject<-GenomicRanges::promoters(gr.ensg, upstream=1000, downstream=300)
	}else if(i=="ref.cpgislands"){
		gl.cpg.obj=methylKit::read.feature.flank(location=my.cgi.file, remove.unsual=TRUE, flank=2000, feature.flank.name=c(i,"ref.cpgisland.shores"))
		subject<-gl.cpg.obj[[i]]
	}else if(i=="ref.cpgisland.shores"){
		gl.cpg.obj=methylKit::read.feature.flank(location=my.cgi.file, remove.unsual=TRUE, flank=2000, feature.flank.name=c("ref.cpgislands",i))
		subject<-gl.cpg.obj[[i]] # gl.cpg.obj[[i]] isa GRanges for 'ref.cpgisland.shores'
	}else{
		stop(paste0(i, " Not supported\n"))
	}
	return(subject)
}

# called by plot.roadmap.gender.meth.diff.by.XCI.Carrel.science.2005.R
# input1: my.tissue (e.g. "PT")
# input2: my.xci (e.g. "XIE") xci.list[[my.xci]]$ensembl_gene_id
# input3: my.context (e.g. "Genes")
# output: data.table of meth.diff (female-male) of this cpg region and this tissue type
# NOTE1: dt.query should have been initialised
# NOTE2: xci.list should have been initialised
get.meth.for.xci<-function(my.tissue, my.xci, my.context, with.ensg=FALSE){
	cat(paste0("Processing ", my.context,"...\n"))

	# strand-aware (also "ensembl_gene_id" available by default)
	if(grepl('^Gene',my.context) || grepl('^Promo',my.context)){
		my.region<-as.data.frame(get.region(my.context)) # isa 'data.frame'
		my.region$ensg_id=rownames(my.region) # isa 'data.frame'
		# chrX-wide
		if(my.xci=="chrX"){
			select<-my.region$seqnames=="chrX"
		# selected genes from Carrel et al.
		}else{
			select<-rownames(my.region) %in% xci.list[[my.xci]][,ensembl_gene_id]
		}
		if(with.ensg){
			subject<-with(my.region[select,], data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand, ENSG=ensg_id))
			group.keys=c("ENSG")               # V3: strand
			group.keys.list="ENSG"                   # same as above 
		}else{
			subject<-with(my.region[select,], data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand))
			group.keys=c("V1","V3","Start","End")               # V3: strand
			group.keys.list="V1,V3,Start,End"                   # same as above 
		}
		#rownames(subject)<-rownames(my.region[select,])

		subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
		query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)

	# non-strand-aware (Tiling5K, CPGi, CPGi-shores, CPGi_Genes_5K, Enhancer, etc)
	}else{
		my.region<-get.region(my.context) # isa 'GRanges' (no ensg)
		# chrX-wide
		if(my.xci=="chrX"){
			select<-seqnames(my.region)=="chrX"
			if(with.ensg){
				subject<-with(as.data.frame(my.region[select,]), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand, ENSG=ensg)) # for CPGi_Genes_5K
			}else{
				subject<-with(as.data.frame(my.region[select,]), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand))
			}
		# XCI genes
		}else{
			if(my.context=="Tiling5K"){
				# start with the gene-body regions of my.xci genes
				gr.genes<-get.region("Genes") # isa 'GRanges'
			#CPGi &  Enhancer
			}else{
				# start with extended gene-body regions of my.xci genes
				gr.genes<-get.region("Genes_20.20") # isa 'GRanges'
				#gr.genes<-get.region("Genes_50.50") # isa 'GRanges' # 1/Aug/2016
			}
			select<-names(gr.genes) %in% xci.list[[my.xci]][,ensembl_gene_id]
			my.overlap=as.matrix( findOverlaps(gr.genes[select,], my.region, ignore.strand=TRUE) ) #By default, any overlap is accepted
			subject<-with(my.region[my.overlap[,"subjectHits"]], data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand))
			rm(gr.genes,my.overlap)
		}
		# cpgislands, tiling, genomeTiling500
		subject.keys=c("Chromosome","Start","End")
		query.key=c("V1","V2","End")
		if(with.ensg){
			group.keys=c("ENSG")
			group.keys.list="ENSG"
		}else{
			group.keys=c("V1","Strand","Start","End")
			group.keys.list="V1,Strand,Start,End"
		}
	}

	cat("\tRemoving leading 'chr'...\n")
	#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
	subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X, chr13=>13
	cat("\tconvert to data.table...\n")
	dt.subject<-as.data.table(subject) # isa data.table
	rm(my.region, subject)

	setkeyv(dt.subject,subject.keys)

	cat("\tFinding CpGs overlapping with ",my.context,"...\n")
	print(system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L)))
	#> dt.overlap
	#   V1    Start      End Strand       V2 V3  V4 V7 V5.x V6.x V5.y V6.y    i.End
	#1:  X  2746554  2748235      *  2746590  - CGA  1    6   17    2   10  2746590
	#2:  X  2835950  2836236      *  2836235  + CGG  1   14   14    9   10  2836235

	cat("\tAggregating %met for this region...\n")
	if(with.ensg){
		dt.cpg.region<-dt.overlap[,list(f.met=sum(V5.x)/sum(V6.x), m.met=sum(V5.y)/sum(V6.y), diff.met=sum(V5.x)/sum(V6.x)-sum(V5.y)/sum(V6.y), num.sites=.N, xci.type=my.xci, tissue=my.tissue, cpg.region=my.context),by=group.keys.list]
		setnames(dt.cpg.region, c("ensembl_gene_id","f.met", "m.met", "diff.met","num.sites","XCI","tissue","cpg.region"))
	}else{
		dt.cpg.region<-dt.overlap[,list(f.met=sum(V5.x)/sum(V6.x), m.met=sum(V5.y)/sum(V6.y), diff.met=sum(V5.x)/sum(V6.x)-sum(V5.y)/sum(V6.y), num.sites=.N, xci.type=my.xci, tissue=my.tissue, cpg.region=my.context),by=group.keys.list]
		setnames(dt.cpg.region, c("chr","strand","start","end","f.met", "m.met","diff.met","num.sites","XCI","tissue","cpg.region"))
	}
	#> dt.cpg.region
	#   chr strand    start      end    diff.met num.sites XCI tissue cpg.region
	#1:   X      *  2746554  2748235  0.09930111         6 XIE     SB       CPGi
	#2:   X      *  2835950  2836236  0.10000000         1 XIE     SB       CPGi
	return(dt.cpg.region)
} #end of get.meth.for.xci

############################################################
## JOB2: make RData for chrX genes at various CPG contexts #
## This stores 'dt.meth.chrX.per.tissue'
## This is for the preprocessing of JOB3 below
## dt.meth.chrX.per.tissue.CG.genes.20.20.RData
## JOB3: in bin/R/RoadMap/plot.roadmap.gender.meth.diff.by.XCI.Carrel.science.2005.R
############################################################
prep.dt.meth.chrX.per.tissue<-function(my.cpg.type="CG", my.cpg.contexts=c("Genes","Promo_15.05", "Promo_10.10","CPGi","CPGi-shores","Enhancer"), gene_size="Genes_20.20"){
	cat(paste("C-context:",my.cpg.type,"...\n"))
	cat(paste("Gene size:",gene_size,"...\n"))
	gr.genes<-get.region(gene_size)
	mcols(gr.genes)<-DataFrame(mcols(gr.genes), ensembl_gene_id=names(gr.genes)) # add ensembl_gene_id as a column

	dt.meth.chrX.per.tissue<-list() # data table to store
	#for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){ # for all tissue types
	for(my.tissue in c("PT","PT.SS.Tech","PT.SS.Bio",avail.tissues[!avail.tissues %in% "PA"])){ # for all tissue types
		## Load this RoadMap tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
		cat("Filtering out autosomes...\n")
		dt.query=dt.query[V1=="X"] #chrX only

		# "ensembl_gene_id" for all region
		dt.meth.per.region=list() # to save %met by cpg contexts
		for(i in my.cpg.contexts){
			cat(paste0("Processing ", i, "...\n"))
			gr.subject<-get.region(i)
			gr.subject<-gr.subject[seqnames(gr.subject)=="chrX",] # chrX only
			subject<-with(as.data.frame(gr.subject), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'

			# strand-aware (also "ensembl_gene_id" available by default)
			if(grepl('Gene',i) || grepl('Promo',i)){
				subject$ensembl_gene_id<-names(gr.subject) # isa 'data.frame'

				subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
				query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)
				group.keys.list="ensembl_gene_id,V1,V3,Start,End"   # same as above 
				#> head(subject)
				#  Chromosome Start   End Strand ensembl_gene_id
				#1       chr1 10369 12368      + ENSG00000223972
				#2       chr1 29307 31306      - ENSG00000227232
			# non-strand-aware (Tiling5K, CPGi, CPGi-shores, CPGi_Genes_5K, Enhancer, etc)
			}else{
				subject.keys=c("Chromosome","Start","End")
				query.key=c("V1","V2","End")
				group.keys.list="ensembl_gene_id,V1,Strand,Start,End"

				##################################################################################################
				## Overlapped region between the extended genes (+-2Kb) and the regions of interests (e.g. CPGi) #
				##################################################################################################
				gr.overlap<-mergeByOverlaps(gr.genes, gr.subject) # any overlap size>0
				#> gr.overlap
				#DataFrame with 31141 rows and 9 columns
				#gr.genes                                    symbol                                          entrezID       CpG        GC ensembl_gene_id                  gr.subject     CpG.1      GC.1
				#<GRanges>                               <character>                                       <character> <integer> <integer>     <character>                   <GRanges> <integer> <integer>
				#1         chr1:-:[ 14363,  29806]                                    WASH7P                                            653635       380      8441 ENSG00000227232     chr1:*:[ 28736,  29810]       116       787
				#2         chr1:+:[ 29554,  31109] MIR1302-10;MIR1302-11;MIR1302-2;MIR1302-9           100302278;100422831;100422834;100422919        40       760 ENSG00000243485     chr1:*:[ 28736,  29810]       116       787
				#########################################################################################
				## Regions of interests (e.g. CPGi) where they overlap with the extended genes (+-2Kb) ##
				#########################################################################################
				# 'Start' and 'End' from gr.subject (e.g. CPGi)
				subject.with.gene<-with(gr.overlap, data.frame(Chromosome=seqnames(gr.overlap$gr.subject), Start=start(gr.overlap$gr.subject), End=end(gr.overlap$gr.subject), Strand=strand(gr.overlap$gr.subject), ensembl_gene_id=ensembl_gene_id, Strand.gene=strand(gr.overlap$gr.genes)))
				#######################################################################
				## Left join with the 'original' subject with the 'annotated' subject #
				## Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
				#######################################################################
				subject<-merge(subject, subject.with.gene, all.x =TRUE) 
				#> head(subject, n=20)
				#   Chromosome  Start    End Strand ensembl_gene_id Strand.gene
				#4        chrX 166505 167721      *            <NA>        <NA>
				#5        chrX 170415 170686      * ENSG00000228572           +
				#6        chrX 171660 171867      * ENSG00000228572           +
				#7        chrX 172619 172915      * ENSG00000228572           +
				#14       chrX 221535 221828      * ENSG00000182378           +
				#15       chrX 221535 221828      * ENSG00000178605           -
				rm(subject.with.gene)
			}
			cat("\tRemoving leading 'chr'...\n")
			#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
			subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X, chr13=>13
			cat("\tconvert to data.table...\n")
			dt.subject<-as.data.table(subject) # isa data.table

			setkeyv(dt.subject,subject.keys)

			cat("\tFinding CpGs overlapping with ",i,"...\n")
			system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L)) # chrX only

			cat("\tAggregating %met...\n")
			dt.meth<-dt.overlap[,list(sum(V5.x),sum(V6.x),sum(V5.y),sum(V6.y),num.sites=.N,i,my.tissue),by=group.keys.list]
			#> dt.meth 
			#   ensembl_gene_id V1 V3     Start       End   V1   V2  V3  V4 num.sites           i my.tissue
			#1: ENSG00000251848  X  +   2728441   2730440   22   29   8  12         1 Promo_15.03        SB
			#2: ENSG00000235483  X  -   2771302   2773301   57   58  37  41         4 Promo_15.03        SB
			#3: ENSG00000229851  X  +   2821445   2823444   24   25  20  22         2 Promo_15.03        SB
			setnames(dt.meth,c("ensembl_gene_id","chr","strand","start","end","female.c","female.cov","male.c","male.cov","num.sites","region","tissue"))
			rm(subject,gr.subject,dt.overlap)

			dt.meth.per.region[[i]]=dt.meth # save this cpg context
		}# end of cpg.contexts
		dt.meth.chrX.per.tissue[[my.tissue]]<-data.table::rbindlist(dt.meth.per.region) # merged all cpg contexts of this tissue
		cat(paste0("Done for ", my.tissue, "...\n"))
	}# end of my.tissue 
	save(dt.meth.chrX.per.tissue, file=file.path("~/results/RoadMap/X-inactivation/RData",paste("dt.meth.chrX.per.tissue",my.cpg.type,gene_size,"RData",sep=".")))
}#prep.dt.meth.chrX.per.tissue

####################
## Same as above  ##
## Only for miR   ##
####################
prep.dt.meth.chrX.mir.per.tissue<-function(my.cpg.type="CG", my.cpg.contexts=c("miR_Genes", "miR_Genes_01.01","miR_Promo", "miR_Promo_01.01", "CPGi"), gene_size="miR_Genes_10.10"){
	cat(paste("C-context:",my.cpg.type,"...\n"))
	cat(paste("Gene size:",gene_size,"...\n"))
	gr.genes<-get.region(gene_size)
	mcols(gr.genes)<-DataFrame(mcols(gr.genes), mirbase_id=gr.genes$Name) # add mirbase_id as a column

	dt.meth.chrX.mir.per.tissue<-list() # data table to store
	for(my.tissue in c("PT","PT.SS.Tech","PT.SS.Bio",avail.tissues[!avail.tissues %in% "PA"])){ # for all tissue types
		## Load this RoadMap tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
		cat("Filtering out autosomes...\n")
		dt.query=dt.query[V1=="X"] #chrX only

		# "mirbase_id" for all region
		dt.meth.per.region=list() # to save %met by cpg contexts
		for(i in my.cpg.contexts){
			cat(paste0("Processing ", i, "...\n"))
			gr.subject<-get.region(i)
			gr.subject<-gr.subject[seqnames(gr.subject)=="chrX",] # chrX only
			subject<-with(as.data.frame(gr.subject), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'

			# strand-aware (also "mirbase_id" available by default)
			if(grepl('Gene',i) || grepl('Promo',i)){
				if(grepl('Gene',i)){
				subject$mirbase_id<-gr.subject$Name # isa 'data.frame'
				}else if(grepl('Promo',i)){
				subject$mirbase_id<-gr.subject$name # isa 'data.frame'
				}else{
					stop(paste0(i, " Not supported\n"))
				}
				subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
				query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)
				group.keys.list="mirbase_id,V1,V3,Start,End"   # same as above 
				#> head(subject)
				#  Chromosome     Start       End Strand    mirbase_id
				#1       chrX 151596961 151599160      - hsa-mir-105-1
				#2       chrX 151596484 151598683      - hsa-mir-105-1
				#3       chrX 151596484 151598683      - hsa-mir-105-2
			# non-strand-aware (Tiling5K, CPGi, CPGi-shores, CPGi_Genes_5K, Enhancer, etc)
			}else{
				subject.keys=c("Chromosome","Start","End")
				query.key=c("V1","V2","End")
				group.keys.list="mirbase_id,V1,Strand,Start,End"
				##################################################################################################
				## Overlapped region between the extended genes (+-2Kb) and the regions of interests (e.g. CPGi) #
				##################################################################################################
				gr.overlap<-mergeByOverlaps(gr.genes, gr.subject) # any overlap size>0
				#> gr.overlap
				#DataFrame with 31141 rows and 9 columns
				#gr.genes                                    symbol                                          entrezID       CpG        GC ensembl_gene_id                  gr.subject     CpG.1      GC.1
				#<GRanges>                               <character>                                       <character> <integer> <integer>     <character>                   <GRanges> <integer> <integer>
				#1         chr1:-:[ 14363,  29806]                                    WASH7P                                            653635       380      8441 ENSG00000227232     chr1:*:[ 28736,  29810]       116       787
				#2         chr1:+:[ 29554,  31109] MIR1302-10;MIR1302-11;MIR1302-2;MIR1302-9           100302278;100422831;100422834;100422919        40       760 ENSG00000243485     chr1:*:[ 28736,  29810]       116       787
				#########################################################################################
				## Regions of interests (e.g. CPGi) where they overlap with the extended genes (+-2Kb) ##
				#########################################################################################
				# 'Start' and 'End' from gr.subject (e.g. CPGi)
				subject.with.gene<-with(gr.overlap, data.frame(Chromosome=seqnames(gr.overlap$gr.subject), Start=start(gr.overlap$gr.subject), End=end(gr.overlap$gr.subject), Strand=strand(gr.overlap$gr.subject), mirbase_id=mirbase_id, Strand.gene=strand(gr.overlap$gr.genes)))
				#######################################################################
				## Left join with the 'original' subject with the 'annotated' subject #
				## Note that ensembl_gene_id can be null <NA> because of left-join (all.x=TRUE)
				#######################################################################
				subject<-merge(subject, subject.with.gene, all.x =TRUE) 
				#> head(subject, n=20)
				#   Chromosome  Start    End Strand ensembl_gene_id Strand.gene
				#4        chrX 166505 167721      *            <NA>        <NA>
				#5        chrX 170415 170686      * ENSG00000228572           +
				#6        chrX 171660 171867      * ENSG00000228572           +
				#7        chrX 172619 172915      * ENSG00000228572           +
				#14       chrX 221535 221828      * ENSG00000182378           +
				#15       chrX 221535 221828      * ENSG00000178605           -
				rm(subject.with.gene)
			}
			cat("\tRemoving leading 'chr'...\n")
			#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
			subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X, chr13=>13
			cat("\tconvert to data.table...\n")
			dt.subject<-as.data.table(subject) # isa data.table

			setkeyv(dt.subject,subject.keys)

			cat("\tFinding CpGs overlapping with ",i,"...\n")
			# Usually, x is a very large data.table with small interval ranges, and y is much smaller keyed data.table with relatively larger interval spans
			# If x is keyed, by.x is equal to key(x), else key(y). by.y defaults to key(y).
			system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L)) # chrX only

			cat("\tAggregating %met...\n")
			dt.meth<-dt.overlap[,list(sum(V5.x),sum(V6.x),sum(V5.y),sum(V6.y),num.sites=.N,i,my.tissue),by=group.keys.list]
			#> dt.meth 
			#   ensembl_gene_id V1 V3     Start       End   V1   V2  V3  V4 num.sites           i my.tissue
			#1: ENSG00000251848  X  +   2728441   2730440   22   29   8  12         1 Promo_15.03        SB
			#2: ENSG00000235483  X  -   2771302   2773301   57   58  37  41         4 Promo_15.03        SB
			#3: ENSG00000229851  X  +   2821445   2823444   24   25  20  22         2 Promo_15.03        SB
			setnames(dt.meth,c("mirbase_id","chr","strand","start","end","female.c","female.cov","male.c","male.cov","num.sites","region","tissue"))
			rm(subject,gr.subject,dt.overlap)

			dt.meth.per.region[[i]]=dt.meth # save this cpg context
		}# end of cpg.contexts
		dt.meth.chrX.mir.per.tissue[[my.tissue]]<-data.table::rbindlist(dt.meth.per.region) # merged all cpg contexts of this tissue
		cat(paste0("Done for ", my.tissue, "...\n"))
	}# end of my.tissue 
	save(dt.meth.chrX.mir.per.tissue, file=file.path("~/results/RoadMap/X-inactivation/RData",paste("dt.meth.chrX.mir.per.tissue",my.cpg.type,gene_size,"RData",sep=".")))
}#prep.dt.meth.chrX.per.tissue


#my.cpg.type: CHH or CHG
#see 'bin/R/Healty_Placenta/chrX.methylation.WGBS.R' for PT.CG.dt.merged (results from legacy RnBeads)
#see 'bin/R/PT.WGoxBS/local.R' for pooled data without 10x 
process.PT.bed<-function(my.cpg.type){
	cat(paste0("Processing PT:", my.cpg.type,"\n"))

	dt.meta<-fread("~/data/RoadMap/BS-Seq/MJ.CH.meta.csv", stringsAsFactor=FALSE)
	############################################
	# 1. 2 female AGA samples (A001 and A005) ##
	############################################
	cat("\tinit female1 and female2...\n")
	print(system.time(female1<-lapply(dt.meta[Context==my.cpg.type & Sex=="female" & Chr!="chrY" & Barcode=="A001",File], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)
	print(system.time(female2<-lapply(dt.meta[Context==my.cpg.type & Sex=="female" & Chr!="chrY" & Barcode=="A005",File], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) 

	names(female1)<-dt.meta[Context==my.cpg.type & Sex=="female" & Chr!="chrY" & Barcode=="A001",Chr]
	names(female2)<-dt.meta[Context==my.cpg.type & Sex=="female" & Chr!="chrY" & Barcode=="A005",Chr]
	# V1: Chr
	# V2: Start (0-based)
	# V3: End (1-based)
	# V4: %met
	# V5: Depth-of-coverage 
	# V6: Strand
	cat("\tmerging female1 and female2...\n")
	print(system.time(dt.merged.female.list<-lapply(names(female1), function(i){
											merge(
												female1[[i]][V5>=min.doc], 
												female2[[i]][V5>=min.doc], 
												by=c("V1","V2","V3","V6")
											)
										}
		) # isa 'list' of data.table (by chr)
	))
	dt.merged.female<-rbindlist(dt.merged.female.list)
	#V4.x: %met
	#V5.x: depth-of-cov
	#> dt.merged.female
	#      V1        V2        V3 V6 V4.x V5.x V4.y V5.y
	#1: chr10     60762     60763  +    0   11    0   11
	#2: chr10     62006     62007  -    0   13    0   11

	# update data.table as following
	# V1: chrX=>X
	# V2: remove 0-based start position
	# V4: %metC to count of metC
	dt.merged.female[,c("V1","V2","V4.x","V4.y"):=list(substr(V1,4,5), NULL, round(V4.x*V5.x/100,2), round(V4.y*V5.y/100,2) )]
	rm(female1,female2,dt.merged.female.list)

	##########################################
	# 2. 2 male AGA samples (A002 and A008) ##
	##########################################
	cat("\tinit male1 and male2...\n")
	print(system.time(male1<-lapply(dt.meta[Context==my.cpg.type & Sex=="male" & Chr!="chrY" & Barcode=="A002",File], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)
	print(system.time(male2<-lapply(dt.meta[Context==my.cpg.type & Sex=="male" & Chr!="chrY" & Barcode=="A008",File], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) 

	names(male1)<-dt.meta[Context==my.cpg.type & Sex=="male" & Chr!="chrY" & Barcode=="A002",Chr]
	names(male2)<-dt.meta[Context==my.cpg.type & Sex=="male" & Chr!="chrY" & Barcode=="A008",Chr]

	cat("\tmerging male1 and male2...\n")
	print(system.time(dt.merged.male.list<-lapply(names(male1), function(i){
											merge(
												male1[[i]][V5>=min.doc], 
												male2[[i]][V5>=min.doc], 
												by=c("V1","V2","V3","V6")
											)
										}
		) # isa 'list' of data.table (by chr)
	))
	dt.merged.male<-rbindlist(dt.merged.male.list)
	dt.merged.male[,c("V1","V2","V4.x","V4.y"):=list(substr(V1,4,5), NULL, round(V4.x*V5.x/100,2), round(V4.y*V5.y/100,2) )]
	rm(male1,male2,dt.merged.male.list)
	###########################
	## Merge Female and Male ##
	###########################
	cat("\tmerging female and male...\n")
	print(system.time(dt.merged<-merge(
				# 2 male merged by Chr, Position, Strand
				# V4: mean count-metC of two female
				# V5: mean depth-cov of two female
				dt.merged.female[,list(V4=round(mean(c(V4.x,V4.y)),2), V5=round(mean(c(V5.x,V5.y)),2)), by="V1,V3,V6"],
				dt.merged.male[,list(V4=round(mean(c(V4.x,V4.y)),2), V5=round(mean(c(V5.x,V5.y)),2)), by="V1,V3,V6"],
				by=c("V1","V3","V6")
			)
	))
	#V1:chr, V2:Pos, V3:Str, V5.x:cnt-met-female, V6.x:depth-cov-female, V5.y:cnt-met-male. V6.y:depth-cov-male
	setnames(dt.merged, c("V1","V2","V3","V5.x","V6.x","V5.y","V6.y"))
	assign(paste("PT",my.cpg.type,"dt.merged",sep="."), dt.merged)

	cat("\tsaving dt.merged...\n")
	my.new.RData<-file.path(paste("~/data/RoadMap/BS-Seq/Schultz2015/PT",my.cpg.type, "dt.merged.RData",sep="."))
	save(list=paste("PT",my.cpg.type,"dt.merged",sep="."), file=my.new.RData)

	rm(dt.merged)
	eval(parse(text=paste0('rm(PT.',my.cpg.type,'.dt.merged)'))) # remove the original object from memory
}# end of read.PT.bed

# used within bin/R/RoadMap/find.overlap.promoter.cpgi.R
plot.intersect<-function(gr.first, gr.second, my.names, is.ignore.strand=TRUE, is.plot=FALSE){
	# intersection betweeen the two regions
	gr.overlap<-intersect(gr.first, gr.second, ignore.strand=is.ignore.strand)
	# the size of overlapped regions, cpgi, and promoter region
	overlap.size<-sapply(split(reduce(gr.overlap, ignore.strand=is.ignore.strand), seqnames(reduce(gr.overlap, ignore.strand=is.ignore.strand))), function(i) sum(width(i)))
	first.size<-sapply(split(reduce(gr.first, ignore.strand=is.ignore.strand), seqnames(reduce(gr.first, ignore.strand=is.ignore.strand))), function(i) sum(width(i)))    # Promoter
	second.size<-sapply(split(reduce(gr.second, ignore.strand=is.ignore.strand), seqnames(reduce(gr.second, ignore.strand=is.ignore.strand))), function(i) sum(width(i))) # CPGi
	df.overlap<-data.frame(overlap.size, first.size[names(overlap.size)], second.size[names(overlap.size)]) # isa 'data.frame'
	colnames(df.overlap)=c("overlap.size",my.names)

	if(is.plot){
		dummy<-rbind(data.frame(`Chromosome`=rownames(df.overlap), `Type`=paste0("overlap/",colnames(df.overlap)[1]), ratio=df.overlap[[1]] / df.overlap[[2]]), 
					data.frame(`Chromosome`=rownames(df.overlap), `Type`=paste0("overlap/",colnames(df.overlap)[2]), ratio=df.overlap[[1]] / df.overlap[[3]]))
		print(ggplot2::ggplot(dummy, aes(Chromosome, ratio*100, fill=Type)) + geom_bar(stat = "identity", position="dodge") + labs(x="Chromosome",y="% overlap") + scale_fill_manual(values=cbPalette))
	}

	#> df.overlap
	#      overlap.size promo.size cpgi.size
	#chr1        963059   10374261   1881629
	#chr2        683829    7766303   1379397
	return(df.overlap)
}

load.my.tissue.dt.merged<-function(my.tissue,my.cpg.type){
	## Load this RoadMap tissue
	my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,my.cpg.type,"dt.merged.RData",sep="."))
	cat(paste0("Loading ", my.tissue, " ", my.cpg.type, "...\n"))
	print(system.time(load(my.new.RData))) # load list of data.table
	cat(paste(my.tissue,my.cpg.type,"dt.merged loaded\n",sep="."))
	this.tissue<-get(paste(my.tissue,my.cpg.type,"dt.merged",sep=".")) # isa list of data.table (STL2, STL3)
	eval(parse(text=paste0('rm(', my.tissue,'.',my.cpg.type,'.dt.merged)'))) # remove the original object from memory

	## Merge all chromosomes
	if(!grepl("PT",my.tissue)){
		this.tissue[["chrY"]] <- NULL # remove chrY
		this.tissue[["chrL"]] <- NULL # remove chrL
		cat("rbindlist for this.tissue...\n")
		print(system.time(dt.query<-rbindlist(this.tissue))) # collapse all chromosome which is 'listed'
	}else{
		dt.query<-this.tissue
	}
	rm(this.tissue)
	#>dt.query
	#   V1    V2 V3  V4 V7 V5.x V6.x V5.y V6.y
	#1: 10 66269  + CGG  1   12   12   44   51
	#2: 10 66273  + CGG  1   14   14   36   53
	#3: 10 71723  + CGA  1    3   10    6   16
	dt.query[,End:=V2] # add this column
	#>dt.query
	#   V1    V2 V3  V4 V7 V5.x V6.x V5.y V6.y   End
	#1: 10 66269  + CGG  1   12   12   44   51 66269
	#2: 10 66273  + CGG  1   14   14   36   53 66273
	#3: 10 71723  + CGA  1    3   10    6   16 71723

	cat("some integer columns as numeric...\n") # why?
	this.col=c("V5.x","V6.x","V5.y","V6.y")
	print(system.time(dt.query[, (this.col) := lapply(.SD,as.numeric), .SDcols=this.col]))
	return(dt.query)
}

prep.dmr<-function(my.tissue,my.cpg.type){
	# this will be made
	my.RData=file.path("~/results/RoadMap/CH/RData",paste(my.tissue,my.cpg.type,"dt.dmr.RData",sep="."))
	if(file.exists(my.RData)){
		cat(paste0(my.RData," already exists\n"))
	}else{
		## Load this RoadMap tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

		if(my.cpg.type=="CG"){
			cpg.contexts=c("Sites","Tiling500","Tiling5K","Genes","Promo_15.05","CPGi","CPGi-shores","Enhancer")
		}else{
			cpg.contexts=c("Tiling500","Tiling5K","Genes","Promo_15.05","CPGi","Enhancer")
		}
		dt.dmr.per.region=list() # to save %met by cpg contexts
		for(i in cpg.contexts){
			cat(paste0("Processing ", i, "...\n"))
			if(i=='Sites'){
				print(system.time(dt.dmr<-dt.query[,.(c.f=V5.x,
													t.f=V6.x-V5.x,
													c.m=V5.y,
													t.m=V6.y-V5.y,
													met=round((V5.x+V5.y)/(V6.x+V6.y),5), 
													met.f=round(V5.x/V6.x,5), 
													met.m=round(V5.y/V6.y,5), 
													p.value=fisher.test(matrix(c(round(V5.x),round(V6.x-V5.x),round(V5.y),round(V6.y-V5.y)),nrow=2),conf.int=F)$p.value,
													#p.value=fast.fisher(matrix(c(round(V5.x),round(V6.x-V5.x),round(V5.y),round(V6.y-V5.y)),nrow=2),conf.int=F)$p.value,
													num.sites=1,
													region=i,
													tissue=my.tissue),
													by="V1,V3,V2,V2"]
				))
			}else{
				subject<-with(as.data.frame(get.region(i)), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'

				cat("\tRemoving leading 'chr'...\n")
				#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
				subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X, chr13=>13
				cat("\tconvert to data.table...\n")
				dt.subject<-as.data.table(subject) # isa data.table
				rm(subject)

				if(grepl('^Gene',i) || grepl('^Promo',i)){
					subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
					query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)
					group.keys=c("V1","V3","Start","End")               # V3: strand
					group.keys.list="V1,V3,Start,End"                   # same as above 
				# cpgislands, Tiling, genomeTiling500
				}else{
					subject.keys=c("Chromosome","Start","End")
					query.key=c("V1","V2","End")
					group.keys=c("V1","Strand","Start","End")
					group.keys.list="V1,Strand,Start,End"
				}
				setkeyv(dt.subject,subject.keys)

				cat("\tFinding CpGs overlapping with ",i,"...\n")
				print(system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L)))

				cat("\tAggregating %met...\n")
				#print(system.time(dt.dmr<-dt.overlap[,list(f.c=sum(V5.x),f.t=sum(V6.x)-sum(V5.x),m.c=sum(V5.y),m.t=sum(V6.y)-sum(V5.y),num.sites=.N,i,my.tissue),by=group.keys.list]))
				print(system.time(dt.dmr<-dt.overlap[,list(c.f=sum(V5.x),
															t.f=sum(V6.x)-sum(V5.x),
															c.m=sum(V5.y),
															t.m=sum(V6.y)-sum(V5.y),
															met=round(sum(V5.x,V5.y)/sum(V6.x,V6.y),3), 
															met.f=round(sum(V5.x)/sum(V6.x),3), 
															met.m=round(sum(V5.y)/sum(V6.y),3), 
															p.value=fisher.test(matrix(c(sum(V5.x),sum(V6.x)-sum(V5.x),sum(V5.y),sum(V6.y)-sum(V5.y)),nrow=2))$p.value,
															#p.value=fast.fisher(matrix(c(sum(V5.x),sum(V6.x)-sum(V5.x),sum(V5.y),sum(V6.y)-sum(V5.y)),nrow=2),conf.int=F)$p.value,
															num.sites=.N,
															region=i,
															tissue=my.tissue),
													by=group.keys.list]
				))
				#     user    system   elapsed
				#	 10936.935     0.048 10936.757
				#>dt.dmr
				#          V1 Strand     Start       End  met met.f met.m   p.value num.sites         i my.tissue
				#       1:  1      *     10501     11000 0.00  0.00  0.00 0.4921348        11 Tiling500        PT
				#		2:  1      *     12501     13000 0.00  0.00  0.00 1.0000000         5 Tiling500        PT
				#		3:  1      *     13001     13500 0.01  0.00  0.01 0.1970864         6 Tiling500        PT
				rm(dt.overlap)
			}# end if i=='sites'
			setnames(dt.dmr,c("chr","strand","start","end","c.f","t.f","c.m","t.m","met","met.f","met.m","p.value","num.sites","region","tissue"))

			# avoid storing single CH sites
			# this is to avoid 'too large for hashing' error of rbindlist
			dt.dmr.per.region[[i]]=dt.dmr
		}# end of for i cpg.contexts
		save(dt.dmr.per.region, file=my.RData)
		cat("dt.dmr.per.region saved\n")
	} 
} # end of prep.dmr


update.dmr.rank<-function(my.tissue,my.cpg.type){
	# read this one
	my.RData=file.path("~/results/RoadMap/CH/RData",paste(my.tissue,my.cpg.type,"dt.dmr.RData",sep="."))
	if(!file.exists(my.RData)){
		stop(paste0(my.RData," not found\n"))
	}else{
		print(system.time(load(my.RData)))
		cat(paste0(my.RData," read\n"))
		cat("dt.dmr.per.region loaded\n")
		# for each region
		for(i in names(dt.dmr.per.region)){
			cat(paste0("Processing ", i, "...\n"))
			dt.dmr<-dt.dmr.per.region[[i]]

			# rank of 1) absolute diff, 2) relative diff, 3) p-value
			col1<-nrow(dt.dmr)+1-rank(abs(dt.dmr[,met.f-met.m]), ties.method="max") # bigger is the winnder (absolute meth diff)
			col2<-nrow(dt.dmr)+1-rank(dt.dmr[,log2((met.f+0.01)/(met.m+0.01))], ties.method="max") # relative meth diff (quotient)
			col3<-rank(dt.dmr[,p.value], ties.method="min") # smaller is the winner (high rank) (by P-value)

			# add combined rank based on above
			dt.dmr.per.region[[i]][,rank:=mapply(max,col1,col2,col3)]
			dt.dmr.per.region[[i]]<-dt.dmr.per.region[[i]][order(rank)] # order by rank
		}
		# save new with 'combinedRank' 
		my.new.RData=file.path("~/results/RoadMap/CH/RData",paste(my.tissue,my.cpg.type,"dt.dmr.rank.RData",sep="."))
		save(dt.dmr.per.region, file=my.new.RData) 
		cat(paste0(my.new.RData," saved\n"))
	}
}# end of update.dmr.rank


# used by bin/R/RoadMap/find.dmr.by.gender.R
# myMart, gr.genes should be pre-defined
find.nearest.genes<-function(gr.target){
	fields <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description") # for biomart purpose
	# Method 1: overlap with gr.genes
	if(FALSE){
		gr.overlap<-mergeByOverlaps(gr.genes, gr.target) # isa 'DataFrame' (any overlap size>0)
		dt.target.with.gene<-as.data.table(with(gr.overlap, data.frame(chr=seqnames(gr.overlap$gr.target), 
													strand=strand(gr.overlap$gr.target), 
													start=start(gr.overlap$gr.target), 
													end=end(gr.overlap$gr.target), 
													met=gr.overlap$met, 
													ensembl_gene_id=ensembl_gene_id
													)))
		dt.target.with.gene<-merge(dt.target.with.gene, 
								getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.target.with.gene[,unique(ensembl_gene_id)], mart = myMart),
								by="ensembl_gene_id")
	}

	# Method 2: nearest
	# IRanges::distanceTonearest(gr.target, gr.genes)
	dt.nearest=as.data.table(IRanges::distanceToNearest(gr.target, gr.genes)) # IRanges
	mcols(gr.target)<-DataFrame(mcols(gr.target), ensembl_gene_id=gr.genes[dt.nearest[,subjectHits]]$ensembl_gene_id, distance=dt.nearest[,distance])

	dt.target.with.gene<-merge(as.data.frame(gr.target), 
								getBM(attributes = fields, filters = "ensembl_gene_id", values = unique(gr.target$ensembl_gene_id) , mart = myMart),
								by="ensembl_gene_id",
								all.x=TRUE) # left-join
	return(as.data.table(dt.target.with.gene))
}

#used by 'bin/R/RoadMap/plot.met.exp.by.gender.R'
#http://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame
is.nan.data.frame <- function(x)do.call(cbind, lapply(x, is.nan))
is.inf.data.frame <- function(x)do.call(cbind, lapply(x, is.infinite))

# used by bin/R/RoadMap/plot.met.exp.by.gender.R
# dt.meth.exp.fg must be pre-defined
PlotMethExp<-function(my.padj,my.target=c("CPGi","Genes")){
	# one-to-many: (1-ensg-to-2-hgnc)
	fields <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "gene_biotype", "description")
	if(my.padj<=0.5){
		# get hgnc_symbol (gene name)
		my.dt.meth.exp<-merge(dt.meth.exp.fg[new.padj<=my.padj], getBM(attributes = fields, filters = "ensembl_gene_id", values = dt.meth.exp.fg[new.padj<my.padj,unique(ensembl_gene_id)], mart = myMart), by="ensembl_gene_id")
	}else{
		my.dt.meth.exp=dt.meth.exp.fg
	}

	if(my.padj==1){
		p<-ggplot(my.dt.meth.exp[!is.na(eval(parse(text=my.target[1]))) & !is.na(eval(parse(text=my.target[2])))], aes_string(x=my.target[1], y=my.target[2]))
	}else{
		p<-ggplot(my.dt.meth.exp[new.padj<=my.padj & !is.na(eval(parse(text=my.target[1]))) & !is.na(eval(parse(text=my.target[2])))], aes_string(x=my.target[1], y=my.target[2]))
	}
	p<-p +
		geom_point(aes(col=new.padj,shape=`Expression?`,size=2^abs(log2FoldChange)),alpha=.9) + 
		ggtitle(paste("% Diff Meth (chrX & padj<=",my.padj,")")) +
		labs(my.target) +
		geom_hline(yintercept=0) + geom_vline(xintercept=0)

	if(my.padj<=0.5){
		p<-p+geom_text(aes(label=my.dt.meth.exp[!is.na(eval(parse(text=my.target[1]))) & !is.na(eval(parse(text=my.target[2]))),hgnc_symbol]), hjust =1.2, nudge_x=0.05, size=4)
	}
	print(p+theme_Publication())
}

