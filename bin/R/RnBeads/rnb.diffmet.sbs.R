#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

library(RnBeads)

logger.start("Logger started", fname = NA )
##########
# Config #
##########
cat("Loading configuation...\n")
source("~/Pipelines/config/RnBeads.sbs.R")
#source("~/Pipelines/lib/methylkit.R") # defines fast.fisher (for 5hmC only)
cat("Loading configuration done...\n")

hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

# pdf output filename 
#pdf.file.name<- file.path(run.dir,paste0('rnb.met.diff.plot.',time.stamp,'.pdf'))
#pdf(file=pdf.file.name)

#####################################
## 1. Load RnBeads Diff Met Result ##
#####################################
# ~/results/RnBeads/5hmC.SGA.AGA/CpG/RData/new.rnb.meth.diff.all.RData
# ~/results/RnBeads/5hmC.SGA.AGA/CpG/RData/rnb.meth.diff.all.RData
if(my.project=="5hmC.SGA.AGA"){
	rnb.dm.rdata<-file.path(rdata.dir,paste0("new.rnb.meth.diff.",my.region,".RData")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
}else{
	#rnb.dm.rdata<-file.path(rdata.dir,paste("rnb.meth.diff",my.region,my.project,"RData",sep=".")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
	rnb.dm.rdata<-file.path(rdata.dir,paste("rnb.meth.diff",my.region,my.project,my.run.date,"RData",sep=".")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
}

if(file.exists(rnb.dm.rdata)){
	cat(paste0("Loading rnb.dm.rdata: ",rnb.dm.rdata,"...\n"))
	load(rnb.dm.rdata) # loading "rnb.diff.met"
	cat("rnb.diff.met loaded\n")
}else{ rnb.diff.met<-list() # this will be set at the end 
	# load RnBead Objects (preprocessed.rnb.set.dir defined in 'config/RnBeads.sbs.R'
	if(file.exists(preprocessed.rnb.set.dir)){ 
		# 1. load RnBSet object
		cat(paste0("Importing rnb.set (preprocessed RnBeads data) from ",preprocessed.rnb.set.dir,"...\n"))
		# rnb.region.types.for.analysis(rnb.set) returns the region(s) that you intend to analyse
		#rnb.set<-load.rnb.set(path=preprocessed.rnb.set.dir)  # isa 'RnBSet'
		rnb.set<-load.rnb.set(path=preprocessed.rnb.set.dir) #isa 'RnBSet'
	}else{
		stop("preprocessed.rnb.set.dir not foundi. Stopped!")
	}

	# my.run.date defined in config/RnBeads.sbs.R
	# ~/results/RnBeads/SGA.AGA/CpG/all/reports-2015-01-21_04_30_PM/differential_rnbDiffMeth
	rnb.dm.input<-file.path(run.dir, paste0(my.run.date,"/differential_rnbDiffMeth")) 
	if(file.exists(rnb.dm.input)){ 
		###########################
		# User-defined annotation #
		###########################
		# function 'set.user.annotation' defined
		if(!grepl("FLD",my.project)){
			source("~/Pipelines/bin/R/RnBeads/set.user.annotation.sbs.R") 
		}

		# 2. load annotation data
		rnb.anno<-list()
		for(i in c("sites",rnb.region.types.for.analysis(rnb.set))){
			cat(paste0("Loading rnb.anno (annotation) for ",i,"...\n"))
			rnb.anno[[i]] <-annotation(rnb.set, type=i, add.names=TRUE) # isa list of data.frame

		}
		# 3. load methylation data
		if(grepl("FLD",my.project)){
			rnb.meth<-list()
			for(i in c("sites",rnb.region.types.for.analysis(rnb.set))){
				cat(paste0("Loading rnb.meth (methylation data) for ",i,"...\n"))
				rnb.meth[[i]] <-meth(rnb.set, type=i, row.names=TRUE)       # methylation level of my.region
			}
		}
		# 4. load diff meth data
		if(file.exists(file.path(rnb.dm.input,"rnbDiffMeth_tables.ffData"))){
			cat(paste0("Loading ",rnb.dm.input, " (diff methylation RnBeads object) for ",my.region,"...\n"))
			rnb.diff.obj<-load.rnb.diffmeth(path=rnb.dm.input) # isa 'RnBDiffMeth'

			#cat("rnb.diff.obj:\n")
			#str(rnb.diff.obj)

			# 5. init rnb.dm by regional types
			rnb.dm<-list()
			for(i in c("sites",get.region.types(rnb.diff.obj))){
				cat(paste0("Loading rnb.diff.obj (diff methylation data) for ",i,"...\n"))
				rnb.dm[[i]]<-get.table(rnb.diff.obj, get.comparisons(rnb.diff.obj), i, return.data.frame=TRUE) # isa list data.frame (or matrix if return.data.frame=FALSE)
				#summary(rnb.dm[[i]])
			}
		}

		# check length(get.region.types(rnb.diff.obj)) length(rnb.region.types.for.analysis(rnb.set)) same?
		# for RnBeads, region.type="sites@
		# the coordinate is *NOT* a 'BED' format, but 1-based 'CG' sites as shown below:
		# start: 1-based 'C' location
		# end: 1-based 'G' location
		# 'rnb.RnBSet.to.GRangesList(rnb.set)' also contains 1-based 'CG' position for sites
		# 6. init rnb.diff.met 
		for(i in c("sites",rnb.region.types.for.analysis(rnb.set))){
			cat(paste0("Joining diff methylation data & annotation for ",i,"...\n"))
			if(exists("rnb.meth")){ # from FLD data
				rnb.diff.met[[i]]<-data.frame(rnb.anno[[i]], rnb.meth[[i]])
				# rnb.diff.met stores BED format coordinate as shown below
				# this is to merge with top10 CpG (my.top.cpg.bed)
				# done for Fludigm and SGA.AGA
				# but not for Boy.Girl yet!!
				if(i=="sites"){  # only 'sites' for FLD data
					rnb.diff.met[[i]]$Start=ifelse(rnb.diff.met[[i]]$Strand=="-", rnb.diff.met[[i]]$Start, rnb.diff.met[[i]]$Start-1)
					rnb.diff.met[[i]]$End=ifelse(rnb.diff.met[[i]]$Strand=="-", rnb.diff.met[[i]]$End, rnb.diff.met[[i]]$End-1)
				}
			}else{
				if(exists("rnb.dm")){
					rnb.diff.met[[i]]<-data.frame(rnb.anno[[i]], rnb.dm[[i]])
					cat(paste0("rnb.diff.met ordered by combinedRank and fdr for ",i,"...\n"))
					if(i=="sites"){
						rnb.diff.met[[i]]<-rnb.diff.met[[i]][order(rnb.diff.met[[i]]$combinedRank, rnb.diff.met[[i]]$diffmeth.p.adj.fdr),]
					}else{
						rnb.diff.met[[i]]<-rnb.diff.met[[i]][order(rnb.diff.met[[i]]$combinedRank, rnb.diff.met[[i]]$comb.p.adj.fdr),]
					}
				}
			}
			cat(paste0("rnb.diff.met created for ",i,"\n\n"))
		}

		cat(paste0("saving rnb.diff.met as a RData object to ",rnb.dm.rdata,"...\n"))
		save(rnb.diff.met, file=rnb.dm.rdata)
		cat(paste0(rnb.dm.rdata," created\n"))
	# load .csv.gz file directly
	}else{
		# Below tries to read methylation result (.csv.gz) directly
		#~/results/RnBeads/SGA.AGA.Girl/CpG/all/reports-2015-05-20_11_29_AM/differential_methylation_data/diffMethTable_site_cmp1.csv.gz
		#~/results/RnBeads/SGA.AGA.Girl/CpG/all/reports-2015-05-20_11_29_AM/differential_methylation_data/diffMethTable_region_cmp1_genomeTiling500.csv.gz
		for(i in c("sites",rnb.region.types.for.analysis(rnb.set))){
			cat(paste0("Reading csv.gz of ",i,"...\n"))
			if(i=="sites"){
				rnb.dm.input<-file.path(run.dir, paste0(my.run.date,'/differential_methylation_data/diffMethTable_site_cmp1.csv.gz')) 
			}else{
				rnb.dm.input<-file.path(run.dir, paste0(my.run.date,'/differential_methylation_data/diffMethTable_region_cmp1_',i,'.csv.gz')) 
			}
			rnb.diff.met[[i]]<-read.csv(file=rnb.dm.input, header=TRUE, stringsAsFactors=TRUE ) # a dataframe
			if(i=="sites"){
				rnb.diff.met[[i]]<-rnb.diff.met[[i]][order(rnb.diff.met[[i]]$combinedRank, rnb.diff.met[[i]]$diffmeth.p.adj.fdr),]
			}else{
				rnb.diff.met[[i]]<-rnb.diff.met[[i]][order(rnb.diff.met[[i]]$combinedRank, rnb.diff.met[[i]]$comb.p.adj.fdr),]
			}
		}
		cat(paste0("saving rnb.diff.met as a RData object to ",rnb.dm.rdata,"...\n"))
		save(rnb.diff.met, file=rnb.dm.rdata)
		cat(paste0(rnb.dm.rdata," created\n"))
	}# end of if(file.exists(rnb.dm.input)){ 
} # end of if(file.exists(rnb.dm.rdata))


###################################################
# if no p.value returned from RnBeads             #
# which is the case of one sample within a group  #
###################################################
if(my.project=="5hmC.SGA.AGA"){
	if(FALSE){
		#################
		## Get P.value ##
		## data.frame  ##
		#################
		query<-rnb.diff.met[["sites"]] # isa data.frame
		# x is a data.frame from the subject with only one element
		get.pvalue<-function(x,i){ 
			if(i=='genes' || i=='promoters' || i=='enstTiling500'){
				select<-query$Chromosome==x$Chromosome & query$Strand==x$Strand & query$Start>=x$Start & query$End<=x$End
			}else{#tiling, genomeTiling500
				select<-query$Chromosome==x$Chromosome & query$Start>=x$Start & query$End<=x$End
			}
			dummy<-query[select,]
			if(nrow(dummy)>=1){
				my.pvalue<-wilcox.test(dummy$mean.g1, dummy$mean.g2, paired=TRUE)$p.value
				#based on weighted mean 
				#my.pvalue<-wilcox.test(dummy$mean.g1*dummy$mean.covg.g1, dummy$mean.g2*dummy$mean.covg.g2, paired=TRUE)$p.value
				return(ifelse(is.nan(my.pvalue),1,my.pvalue))
			}else{
				return(1)
			}
		} #end of function get.pvalue
		#################
		## Get P.value ##
		## data.table  ##
		#################
		library(data.table)
		#dt.diff.met<-lapply(rnb.diff.met, as.data.table) # isa a list of data.table
		query<-as.data.table(rnb.diff.met[["sites"]]) # isa data.table
		# x is a data.frame from the subject with only one element
		get.dt.pvalue<-function(x,i){ 
			if(i=='genes' || i=='promoters' || i=='enstTiling500'){
				by.this=c("Chromosome","Strand","Start","End")
			}else{#tiling, genomeTiling500
				by.this=c("Chromosome","Start","End")
			}
			my.subject<-as.data.table(x) # 
			setkeyv(my.subject,by.this)
			dt.overlap<-foverlaps(query, my.subject, by.x=by.this, by.y=by.this, type="within", nomatch=0L)
			#dt.overlap<-foverlaps(query, my.subject, by.x=by.this, by.y=by.this, type="within", nomatch=0L, which=TRUE)
			#dt.overlap[,.(.N), by="yid"] # count query by subject
			#query[dt.overlap[yid==1,xid]] # get query of yid=1
			#dt.overlap[,get.dt.pvalue(.SD, yid), by="yid"]
			#dt.overlap[,lapply(.SD, get.dt.pvalue,yid), by="yid"]
			if(nrow(dt.overlap)>=1){
				my.pvalue<-wilcox.test(dt.overlap[,mean.g1], dt.overlap[,mean.g2], paired=TRUE)$p.value
				#based on weighted mean 
				#my.pvalue<-wilcox.test(dummy$mean.g1*dummy$mean.covg.g1, dummy$mean.g2*dummy$mean.covg.g2, paired=TRUE)$p.value
				return(ifelse(is.nan(my.pvalue),1,my.pvalue))
			}else{
				return(1)
			}
		} #end of function get.dt.pvalue
		#################
		## Get P.value ##
		## sqlite      ##
		#################
		library("RSQLite")
		sqlite.db="~/Pipelines/5hmC_new.db"
		#dbd <- dbDriver("SQLite")
		dbi <- dbConnect(RSQLite::SQLite(), sqlite.db)
		get.sqlite.pvalue<-function(x){ 
			chr<-as.character(x$Chromosome)
			start<-x$Start
			end<-x$End

			sql<-paste0("SELECT mean_g1, mean_g2 FROM \"5hmcCpGs\" WHERE chr='",chr,"' AND start BETWEEN ",start," AND ",end)
			overlap<-dbGetQuery(dbi,sql)
			if(nrow(overlap)>=1){
				my.pvalue<-wilcox.test(overlap$mean_g1, overlap$mean_g2, paired=TRUE)$p.value
				return(ifelse(is.nan(my.pvalue),1,my.pvalue))
			}else{
				return(1)
			}
		} #end of function get.sqlite.pvalue
	}# end of if FALSE
	#################
	## Get P.value ##
	## MySQL       ##
	#################
	#write.table(cbind(0, rnb.diff.met[["sites"]][,c("Chromosome","Start","End","Strand","mean.g1","mean.g2")]), file="~/results/RnBeads/5hmC.SGA.AGA/CpG/5hmcCpGs.txt", sep="\t", row.names=F, col.names=F, quote=F)
	library("RMySQL")
	#dbd <- dbDriver("MySQL", fetch.default.rec=10000) # create a MySQL instance and create one connection.
	#dbi <- dbConnect(dbd, dbname='PLACENTOME') # configuration file \file{\$HOME/.my.cnf}
	dbi <- dbConnect(RMySQL::MySQL(), dbname='PLACENTOME')
	get.mysql.pvalue<-function(x){ 
		chr<-as.character(x$Chromosome)
		start<-x$Start
		end<-x$End

		sql<-paste0("SELECT mean_g1, mean_g2 FROM 5hmcCpGs WHERE chr='",chr,"' AND start BETWEEN ",start," AND ",end)
		overlap<-dbGetQuery(dbi,sql)
		if(nrow(overlap)>=1){
			my.pvalue<-wilcox.test(overlap$mean_g1, overlap$mean_g2, paired=TRUE)$p.value
			return(ifelse(is.nan(my.pvalue),1,my.pvalue))
		}else{
			return(1)
		}
	} #end of function get.mysql.pvalue
	# for each regional analysis
	#my.target.region<-c("sites","tiling","genes","promoters","cpgislands","enstTiling500","genomeTiling500")
	my.target.region<-c("tiling","genes","promoters","cpgislands","genomeTiling500")
	for(i in my.target.region){
		cat("region type=",i,"\n")
		cat("calculating pvalue...\n")
		if(i=="sites"){
			query<-with(rnb.diff.met[[i]],data.frame(
				g1.numCs=round(mean.g1*mean.covg.g1),
				g1.numTs=round(mean.covg.g1 - mean.g1*mean.covg.g1),
				g2.numCs=round(mean.g2*mean.covg.g2),
				g2.numTs=round(mean.covg.g2 - mean.g2*mean.covg.g2)))
			# with single-core
			#my.pvalue.list<-apply(query, 1, function(k){round(fisher.test(matrix(k, nrow=2))$p.value,5)})

			# with multi-core
        	f.fisher=function(x){ fast.fisher(matrix( x ,ncol=2,byrow=T),conf.int = F)$p.value }
        	my.pvalue.list<-simplify2array( parallel::mclapply(split(query,1:nrow(query)), f.fisher, mc.cores=num.cores) ) # apply fisher test
			
			# P-value
			# NKX1-2: Chromosome 10: 126,135,592-126,138,753 reverse stran
			#select<-query$Chromosome=="chr10" & query$Start>=126135592 & query$End<=126140753 # TSS2000 + NKX1-2 
			#hist(log(query[select,]$combinedRank), freq=FALSE, main=paste0("NKX1-2 (n=",table(select)["TRUE"],")"), xlab="log(combinedRank)")
			#hist(log(query[!select,]$combinedRank), freq=FALSE, main=paste0("Elsewhere (other than NKX1-2, n=",table(select)["FALSE"],")"), xlab="log(combinedRank)")
			#wilcox.test(log(query[select,]$combinedRank), log(query[!select,]$combinedRank))$p.value
		}else{
			subject<-rnb.diff.met[[i]] # isa data.frame
			subject.list<-split(subject, 1:nrow(subject))
			# NKX1-2: Chromosome 10: 126,135,592-126,138,753 reverse stran
			#select<-subject$Chromosome=="chr10" & subject$Start>=126135592 & subject$End<=126140753 # TSS2000 + NKX1-2 
			#subject<-dt.diff.met[[i]][select,]
			if(FALSE){
				####################
				## data.frame way ##
				####################
				#my.pvalue.list<-simplify2array(lapply(split(subject, 1:nrow(subject)), get.pvalue, i)) # isa array (single-core)
				my.pvalue.list<-simplify2array(parallel::mclapply(subject.list, get.pvalue, i, mc.cores=num.cores)) # isa array (multi-core)
				####################
				## data.table way ##
				####################
				my.pvalue.list<-simplify2array(parallel::mclapply(subject.list, get.dt.pvalue, i, mc.cores=num.cores)) # isa array 
				################
				## sqlite-way ##
				################
				my.pvalue.list<-simplify2array(parallel::mclapply(subject.list, get.sqlite.pvalue, mc.cores=num.cores)) # isa array 
				#my.pvalue.list<-sapply(subject.list, get.sqlite.pvalue)
				dbDisconnect(dbi) # Disconnect from the database
			}#end of if FALSE
			###############
			## MySQL-way ##
			###############
			#my.pvalue.list<-simplify2array(parallel::mclapply(subject.list, get.mysql.pvalue, mc.cores=num.cores)) # isa array 
			my.pvalue.list<-sapply(subject.list, get.mysql.pvalue)
			#dbDisconnect(dbi) # Disconnect from the database
		}
		cat("updating adjusted pvalue...\n")
		# update pvalue of the rnb.diff.met
		rnb.diff.met[[i]]$comb.p.val<-my.pvalue.list
		rnb.diff.met[[i]]$comb.p.adj.fdr<-p.adjust(rnb.diff.met[[i]]$comb.p.val, method="BH") # update adjusted pvalue of the rnb.diff.met
		cat("updating combined rank...\n")
		# update the rank
		if(i=="sites"){
			col1<-nrow(rnb.diff.met[[i]])+1-rank(abs(rnb.diff.met[[i]]$mean.diff), ties.method="max") # rank.abs (absolute meth diff)
			col2<-nrow(rnb.diff.met[[i]])+1-rank(abs(rnb.diff.met[[i]]$mean.quot.log2), ties.method="max") # rank.rel (relative meth diff (quotient))
			col3<-rank(rnb.diff.met[[i]]$diffmeth.p.val, ties.method="min") # smaller is the winner (high rank) (by P-value)
		}else{
			col1<-nrow(rnb.diff.met[[i]])+1-rank(abs(rnb.diff.met[[i]]$mean.mean.diff), ties.method="max") # bigger is the winnder (absolute meth diff)
			col2<-nrow(rnb.diff.met[[i]])+1-rank(abs(rnb.diff.met[[i]]$mean.mean.quot.log2), ties.method="max") # relative meth diff (quotient)
			col3<-rank(rnb.diff.met[[i]]$comb.p.val, ties.method="min") # smaller is the winner (high rank) (by P-value)
		}
		rnb.diff.met[[i]]$combinedRank=mapply(max,col1,col2,col3) # update rank (single-core)
		rnb.diff.met[[i]]$combinedRank2=mapply(max,col1,col3) # update rank (single-core)
		rnb.diff.met[[i]]$combinedRank3=mapply(max,col2,col3) # update rank (single-core)
		#rnb.diff.met[[i]]$combinedRank=parallel::mcmapply(max,col1,col2,col3, mc.cores=num.cores) # update rank (multi-core)
		# re-order by the updated combinedRank
		rnb.diff.met[[i]]<-rnb.diff.met[[i]][order(rnb.diff.met[[i]]$combinedRank),]
	} #end of foreach region
	# save updated rnb.diff.met
	new.rnb.dm.rdata<-file.path(rdata.dir,paste0("new.rnb.meth.diff.",my.region,".RData")) # this will store 'rnb.diff.met' which is a list of data.frame for sites, promoters...
	cat(paste0("saving new.rnb.diff.met as a RData object to ",new.rnb.dm.rdata,"...\n"))
	save(rnb.diff.met, file=new.rnb.dm.rdata)
	cat(paste0(new.rnb.dm.rdata," created\n"))
} #end of if(my.project=="5hmC.SGA.AGA")

############################
## Write Top Diff by Rank ##
##### #######################
# FLD assay does not contain 'combinedRank'
# Following for AGA.Boy.Girl (WGoxBS, n=4)
## tiling5K: median(my.diff.met$num.sites)=17, length(my.diff.met$num.sites)=555791 (after filter:282567)
## genes: median(my.diff.met$num.sites)=61, length:47550 
## promoters: median:18, length:54094
## cpgislands: 51, 26620
## tiling500: 2, 3502611
if(FALSE){
	if(length(grep("combinedRank",colnames(rnb.diff.met$sites)))){
		for(i in names(rnb.diff.met)){
			cat(paste0("writing top100 for ",i,"...\n"))
			my.diff.met<-rnb.diff.met[[i]][order(rnb.diff.met[[i]]$combinedRank),]
			###########################
			# filtering based on Rank
			###########################
			if(i=="sites"){select<-TRUE}else{select<-my.diff.met$num.sites >= median(my.diff.met$num.sites)}
			top.diff.met<-head(my.diff.met[select,],n=100)
			if(my.project=="5hmC.SGA.AGA"){
				write.csv(top.diff.met, file=file.path(run.dir,paste0(my.project,".",i,".top100.by.new.rank.csv")))
			}else{
				write.csv(top.diff.met, file=file.path(run.dir,paste0(my.project,".",i,".top100.by.rank.csv")))
				#write.csv(my.diff.met, file=gzfile(file.path(run.dir,paste0(my.project,".",i,".all.by.rank.csv.gz"))))
			}
		}
	}
}

	# DEBUG mode
	if(FALSE){
		subject<-rnb.diff.met[[i]]
		select<-subject$num.sites>=3
		subject<-head(subject[select,],n=30)
		# get pvalue of diff.met from the function 'get.pvalue' defined and shown above
		my.pvalue.list<-lapply(split(subject, rownames(subject)), get.pvalue, i) # isa list
		# update pvalue of the rnb.diff.met
		subject$comb.p.val<-as.numeric(my.pvalue.list[rownames(subject)]) # ia a numeric vector
		# update adjusted pvalue of the rnb.diff.met
		subject$comb.p.adj.fdr<-p.adjust(subject$comb.p.val, method="BH")
		# update rank
		subject$combinedRank=mapply(
			max,
			nrow(subject)+1-rank(abs(subject$mean.mean.diff), ties.method="max"),
			nrow(subject)+1-rank(abs(subject$mean.mean.quot.log2), ties.method="max"),
			rank(subject$comb.p.val, ties.method="max")
		)
	}#end of DEBUG
	if(FALSE){
		###################
		# list of GRanges #
		###################
		# Note that 'rnb.diff.met' is based on .bed format which is 0-based
		# whereas GRange is supposed to have 1-based 
		library(GenomicRanges)
		if(TRUE){
			gr.diff.met=lapply(rnb.diff.met, makeGRangesFromDataFrame,keep.extra.columns=TRUE)
		}else{
			gr.diff.met<-lapply(rnb.diff.met, function(dummy) 
				GRanges(
					seqnames<-as.character(dummy$Chromosome), 
					ranges<-IRanges(start=dummy$Start, end=dummy$End), 
					strand<-as.character(dummy$Strand), 
					data.frame(dummy[,5:ncol(dummy)])
				)
			)
		}
		#################
		## Get P.value ##
		## GRanges     ##
		#################
		for(i in names(gr.diff.met)[2:length(names(gr.diff.met))]){
			subject<-gr.diff.met[[i]] # isa GRanges
			subject<-split(subject,names(my.gr.region)) # isa GRangesList

			lapply(subject, get.gr.pvalue)
		}

		query<-gr.diff.met[["sites"]] # isa GRanges
		get.gr.pvalue<-function(x){
			dummy<-subsetByOverlaps(query, x) # isa GRanges
			#get p-value
			#based on mean 
			#t.test(dummy$mean.g1, dummy$mean.g2, paired=TRUE)$p.value
			#wilcox.test(dummy$mean.g1, dummy$mean.g2, paired=TRUE)$p.value

			#based on weighted mean 
			#t.test(dummy$mean.g1*dummy$mean.covg.g1, dummy$mean.g2*dummy$mean.covg.g2, paired=TRUE)$p.value
			#return(wilcox.test(dummy$mean.g1*dummy$mean.covg.g1, dummy$mean.g2*dummy$mean.covg.g2, paired=TRUE)$p.value)
			return(ifelse(is.nan(wilcox.test(dummy$mean.g1*dummy$mean.covg.g1, dummy$mean.g2*dummy$mean.covg.g2, paired=TRUE)$p.value), 1, wilcox.test(dummy$mean.g1*dummy$mean.covg.g1, dummy$mean.g2*dummy$mean.covg.g2, paired=TRUE)$p.value))
		}
	} #end of if FALSE

##############################
# manhattan plot using qqman #
##############################
#my.target.region<-c("sites", "tiling","genes","promoters","cpgislands","genomeTiling500")
	library(qqman)

	my.diff.met<-rnb.diff.met[["sites"]]; select<-!is.na(my.diff.met$combinedRank)
	upper<-log(max(my.diff.met[select,]$combinedRank))
	top.diff.met<-paste(my.diff.met[select,]$Chromosome, my.diff.met[select,]$Start, my.diff.met[select,]$Strand, sep=".")[1:10]
	p.site<-with(my.diff.met[select,], 
			data.frame(
				SNP=paste(my.diff.met[select,]$Chromosome, my.diff.met[select,]$Start, my.diff.met[select,]$Strand, sep="."),
				CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
					function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
				BP=Start,
				RANK=upper-log(combinedRank),
				NUM.SITE=1
			))
	manhattan(p.site, p="RANK", logp=FALSE, ylab="", genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X")) 


#for(i in my.target.region)
if(FALSE){
	library(qqman)
	for(i in names(rnb.diff.met)){
		cat("type=",i,"\n")
		my.diff.met<-rnb.diff.met[[i]]
		###########################
		# filtering based on Rank
		###########################
		# combinedRank coule be 'NA' for Boy.Girl
		if(i=="sites"){select<-!is.na(my.diff.met$combinedRank)}else{select<-!is.na(my.diff.met$combinedRank) & my.diff.met$num.sites >= median(my.diff.met$num.sites)}
		upper<-log(max(my.diff.met[select,]$combinedRank))

		cat("making a data.frame for manhattan plot\n")
		if(i=="sites"){
			main.title = paste0("Manhattan plot of ",nrow(my.diff.met[select,])," CpG ",i)
			top.diff.met<-paste(my.diff.met[select,]$Chromosome, my.diff.met[select,]$Start, my.diff.met[select,]$Strand, sep=".")[1:10]
			dummy<-with(my.diff.met[select,], 
					data.frame(
						SNP=paste(my.diff.met[select,]$Chromosome, my.diff.met[select,]$Start, my.diff.met[select,]$Strand, sep="."),
						CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
							function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
						BP=Start,
						P=diffmeth.p.val,
						P.ADJ=diffmeth.p.adj.fdr,
						RANK=upper-log(combinedRank),
						NUM.SITE=1
					))
		}else if(i=="enstTiling500" || i=="genomeTiling500"){
			main.title = paste0("Manhattan plot of ",i," (",nrow(my.diff.met[select,])," regions of >=",median(my.diff.met$num.sites),"CpGs)")
			top.diff.met<-as.character(head(my.diff.met[select,c("id")],n=10))
			dummy<-with(my.diff.met[select,], 
					data.frame(
						SNP=as.character(id),
						CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
							function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
						BP=Start,
						P=comb.p.val,
						P.ADJ=comb.p.adj.fdr,
						RANK=upper-log(combinedRank),
						NUM.SITE=num.sites
					))
		}else{
			main.title = paste0("Manhattan plot of ",i," (",nrow(my.diff.met[select,])," regions of >=",median(my.diff.met$num.sites),"CpGs)")
			top.diff.met<-paste0(my.diff.met[select,]$Chromosome, ".", my.diff.met[select,]$Start)[1:10]
			dummy<-with(my.diff.met[select,], 
					data.frame(
						SNP=paste0(my.diff.met[select,]$Chromosome, ".", my.diff.met[select,]$Start),
						CHR=sapply(substr(as.character(Chromosome),4,nchar(as.character(Chromosome))), # chromosome (e.g. chr1, chr2, chrX, chrM)
							function(i) if(i=="X"){i<-23}else if(i=="Y"){i<-24}else if(i=="M"){i<-25}else{i<-as.numeric(i)}),
						BP=Start,
						P=comb.p.val,
						P.ADJ=comb.p.adj.fdr,
						RANK=upper-log(combinedRank),
						NUM.SITE=num.sites
					))
		} # end of if i 
		cat("plotting a manhattan plot\n")
		if(my.project=="5hmC.SGA.AGA"){
			my.filename=file.path(run.dir,paste0(my.project,'.',i,".manhattan.new.rank"))
		}else{
			my.filename=file.path(run.dir,paste0(my.project,'.',i,".manhattan"))
		}
		tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
		if(grepl("Boy.Girl",my.project)){
			manhattan(dummy, p="RANK", logp=FALSE, ylab=paste0(round(upper,2),"-log(combinedRank)"), genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X")) 
			#manhattan(dummy, genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X"), cex=0.6) 
		}else{
			manhattan(dummy, p="RANK", logp=FALSE, ylab=paste0(round(upper,2),"-log(combinedRank)"), genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X","Y")) 
			#manhattan(dummy, genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X","Y"), cex=0.6) 
		}
		dev.off()
	} # end of Manhattan plot 
} # end of FALSE

#########################################################
## Plot Methylation Difference in the Gene of Interest 
## via rnb.plot.locus.profile() of RnBeads
#########################################################
if(!exists("rnb.set")){
	# 1. load RnBSet object
	cat(paste0("Importing rnb.set (preprocessed RnBeads data) from ",preprocessed.rnb.set.dir,"...\n"))
	# rnb.region.types.for.analysis(rnb.set) returns the region(s) that you intend to analyse
	rnb.set<-load.rnb.set(path=preprocessed.rnb.set.dir) #isa 'RnBSet'
}
library(GenomicFeatures)
#up and down flank size
my.upflank=2000 #up from TSS
my.downflank=200 # down form TES

#https://support.bioconductor.org/p/62347/
fields <- c("ensembl_gene_id", "hgnc_symbol","chromosome_name","strand","start_position","end_position") # get the gene names

#my.gene.list<-c("NKX1-2","NBEA","MAB21L1","GRM7","BHLHB9","ABCB10P1","CD14")
my.gene.list<-c("C1orf216")
for(my.gene in my.gene.list){
	my.ensg=getBM(attributes = fields, filters = "hgnc_symbol", values=my.gene, mart = myMart) # is a data.frame
	my.gr.ensg<-gr.ensg[mcols(gr.ensg)$gene_id %in% my.ensg$ensembl_gene_id] # NKX1-2
	my.gr.ext.ensg<-promoters(my.gr.ensg, upstream=my.upflank, downstream=end(my.gr.ensg)-start(my.gr.ensg)+1+my.downflank) # includes 2000 +TSS and 100 +TES 

	chrom <-as.character(seqnames(my.gr.ext.ensg))
	start <-start(my.gr.ext.ensg)
	end <- end(my.gr.ext.ensg) 
	for(my.group in my.group.list){ # my.group.list defined within the config
		if(my.project=='5hmC.SGA.AGA'){
			# rnb.sample.groups returns NULL
			# as there's only one sample within a group
			# so, manually set up the group definition to display
			sample.grouping<-list(`AGA`=c(2), `SGA`=c(1))
		#}else if(my.project=='FLD.SGA.PAPPA'){
		#	sample.grouping<-list(`Normal.P`=c(2,4,6,8,10,12,14,16,18), `Low.P`=c(1,3,5,7,9,11,13,15,17))
		}else{	
			sample.grouping<-rnb.sample.groups(rnb.set)[[my.group]]
		}
		cat(paste0("making a methylation locus profile for ",my.gene,"\n"))

		my.filename<-file.path(run.dir, my.run.date, paste0(my.project,'.',my.gene, ".methylation.profile.by.",my.group))
		tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
		rnb.plot.locus.profile(rnb.set, chrom, start, end, grps=sample.grouping)
		dev.off()
	}
}

###############################
## Top 10 NKX1-2 CpG Heatmap ##
###############################
if(FALSE){
	#my.top.cpg.bed<-"~/Pipelines/data/NKX1-2/MJ.WGBS.SGA.AGA.top.chr10.cpg.bed" # top ranked CpG from MJ WGBS
	my.top.cpg.bed<-"~/Pipelines/data/NKX1-2/MJ.WGBS.SGA.AGA.top.chr10.cpg.named.bed" # top ranked CpG from MJ WGBS
	wgbs.top.cpg<-read.delim(my.top.cpg.bed, header=F)
	#colnames(wgbs.top.cpg)<-c('chr','start','end','aga','sga','strand','rank')
	#gr.top.cpg<- with(wgbs.top.cpg, GRanges(chr, IRanges(start, end), strand, rank, aga, sga))
	colnames(wgbs.top.cpg)<-c('chr','start','end','name','score','strand')

	# Top CpGs from this FLD 
	if(my.project=="FLD.v3.aga.low.vs.normal.pappa"){
		group1<-samples(rnb.set)[rnb.sample.groups(rnb.set)$Condition$CPA]
		group2<-samples(rnb.set)[rnb.sample.groups(rnb.set)$Condition$LPA]
	}else if(my.project=="FLD.v3.pet.control"){
		group1<-samples(rnb.set)[rnb.sample.groups(rnb.set)$Condition$control]
		group2<-samples(rnb.set)[rnb.sample.groups(rnb.set)$Condition$PET]
	}else{
		group1<-samples(rnb.set)[rnb.sample.groups(rnb.set)$Condition$AGA]
		group2<-samples(rnb.set)[rnb.sample.groups(rnb.set)$Condition$SGA]
	}
	my.diff.met<-rnb.diff.met[["sites"]][,c("Chromosome","Start","End","Strand",group1,group2)] # isa 'data.frame': Diff Meth CpGs from this analysis

	# Top CpGs from this FLD and previous WGBS (top CpGs within NKX1-2)
	my.top.met<-merge(my.diff.met, wgbs.top.cpg, by.x=c("Chromosome","Start","End"), by.y=c("chr","start","end")) # CpGs of 10 diff meth from WGBS
	my.top.met<-my.top.met[order(my.top.met$Start),]
	rownames(my.top.met)<-my.top.met$name

	out.file.name<-file.path(run.dir, my.run.date, paste0(my.project, ".top.cpg.methylation"))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size
	#heatmap.2(t(as.matrix(my.top.met[my.top.met$score<=10,c(group1,group2)]*100)), col=hmcol, trace="none", dendrogram="none", key.xlab="% methylation", keysize=1, cexCol=0.8, cexRow=0.8, srtCol=30, Rowv=F, Colv=F, cellnote=round(t(as.matrix(my.top.met[my.top.met$score<=10,c(group1,group2)]*100)),2), notecol="black", notecex=0.7)
	heatmap.2(t(as.matrix(my.top.met[,c(group1,group2)]*100)), col=hmcol, trace="none", dendrogram="none", key.xlab="% methylation", keysize=1, cexCol=0.8, cexRow=0.8, srtCol=30, Rowv=F, Colv=F, cellnote=round(t(as.matrix(my.top.met[,c(group1,group2)]*100)),2), notecol="black", notecex=0.7)
	dev.off()

	# Top 10 CpGs from this FLD and previous WGBS (top 10 only)
	my.diff.top10.met<-merge(my.diff.met, wgbs.top.cpg[wgbs.top.cpg$score<=10,], by.x=c("Chromosome","Start","End"), by.y=c("chr","start","end")) # CpGs of 10 diff meth from WGBS
	my.diff.top10.met<-my.diff.top10.met[order(my.diff.top10.met$Start),]
	rownames(my.diff.top10.met)<-my.diff.top10.met$name

	out.file.name<-file.path(run.dir, my.run.date, paste0(my.project, ".top10.cpg.methylation"))
	tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size
	heatmap.2(t(as.matrix(my.diff.top10.met[,c(group1,group2)]*100)), col=hmcol, trace="none", dendrogram="none", key.xlab="% methylation", keysize=1, cexCol=0.8, cexRow=0.8, srtCol=30, Rowv=F, Colv=F, cellnote=round(t(as.matrix(my.diff.top10.met[,c(group1,group2)]*100)),2), notecol="black", notecex=0.7)
	dev.off()
}

##########################################################
# Old in-house plotting method for methylatiton profile ##
##########################################################
if(FALSE){
	plot.diff.met<-function(my.region, my.gene){
		my.gene.rnb.diff.met<-rnb.diff.met$my.region[rnb.diff.met$symbol==my.gene,]
		if(as.character(unique(my.gene.rnb.diff.met$Strand))=="+"){
			my.gene.rnb.diff.met<-my.gene.rnb.diff.met[order(my.gene.rnb.diff.met$enst, my.gene.rnb.diff.met$Start, decreasing=FALSE),] # sort by enst, then position
		}else if(as.character(unique(my.gene.rnb.diff.met$Strand))=="-"){
			MY.gene.rnb.diff.met<-my.gene.rnb.diff.met[order(my.gene.rnb.diff.met$enst, my.gene.rnb.diff.met$Start, decreasing=TRUE),] # sort by enst, then position
		}else{
			cat(paste0(my.gene, " has more than one strand. Plot stopped!\n"))
			return(1)
		}

		my.gene.rnb.diff.met$Chromosome<-droplevels(my.gene.rnb.diff.met$Chromosome)
		my.gene.rnb.diff.met$id<-droplevels(my.gene.rnb.diff.met$id)
		my.gene.rnb.diff.met$enst<-droplevels(my.gene.rnb.diff.met$enst)
		my.gene.rnb.diff.met$ensg<-droplevels(my.gene.rnb.diff.met$ensg)
		my.gene.rnb.diff.met$symbol<-droplevels(my.gene.rnb.diff.met$symbol)
		my.gene.rnb.diff.met$enst_biotype<-droplevels(my.gene.rnb.diff.met$enst_biotype)

		diff.met.enst.list<-split(my.gene.rnb.diff.met, my.gene.rnb.diff.met$enst)
		for(i in names(diff.met.enst.list)){
			dummy<-diff.met.enst.list[[i]]
			# get the index for tiling array
			# ENST00000451024.12 => 12
			tiling.index<-sapply(dummy$id, function(x) as.numeric(substr(as.character(x), nchar(as.character(dummy[1,]$enst))+2, nchar(as.character(x)))) ) # get the index
			#data.frame(index=tiling.index, perc.met=dummy$mean.mean.g1, group=as.factor("AGA"), col="red")

			# matplot: Plot the columns of one matrix against the columns of another
			matplot(dummy[,c("mean.mean.g1", "mean.mean.g2")], pch=20, xaxt="n", main=paste0("enstTiling500 across ",my.gene," (",i,")"),
					xlab="Tiling index", ylab="% Methylation", col=my.col, cex=0.5) # xaxt="n" disables x-axis
			lines(stats::lowess(dummy$mean.mean.g1), col=my.col["AGA"], lty=2, lwd=3) # lty=2 (----), #lwd: thickness
			lines(stats::lowess(dummy$mean.mean.g2), col=my.col["SGA"], lty=2, lwd=3)
			axis(1, at=1:nrow(dummy), lab=tiling.index)
			#legend("topright", legend=names(my.col), lty=2, lwd=3, cex=2, col=my.col)
		}
	} # end of 'plot.diff.met' function

	for(i in names(rnb.diff.met)){
		###################
		# diagnostic plot
		###################
		cat("Making diagonstic plots for ",i,"...\n")
		if(i=="sites"){
			plot(my.diff.met$mean.mean.diff, -log10(my.diff.met$diffmeth.p.val+0.01), main=i, xlab="%met(AGA)-%met(SGA)", col="#00000020", cex=0.1)
			plot(my.diff.met$mean.mean.diff, -log10(my.diff.met$diffmet.p.adj.fdr), main=i, xlab="%met(AGA)-%met(SGA)", col="#00000020", cex=0.1)
		}else{
			plot(my.diff.met$mean.mean.diff, -log10(my.diff.met$comb.p.val), main=i, xlab="%met(AGA)-%met(SGA)", col="#00000020", cex=0.1)
			plot(my.diff.met$mean.mean.diff, -log10(my.diff.met$com,b.p.adj.fdr), main=i, xlab="%met(AGA)-%met(SGA)", col="#00000020", cex=0.1)
		}
			
		plot(my.diff.met$mean.mean.diff, my.diff.met$combinedRank, main=i, xlab="%met(AGA)-%met(SGA)", col="#00000020", cex=0.1)

		if(i=="enstTiling500"){
			cat("Making gene methylation profiles...\n")
			plot.diff.met(i, 'NKX1-2')
			plot.diff.met(i, 'RP13-238F13.3')
			plot.diff.met(i, 'NBEA')
			plot.diff.met(i, 'MAB21L1')
			plot.diff.met(i, 'GRM7')
		}

		##########################
		# conventional filtering
		##########################
		#select<-abs(my.diff.met$mean.mean.diff)>=0.25 & my.diff.met$comb.p.adj.fdr<=0.05
		if(i=="sites"){
			select<-my.diff.met$diffmet.p.adj.fdr<=0.05
		}else{
			select<-my.diff.met$comb.p.adj.fdr<=0.05
		}
		top.diff.met<-my.diff.met[select,]
		top.diff.met<-top.diff.met[order(abs(top.diff.met$mean.mean.diff), decreasing=TRUE),]
	}
} # old in-house meth profiles via 'matplot' 

# P-value 
#CpG=data.frame(`nkx1.2`=c(19,171),`elsewhere`=c(7,10481238),row.names=c("binary","non-binary"))
#tiling500=data.frame(`nkx1.2`=c(3,5),`elsewhere`=c(97,2941628),row.names=c("top100","non-top100"))

###################################################
# Relationship between methylation and expression #
###################################################
if(FALSE){
	if(grepl("Boy.Girl",my.project)){
		## diff.emt and diff.exp
		degseq.res<-as.data.frame(res); degseq.res$ID<-rownames(deseq.res)
		dummy<-merge(rnb.diff.met[["promoters"]], deseq.res, by=c("ID"))
		select<-dummy$num.sites>=median(dummy$num.sites) & !is.na(dummy$log2FoldChange) # CpG at least this number & there should be a fold change
		plot(x=dummy[select,c("mean.mean.diff")], y=-dummy[select,c("log2FoldChange")], main="Methylation and Expression Difference between 2F/2M", xlab="mean methylation difference(F-M)", ylab="log2(Fold Change(F-M))", col="#00000020", pch=20, cex=0.5)

		select<-dummy$num.sites>=median(dummy$num.sites) & !is.na(dummy$log2FoldChange)  & dummy$pvalue<=0.1
		plot(x=dummy[select,c("mean.mean.diff")], y=-dummy[select,c("log2FoldChange")], main="Methylation and Expression Difference between 2F/2M", xlab="mean methylation difference(F-M)", ylab="log2(Fold Change(F-M))", col="#00000020", pch=20, cex=0.5)

		select<-dummy$num.sites>=median(dummy$num.sites) & !is.na(dummy$log2FoldChange)  & dummy$combinedRank <= quantile(dummy$combinedRank)[2]
		plot(x=dummy[select,c("mean.mean.diff")], y=-dummy[select,c("log2FoldChange")], main="Methylation and Expression Difference between 2F/2M", xlab="mean methylation difference(F-M)", ylab="log2(Fold Change(F-M))", col="#00000020", pch=20, cex=0.5)
	}
}

##
destroy(rnb.set) #The ff files behind an RnBeads object can be deleted completely from the hard disk by executing the destructor method:
if(exists("rnb.diff.obj")){destroy(rnb.diff.obj)} #The ff files behind an RnBeads object can be deleted completely from the hard disk by executing the destructor method:
warnings() # print out the warnings
logger.completed()
logger.close()

cat("All done\n")
