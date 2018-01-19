library(data.table)
TR_PREFIX='GRCh37'
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R") 

# WGoxBS
process.PT.bed<-function(){
	cat("Processing PT:\n")
	dt.meta<-fread("/scratch/obsgynae/POPS/results/RnBeads/Meta/cpg.sample.annotation.sbs.csv", stringsAsFactor=FALSE)

	# 2 AGA females (merged/aggregated by chr/position/strand)
	cat("processing female...\n")
	print(system.time(dl.female<-lapply(dt.meta[Condition=='AGA' & FetalSex=="F",BedFile], function(i){fread(paste("zcat",i),sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)
	names(dl.female)<-dt.meta[Condition=='AGA' & FetalSex=="F",SampleName]
	dt.female<-rbindlist(dl.female)[,list(`num.c`=round(sum(V4/100*V5),2),`depth`=sum(V5), .N),.(V1,V3,V6)] # V1:chr, V3:pos, V6:str
	#	  V1       V3 V6 num.c depth N
	#1: chrX 21957917  +    25    66 2
	#2: chrX 21957918  -     1     2 1
	#3: chrX 21957958  +     9    26 2
	setnames(dt.female, c("V1","V2","V3","V5","V6","N")) # V1:chr, V2:Pos, V3:Str, V5:cnt-met, V6:depth-cov, N: No. of sample having this CpG
	dt.female[,V1:=substr(V1,4,5)] # chrX=>X chrMT=>MT

	# 2 AGA males (merged/aggregated by chr/position/strand)
	cat("processing female...\n")
	print(system.time(dl.male<-lapply(dt.meta[Condition=='AGA' & FetalSex=="M",BedFile], function(i){fread(paste("zcat",i),sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)
	names(dl.male)<-dt.meta[Condition=='AGA' & FetalSex=="M",SampleName]
	dt.male<-rbindlist(dl.male)[,list(`num.c`=round(sum(V4/100*V5),2),`depth`=sum(V5), .N),.(V1,V3,V6)]
	setnames(dt.male, c("V1","V2","V3","V5","V6","N"))
	dt.male[,V1:=substr(V1,4,5)]

	## Save
	cat("Saving dt.female and dt.male...\n")
	print(system.time(save(dt.female, dt.male, file="~/results/RoadMap/BS-Seq/RData/PT.CG.dt.female.male.RData")))

	PT.CG.dt.merged<-merge(dt.female, dt.male,by=c("V1","V2","V3"))
	#> PT.CG.dt.merged
    #	V1       V2 V3 V5.x V6.x N.x V5.y V6.y N.y
	#1:  1    10469  +    0    7   2    2    3   2
	#2:  1    10470  -    1    9   2    0    8   2
	#3:  1    10471  +    1   14   2    3    8   2
	save(PT.CG.dt.merged, file="~/data/RoadMap/BS-Seq/Schultz2015/PT.CG.dt.merged.RData")
}

##############################
## % methylation @ Promoter ##
##############################
# see ~/Pipelines/bin/R/RoadMap/local.R for details
prep.dt.meth.chrX.per.tissue<-function(my.min.doc=10){
	min.doc=my.min.doc
	my.cpg.type="CG"
	dt.tss<-dt.ensg[chromosome_name=="X",.(ensembl_gene_id,hgnc_symbol,chromosome_name,strand,TSS=ifelse(strand=="+",start_position,end_position))] # TSS of gene
	dt.promoter<-dt.tss[,.(ensembl_gene_id,hgnc_symbol,chromosome_name,strand,start=TSS-1000,end=TSS+1000)] # regions of interest (e.g. promoter)
	dt.gene<-dt.ensg[chromosome_name=="X",.(ensembl_gene_id,hgnc_symbol,chromosome_name,strand,start=start_position,end=end_position)]
	dt.meth.per.region=list() # to save %met by cpg contexts
	counter=1
	for(my.tissue in c("PT","PT.SS.Tech","PT.SS.Bio")){ # for all tissue types
		my.new.RData=file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,my.cpg.type,"dt.merged.RData",sep="."))
		cat(paste0("Loading ", my.tissue, " ", my.cpg.type, "...\n"))
		print(system.time(load(my.new.RData))) # load list of data.table
		cat(paste(my.tissue,my.cpg.type,"dt.merged loaded\n",sep="."))
		this.tissue<-get(paste(my.tissue,my.cpg.type,"dt.merged",sep=".")) # isa list of data.table (STL2, STL3)
		eval(parse(text=paste0('rm(', my.tissue,'.',my.cpg.type,'.dt.merged)'))) # remove the original object from memory
		## Merge all chromosomes
		if(!grepl("PT",my.tissue)){
			dt.query<-this.tissue[["chrX"]]
		}else{
			dt.query<-this.tissue[V1=="X"]
		}
		rm(this.tissue)
		dt.query[,End:=V2] # add this column

		for(my.region in c("Genes","Promo_10.10")){
			# strand-aware (also "ensembl_gene_id" available by default)
			if(grepl('Gene',my.region) || grepl('Promo',my.region)){
				if(grepl('Gene',my.region)){dt.subject=dt.gene}else{dt.subject=dt.promoter}
				subject.keys=c("chromosome_name","strand","start","end") # regions of interests 
				query.key=c("V1","V3","V2","End")                   # columns from dt.query (dt.merged data) V1:chr, V2:pos, V3:str
				group.keys.list="ensembl_gene_id,hgnc_symbol,V1,V3,start,end" 
			# non-strand-aware (CPGi, CPGi-shores,Enhancer, etc)
			}else{
				subject.keys=c("chromosome_name","start","end") # regions of interests 
				query.key=c("V1","V2","End")
				group.keys.list="ensembl_gene_id,hgnc_symbol,V1,start,end"
			}
			setkeyv(dt.subject,subject.keys)

			cat("\tFinding CpGs overlapping with Promoter...\n")
			#by.x, by.y: A vector of column names (or numbers) to compute the overlap joins. 
			#			The last two columns in both by.x and by.y should each correspond to the start and end interval columns in x and y respectively. 
			#			And the start column should always be <= end column. If x is keyed, by.x is equal to key(x), else key(y). by.y defaults to key(y).
			system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L)) # chrX only

			cat("\tAggregating %met...\n")
			dt.meth<-dt.overlap[V6.x>=min.doc & V6.y>=min.doc,list(sum(V5.x),sum(V6.x),round(sum(V5.x)/sum(V6.x)*100,2),sum(V5.y),sum(V6.y),round(sum(V5.y)/sum(V6.y)*100,2),num.sites=.N,my.region,my.tissue),by=group.keys.list]
			#> dt.meth 
			#   ensembl_gene_id V1 V3     Start       End   V1   V2  V3  V4 num.sites           i my.tissue
			#1: ENSG00000251848  X  +   2728441   2730440   22   29   8  12         1 Promo_15.03        SB
			setnames(dt.meth,c("ensembl_gene_id","hgnc_symbol","chr","strand","start","end","female.c","female.cov","female.met","male.c","male.cov","male.met","num.sites","region","tissue"))
			dt.meth.per.region[[counter]]=dt.meth # save this cpg context
			counter=counter+1
		} # end of region
	}# end of tissue
	dt.meth.chrX.per.tissue<-rbindlist(dt.meth.per.region)
	save(dt.meth.chrX.per.tissue, file=file.path("~/results/RoadMap/X-inactivation/RData",paste("dt.meth.chrX.per.tissue",my.cpg.type,min.doc,"RData",sep=".")))
}
prep.dt.meth.chrX.per.tissue(my.min.doc=1)

#########
## SMS ##
#########
#dt.meth[hgnc_symbol=="SMS"]
#PT.CG.dt.merged[V1=="X" & V3=="+" & V2>=dt.tss[hgnc_symbol=="SMS",TSS-1000] & V2<=dt.tss[hgnc_symbol=="SMS",TSS+1000] & V6.x>=min.doc & V6.y>=min.doc,list(sum(V5.x), sum(V6.x), round(sum(V5.x)/sum(V6.x)*100,2), sum(V5.y), sum(V6.y), round(sum(V5.y)/sum(V6.y)*100,2),.N)]
#PT.SS.Bio.CG.dt.merged[V1=="X" & V3=="+" & V2>=dt.tss[hgnc_symbol=="SMS",TSS-1000] & V2<=dt.tss[hgnc_symbol=="SMS",TSS+1000] & V6.x>=min.doc & V6.y>=min.doc,list(sum(V5.x), sum(V6.x), round(sum(V5.x)/sum(V6.x)*100,2), sum(V5.y), sum(V6.y), round(sum(V5.y)/sum(V6.y)*100,2),.N)]
#PT.SS.Tech.CG.dt.merged[V1=="X" & V3=="+" & V2>=dt.tss[hgnc_symbol=="SMS",TSS-1000] & V2<=dt.tss[hgnc_symbol=="SMS",TSS+1000] & V6.x>=min.doc & V6.y>=min.doc,list(sum(V5.x), sum(V6.x), round(sum(V5.x)/sum(V6.x)*100,2), sum(V5.y), sum(V6.y), round(sum(V5.y)/sum(V6.y)*100,2),.N)]
# WGoxBS (SMS only)
process.PT.SMS.bed<-function(my.cpg.type){
	cat(paste0("Processing PT:", my.cpg.type,"\n"))
	dt.meta<-fread("~/Pipelines/data/SMS/meta.WGoxBS.csv", stringsAsFactor=FALSE)

	# 2 AGA females
	print(system.time(dl.female<-lapply(dt.meta[Condition=='AGA' & FetalSex=="F",BedFile], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)
	names(dl.female)<-dt.meta[Condition=='AGA' & FetalSex=="F",SampleName]
	dt.female<-rbindlist(dl.female)[,list(`num.c`=round(sum(V4/100*V5),2),`depth`=sum(V5), .N),.(V1,V3,V6)]
	#	  V1       V3 V6 num.c depth N
	#1: chrX 21957917  +    25    66 2
	#2: chrX 21957918  -     1     2 1
	#3: chrX 21957958  +     9    26 2
	setnames(dt.female, c("V1","V2","V3","V5","V6","N"))
	# V1:chr, V2:Pos, V3:Str, V5:cnt-met, V6:depth-cov, N: No. of sample having this CpG

	# 2 AGA males
	print(system.time(dl.male<-lapply(dt.meta[Condition=='AGA' & FetalSex=="M",BedFile], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)
	names(dl.male)<-dt.meta[Condition=='AGA' & FetalSex=="M",SampleName]
	dt.male<-rbindlist(dl.male)[,list(`num.c`=round(sum(V4/100*V5),2),`depth`=sum(V5), .N),.(V1,V3,V6)]
	setnames(dt.male, c("V1","V2","V3","V5","V6","N"))

	dt.female[V2>=dt.ensg[hgnc_symbol=="SMS",start_position-1000] & V2<=dt.ensg[hgnc_symbol=="SMS",start_position+1000] & V6>=5, list(sum(V5), sum(V6), round(sum(V5)/sum(V6)*100,2),.N)]
	dt.male[V2>=dt.ensg[hgnc_symbol=="SMS",start_position-1000] & V2<=dt.ensg[hgnc_symbol=="SMS",start_position+1000] & V6>=5, list(sum(V5), sum(V6), round(sum(V5)/sum(V6)*100,2),.N)]
}

#########################
## Methylation Profile ##
## via GGBIO           ##
#########################
library(ggbio)
my.point.size=2.5; my.alpha=0.5; min.doc=10
my.cpg.type="CG"
my.target.gene="SMS"

my.upflank=100 # up from TSS
my.downflank=100 # down form TES

## Target Gene
my.dt.ensg=dt.ensg[hgnc_symbol==my.target.gene]
my.gr.ensg<-GenomicRanges::GRanges(seqnames=Rle(my.dt.ensg$chromosome_name),IRanges(start=my.dt.ensg$start,end=my.dt.ensg$end,names="SMS"),strand=my.dt.target$strand)
#my.gr.ensg<-GenomicRanges::GRanges(seqnames=Rle(paste0("chr",my.dt.ensg$chromosome_name)),IRanges(start=my.dt.ensg$start,end=my.dt.ensg$end,names="SMS"),strand=my.dt.target$strand)

## Target Region 
my.dt.target=dt.ensg[hgnc_symbol==my.target.gene, .(chromosome_name,start=ifelse(strand=="+",start_position-1000-my.upflank,end_position+1000+my.upflank),end=ifelse(strand=="+",start_position+1000+my.downflank,end_position-1000-my.downflank),strand)]
my.gr.target<-GenomicRanges::GRanges(seqnames=Rle(my.dt.target$chromosome_name),IRanges(start=my.dt.target$start,end=my.dt.target$end,names="SMS"),strand=my.dt.target$strand)
#my.gr.target<-GenomicRanges::GRanges(seqnames=Rle(paste0("chr",my.dt.target$chromosome_name)),IRanges(start=my.dt.target$start,end=my.dt.target$end,names="SMS"),strand=my.dt.target$strand)
#my.gr.target<-GenomicRanges::promoters(gr.ensg[mcols(gr.ensg)$gene_name==my.target.gene], upstream=1000+my.upflank, downstream=1000+my.downflank)
#seqlevels(my.target)=mapSeqlevels(seqlevels(my.gr.target), style="UCSC") # X=>chrX

plot.meth.profile<-funcion(my.gr.target,min.doc=10){
	cat("\tconvert to data.table...\n")
	dt.target<-data.table(as.data.frame(my.gr.target)) 
	#dt.target[,seqnames:=substr(seqnames,4,5)] # chrX=>X chrMT=>MT
	#> data.table(as.data.frame(my.gr.target))
	#   seqnames    start      end width strand
	#1:     X 21957591 21959791  2201      +
	setkeyv(dt.target,c("seqnames","strand","start","end")) # regions of interests 

	query.key=c("V1","V3","V2","End") # # columns from dt.query (dt.merged data) V1:chr, V3:str,  V2:pos, End:(same as V2)

	# R/3.1.1 (this ggbio (autoplot specifically) conflicts with ggplot which may be >2.0)
	# use R/3.2.1 (this version of ggbio works, but color/fill does does not work :()
	## see bin/R/RnBeads/SureSelect.ggplot.R

	####################
	# transcript track #
	####################
	#p.tr.reduced<-autoplot(hg19ensGene, which=my.gr.target,stat="reduce") + theme_alignment() # hg19ensGene defined in 'Annotation.R'
	p.tr.reduced<-autoplot(my.gr.ensg, which=my.gr.target,stat="reduce")+ theme_alignment()

	## Load Placenta methylation
	dt.pt.locus<-list()
	#for(my.tissue in c("PT","PT.SS.Tech","PT.SS.Bio")){
	for(my.tissue in c("PT")){
		## Load this tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
		cat("\tFinding CpGs overlap...\n")
		system.time(dt.overlap<-foverlaps(dt.query[V1==dt.target$seqnames], dt.target, by.x=query.key, type="any", nomatch=0L))
		#> dt.overlap
		#   V1 V3    start      end width       V2 V5.x V6.x N.x V5.y V6.y N.y      End
		#1:  X  + 21957591 21959791  2201 21958345    1    1   1    0    1   1 21958345
		if(nrow(dt.overlap)){
			dt.meth<-rbind(dt.overlap[V6.x>=min.doc & V6.y>=min.doc,.(V1,V3,V2,round(V5.x/V6.x*100,2),DoC=V6.x,"Female",my.tissue)], 
						   dt.overlap[V6.x>=min.doc & V6.y>=min.doc,.(V1,V3,V2,round(V5.y/V6.y*100,2),DoC=V6.y,"Male",my.tissue)]) # isa 'data.table'
			setnames(dt.meth,c("chr","strand","start","methylation_level","CpG_Depth","Gender","Tissue"))
			dt.pt.locus[[my.tissue]]<-dt.meth
		}
	}# end of my.tissue 
	## Load RoadMap Tissue
	dt.roadmap.locus<-list()
	for(my.tissue in c("AD","AO","EG","FT","GA","PO","SB","SX")){ # except 'PT': placenta
		## Load this tissue
		dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)
		cat("\tFinding CpGs overlap...\n")
		system.time(dt.overlap<-foverlaps(dt.query[V1==dt.target$seqnames], dt.target, by.x=query.key, type="any", nomatch=0L))
		#> dt.overlap
		#   V1    Start      End Strand       V2 V3  V4 V7 V5.x V6.x V5.y V6.y    i.End
		#1:  X  2746554  2748235      *  2746590  - CGA  1    6   17    2   10  2746590
		#2:  X  2835950  2836236      *  2836235  + CGG  1   14   14    9   10  2836235
		if(nrow(dt.overlap)){
			dt.meth<-rbind(dt.overlap[V6.x>=min.doc & V6.y>=min.doc,.(V1,V3,V2,round(V5.x/V6.x*100,2),DoC=V6.x,"Female",my.tissue)], 
						   dt.overlap[V6.x>=min.doc & V6.y>=min.doc,.(V1,V3,V2,round(V5.y/V6.y*100,2),DoC=V6.y,"Male",my.tissue)]) # isa 'data.table'
			setnames(dt.meth,c("chr","strand","start","methylation_level","CpG_Depth","Gender","Tissue"))
			dt.roadmap.locus[[my.tissue]]<-dt.meth
		}
	}# end of my.tissue 

	# Placenta gender-specific
	p.pt<-ggplot(rbindlist(dt.pt.locus), aes(start, methylation_level, col=Gender)) + geom_point(size=my.point.size,alpha=my.alpha) + ylim(0, 100) + 
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="none") 
	p.pt.facet<-ggplot(rbindlist(dt.pt.locus), aes(start, methylation_level, col=Gender)) + geom_point(size=my.point.size,alpha=my.alpha) + ylim(0, 100) + 
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Tissue ~ .)
	# Roadmap gender-specific
	p.roadmap<-ggplot(rbindlist(dt.roadmap.locus), aes(start, methylation_level, col=Gender)) + geom_point(size=my.point.size,alpha=my.alpha) + ylim(0, 100) + 
				geom_smooth(se=FALSE) + 
				scale_colour_manual(name="Sex", values=my.col[["Gender"]]) + 
				labs(y="% Methylation") +
				theme_Publication() + 
				theme(legend.position="none") + 
				facet_grid(Tissue ~ .)

	# PT and RoadMap 
	my.file.name<-file.path("~/results/RoadMap/X-inactivation/MethProfile",paste0(my.target.gene,".meth.profile.tiff"))
	tiff(filename=my.file.name,width=10, height=15 ,units="in",res=300, compression = 'lzw') #A4 size 
	print(tracks(`Tr`=p.tr.reduced, `Placenta`=p.pt, `Schultz et al`=p.roadmap, heights = c(1,3.5,28), scale.height=2.5, label.text.cex = 1.5, xlim=my.gr.target))
	dev.off()
}

load.my.tissue.dt.merged<-function(my.tissue,my.cpg.type){
	## Load this RoadMap tissue
	my.new.RData<-file.path("~/data/RoadMap/BS-Seq/Schultz2015", paste(my.tissue,my.cpg.type,"dt.merged.RData",sep="."))
	cat(paste0("Loading ", my.tissue, " ", my.cpg.type, "...\n"))
	if(file.exists(my.new.RData)){
		print(system.time(load(my.new.RData))) # load list of data.table
	}else{
		stop(paste0(my.new.RData, " not found\n"))
	}
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
	dt.query[,End:=V2] # add this column
	#>dt.query
	#   V1    V2 V3  V4 V7 V5.x V6.x V5.y V6.y   End
	#1: 10 66269  + CGG  1   12   12   44   51 66269
	cat("some integer columns as numeric...\n") # why?
	this.col=c("V5.x","V6.x","V5.y","V6.y")
	print(system.time(dt.query[, (this.col) := lapply(.SD,as.numeric), .SDcols=this.col]))
	return(dt.query)
}

