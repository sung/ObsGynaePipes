#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 3/Mar/2016
# Last modified 23/Aug/2016

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")

################################
## RoadMap                    ##
## Schultz et al. Nature 2015 ##
## Tissue comparison          ##
## Female-Male %Met Difference##
################################
my.cpg.type="CG" # CG, CHH, CHG
# 8 tissues ("AD","AO","EG","FT","GA","PO","SB","SX"): (7/March/2016)
# PA: not available yet
# PT: for our Placenta
#for(my.tissue in avail.tissues){
#for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
for(my.tissue in c("PT")){
	file.name<-file.path("~/results/RoadMap/BS-Seq/MethDiff",paste(my.tissue,my.cpg.type,time.stamp,"pdf",sep="."))
	pdf(file=file.name, width=11.7, height=8.3, title=paste0(my.cpg.type," methylation difference by sex in ",my.tissue))

	## Load this RoadMap tissue
	dt.query=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

	################################
	# 1. aggregate by CPG context ##
	################################
	#cpg.contexts=c("Sites","Tiling5K","Genes","Promo_15.05","Promo_15.00","Promo_10.00","Promo_10.10","CPGi","CPGi-with-promo","CPGi-no-promo","CPGi-shores","Enhancer")
	cpg.contexts=c("Sites","Tiling5K","Genes","Promo_15.05","CPGi","CPGi-shores","Enhancer")
	cpg.chrX.list=list() # to save chrX only (key: cpg context)
	for(i in cpg.contexts){
		cat(paste0("Processing ", i, "...\n"))
		if(i=='Sites'){
			group.keys.list="V1,V3,V2,End"                   # same as above 
			dt.cpg.region<-dt.query[,list(diff.met=V5.x/V6.x-V5.y/V6.y, num.sites=.N),by=group.keys.list]
		}else{
			subject<-with(as.data.frame(get.region(i)), data.frame(Chromosome=seqnames, Start=start, End=end, Strand=strand)) #isa 'data.frame'

			cat("\tRemoving leading 'chr'...\n")
			#subject$Chromosome<-simplify2array(strsplit(as.character(subject$Chromosome), "chr"))[2,] # chrX = > X
			subject$Chromosome<-substr(subject$Chromosome, 4, 5) #chrX=>X, chr13=>13
			cat("\tconvert to data.table...\n")
			dt.subject<-as.data.table(subject) # isa data.table
			rm(subject)

			if(grepl('Gene',i) || grepl('Promo',i) || grepl('enst',i)){
				subject.keys=c("Chromosome","Strand","Start","End") # regions of interests 
				query.key=c("V1","V3","V2","End")                   # columns from dt.query (roadmap merged data)
				group.keys=c("V1","V3","Start","End")               # V3: strand
				group.keys.list="V1,V3,Start,End"                   # same as above 
			# cpgislands, tiling, genomeTiling500
			}else{
				subject.keys=c("Chromosome","Start","End")
				query.key=c("V1","V2","End")
				group.keys=c("V1","Strand","Start","End")
				group.keys.list="V1,Strand,Start,End"
			}
			setkeyv(dt.subject,subject.keys)

			cat("\tFinding CpGs overlapping with ",i,"...\n")
			system.time(dt.overlap<-foverlaps(dt.query, dt.subject, by.x=query.key, type="any", nomatch=0L))

			cat("\tAggregating %met...\n")
			dt.cpg.region<-dt.overlap[,list(diff.met=sum(V5.x)/sum(V6.x)-sum(V5.y)/sum(V6.y), num.sites=.N),by=group.keys.list]
		}# end if i=='Sites'
		setnames(dt.cpg.region, c("chr","strand","start","end","diff.met","num.sites"))

		if(FALSE){
			# density-line by Chr
			a1<-ggplot(dt.cpg.region, aes(diff.met*100)) + 
				geom_density(aes(colour=chr),alpha=.8) + 
				labs(x="% methylation difference (female - male)") + 
				geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
				ggtitle(paste0(my.cpg.type," ",tissue.list[[my.tissue]], ": ",i)) +
				theme_Publication()
			print(a1)

			# density-fill by Chr
			if(my.cpg.type=="CG"){
				a2<-ggplot(dt.cpg.region, aes(diff.met*100)) + 
					geom_density(aes(fill=chr),alpha=.2) + 
					labs(x="% methylation difference (female - male)") + 
					geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
					ggtitle(paste0(my.cpg.type," ",tissue.list[[my.tissue]], ": ",i)) +
					theme_Publication()
				print(a2)
			# chrX only for CH
			}else{
				a2<-ggplot(dt.cpg.region[chr=="X"], aes(diff.met*100)) + 
					geom_density(aes(colour=chr),alpha=.8) + 
					labs(x="% methylation difference (female - male)") + 
					geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
					ggtitle(paste0(my.cpg.type," ",tissue.list[[my.tissue]], ": ",i, "chrX only")) + 
					theme_Publication()
				print(a2)
			}
		}

		# Boxplot by Chr
		min.cpg=38
		a3<-ggplot(dt.cpg.region[num.sites>=min.cpg], aes(chr,diff.met*100)) + 
			geom_boxplot() + 
			labs(y="% methylation difference (female - male)") + 
			geom_hline(yintercept=0) + 
			ylim(c(-100,100)) +
			ggtitle(paste0(my.cpg.type," ",tissue.list[[my.tissue]], ": ",i)) + 
			theme_Publication()
		if(i=="Sites"){
			out.file.name<-file.path("~/results/RoadMap/BS-Seq/MethDiff",paste(my.tissue,my.cpg.type,"meth.diff.boxplot.sites",sep="."))
			tiff(filename=paste0(out.file.name,".tiff"),width=11.7, height=8.27 ,units="in",res=300, compression = 'lzw') #A4 size 
			print(a3)
			dev.off()
		}else{
			print(a3)
		}

		# avoid storing single CH sites
		# this is to avoid 'too large for hashing' error of rbindlist
		if(grepl("CH",my.cpg.type) & i=="Sites"){cat("pass storing dt.cpg.region\n")}else{cpg.chrX.list[[i]]=dt.cpg.region[chr=="X",]}
		#cpg.chrX.list[[i]]=dt.cpg.region[chr=="X",]
	}# end of cpg.contexts

	# add a column (cpg context)
	lapply(names(cpg.chrX.list), function(i) cpg.chrX.list[[i]][,type:=i,])

	dt.cpg.chrX<-rbindlist(cpg.chrX.list)

	# Hitogram-line by CpG context (chrX only) 
	p<-ggplot(dt.cpg.chrX, aes(diff.met*100)) + 
		geom_density(aes(colour=type)) + 
		labs(x="% methylation difference (female - male)") + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle(paste0(my.cpg.type," ",tissue.list[[my.tissue]], ": chrX")) +
		#scale_colour_manual(values=cbPalette, name="CpG\nContext") +
		scale_colour_Publication(name="CpG\nContext") +
		theme_Publication()
	print(p)

	# Hitogram-fill by CpG context (chrX only) 
	q<-ggplot(dt.cpg.chrX, aes(diff.met*100)) + 
		geom_density(aes(fill=type),alpha=.2) + 
		labs(x="% methylation difference (female - male)") + 
		geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
		ggtitle(paste0(my.cpg.type," ",tissue.list[[my.tissue]], ": chrX")) +
		scale_fill_Publication(name="CpG\nContext") +
		theme_Publication()
	print(q)

	# Boxplot by CpG context (chrX only) 
	r<-ggplot(dt.cpg.chrX, aes(type, diff.met*100)) + 
		geom_boxplot(aes(fill=type)) + 
		labs(y="% methylation difference (female - male)") + 
		geom_hline(yintercept=0) + 
		ggtitle(paste0(my.cpg.type," ",tissue.list[[my.tissue]], ": chrX")) +
		#scale_fill_manual(values=cbPalette, name="CpG\nContext") +
		scale_fill_Publication(name="CpG\nContext") +
		theme_Publication()
	print(r)

	dev.off()
	cat(paste0("Done for ", my.tissue, "...\n"))
}# end of my.tissue 


	###############
	## plot_grid ##
	###############
	library("cowplot") # only for R >= 3.1.2
	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Sites"]

	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Genes"]
	dt.meth=dt.meth.per.tissue[[my.tissue]][Region=="Promo_15.05"]

cat("All done\n")
