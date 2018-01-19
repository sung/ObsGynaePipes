#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 3/Mar/2016
# Last modified 9/May/2016

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")

my.cpg.type="CHG" # CG, CHH, CHG

file.name<-file.path("~/results/RoadMap/CH/Figures",paste("meth.diff.by.gender.all.tissue",my.cpg.type,time.stamp,"pdf",sep="."))
pdf(file=file.name, width=11.7, height=8.3, title="Methylation difference by gender across various tissue types")

my.RData=file.path("~/results/RoadMap/CH/RData",paste(my.cpg.type,"dt.dmr.list.RData",sep="."))
if(file.exists(my.RData)){
	load(my.RData)
	cat("dt.dmr.list loaded\n")
}else{
	dt.dmr.list<-list()
	for(my.tissue in c("PT",avail.tissues[!avail.tissues %in% "PA"])){
		cat(paste0("Processing ", my.tissue, "...\n"))
		this.RData=file.path("~/results/RoadMap/CH/RData",paste(my.tissue,my.cpg.type,"dt.dmr.rank.RData",sep="."))
		if(!file.exists(this.RData)){
			stop(paste0(this.RData," not found\n"))
		}else{
			print(system.time(load(this.RData)))
			cat(paste(this.RData," loaded\n")) # dt.dmr.per.region
			cat("dt.dmr.per.region loaded\n")
		}

		#########################
		## Filter by num.sites ##
		#########################
		new.dt.dmr.per.region=lapply(dt.dmr.per.region, function(i) i[num.sites>=i[,median(num.sites)]]) # at least median num.sites

		# Boxplot of % meth.diff by region of this tissue
		p1<-ggplot(rbindlist(new.dt.dmr.per.region), aes(region, (met.f-met.m)*100)) + 
			geom_boxplot(aes(fill=region),outlier.shape=NA) + 
			labs(y="% Methylation Difference (female - male)") + 
			geom_hline(yintercept=0) + 
			coord_cartesian(ylim=c(-1.2,1.2)) + 
			ggtitle(paste(my.tissue, my.cpg.type))
		print(p1)

		# Boxplot of % meth by region of this tissue
		p2<-ggplot(rbindlist(new.dt.dmr.per.region), aes(region, (met)*100)) + 
			geom_boxplot(aes(fill=region),outlier.shape=NA) + 
			labs(y="% Methylation Difference (female - male)") + 
			coord_cartesian(ylim=c(0,1.5)) + 
			ggtitle(paste(my.tissue, my.cpg.type))
		print(p2)

		dt.dmr.list[[my.tissue]]<-rbindlist(new.dt.dmr.per.region)
		save(dt.dmr.list,file=my.RData)
	}
}

# Boxplot of % meth.diff by region across tissue types 
p3<-ggplot(rbindlist(dt.dmr.list), aes(tissue, (met.f-met.m)*100)) + 
	geom_boxplot(aes(fill=region),outlier.shape=NA) + 
	labs(y="% Methylation Difference (female - male)") + 
	geom_hline(yintercept=0) + 
	scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts")) +
	coord_cartesian(ylim=c(-1.2,1.2)) + 
	ggtitle(my.cpg.type)
#print(p3)

# Boxplot of % meth by region across tissue types 
p4<-ggplot(rbindlist(dt.dmr.list), aes(tissue, (met)*100)) + 
	geom_boxplot(aes(fill=region),outlier.shape=NA) + 
	labs(y=paste("%",my.cpg.type,"Methylation")) + 
	scale_fill_manual(values=cbPalette, name=paste0(my.cpg.type,"\nregional contexts")) +
	coord_cartesian(ylim=c(0,1.5)) + 
	ggtitle(my.cpg.type)
print(p4)

dev.off()
cat("All done\n")
