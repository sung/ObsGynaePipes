
TR_PREFIX="GRCh37"
library(data.table)
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")

####################
## Manhattan Plot ##
####################
library(qqman)

for(RNASEQ_TYPE  in c("TOTAL","SMALL")){
	source("~/Pipelines/bin/R/chrX/local.R") # defines time.stamp 

	#dt.exp.pops[chromosome_name %in%  my.chr.order]
	
	#get gene-name
	#fields <- c(my.id, "start_position")
	#dt.gene<-as.data.table(getBM(attributes = fields, filters = my.id, values = dt.exp.pops[chromosome_name %in%  my.chr.order, unique(get(my.id))], mart = myMart)) # is a data.table
	dt.foo<-merge(dt.exp.pops[!is.na(padj) & chromosome_name %in%  my.chr.order], dt.ensg[,.(ensembl_gene_id,start_position)])

	dt.man<-dt.foo[,list(`SNP`=paste0(chromosome_name, ".",start_position), `CHR`=as.numeric(ifelse(chromosome_name=="X",23,chromosome_name)),BP=start_position,P=padj+my.beta,RANK=ifelse(log2FoldChange>0,-log10(padj+my.beta), log10(padj+my.beta)))]
	if(RNASEQ_TYPE=="TOTAL"){
		my.ylim=c(dt.man[!is.na(RANK),min(RANK)]-5, dt.man[!is.na(RANK),max(RANK)]+10)
	}else{
		my.ylim=c(dt.man[!is.na(RANK),min(RANK)]-0.55, dt.man[!is.na(RANK),max(RANK)]+1)
	}

	##
	## Chromosome-wise
	##
	file.name<-file.path("~/results/chrX/Figures/Manhattan",paste("Manhattan",RNASEQ_TYPE,"tiff",sep="."))
	tiff(filename=file.name,width=10, height=7,units="in",res=300, compression = 'lzw') #A4 size 

	#file.name<-file.path("~/results/chrX/Figures/Manhattan",paste("Manhattan",RNASEQ_TYPE,"pdf",sep="."))
	#pdf(file=file.name, width=10, height=7, title="Manhattan")

	par(mar=c(5,6,1,1)) # c(bottom, left, top, right)
	manhattan(dt.man, p="RANK", logp=FALSE, ylab=expression("-log"[10]*"(P-value)"), genomewideline=FALSE, suggestiveline=FALSE, chrlabs=c(1:22, "X"), font.lab=2, ylim=my.ylim, cex=3, cex.lab=3, cex.axis=2, yaxt='n')
	abline(h=my.abline, col=c(cbPalette[2], "black",cbPalette[2]), lty=c(2,1,2),lwd=c(2,1,2))
	axis(2, at=seq(round(my.ylim[1],-1),round(my.ylim[2],-2),by=20),labels=abs(seq(round(my.ylim[1],-1),round(my.ylim[2],-2),by=20)), cex.axis=2)
	# https://gist.github.com/cboettig/709578
	# side = 3 means top margin
	# adj = 0 means left align
	# line= 1.2 moves the text up
	#mtext(my.label, side=3, adj=0, line=1.2, cex=1.5, font=2) # font=2 to bold face
	dev.off()

	##
	## chrX only (padj from DESeq2)
	##
	file.name<-file.path("~/results/chrX/Figures/Manhattan",paste("Manhattan",RNASEQ_TYPE,"chrX.tiff",sep="."))
	tiff(filename=file.name,width=10, height=7,units="in",res=300, compression = 'lzw')
	par(mar=c(5,5,1,1))
	manhattan(dt.man[CHR==23], p="RANK", logp=FALSE, ylab="-log10(P-value)", genomewideline=FALSE, suggestiveline=FALSE, chrlabs=c("X"), font.lab=2, ylim=my.ylim, cex=2, cex.lab=2, cex.axis=1.3, xlab="Chromosome X position (Mb)",yaxt='n')
	abline(h=my.abline, col=c(cbPalette[2], "black",cbPalette[2]), lty=c(2,1,2),lwd=c(2,1,2))
	axis(2, at=seq(round(my.ylim[1],-1),round(my.ylim[2],-2),by=20),labels=abs(seq(round(my.ylim[1],-1),round(my.ylim[2],-2),by=20)), cex.axis=1.5)
	#mtext("D", side=3, adj=0, line=1.2, cex=1.5, font=2) # font=2 to bold face
	dev.off()

	##
	## chrX only (padj corrected based on chrX only)
	##
	fields <- c(my.id, "start_position")
	dt.gene<-as.data.table(getBM(attributes = fields, filters = my.id, values = dt.exp.pops.chrX[, unique(get(my.id))], mart = myMart)) # is a data.table
	dt.foo<-merge(dt.deg.chrX, dt.gene)

	dt.man<-dt.foo[,list(`CHR`=23,BP=start_position,P=new.padj+my.beta,RANK=ifelse(log2FoldChange>0,-log10(new.padj+my.beta), log10(new.padj+my.beta)))]
	if(RNASEQ_TYPE=="TOTAL"){
		my.ylim=c(dt.man[!is.na(RANK),min(RANK)]-5, dt.man[!is.na(RANK),max(RANK)]+10)
	}else{
		my.ylim=c(dt.man[!is.na(RANK),min(RANK)]-0.55, dt.man[!is.na(RANK),max(RANK)]+1)
	}

	file.name<-file.path("~/results/chrX/Figures/Manhattan",paste("Manhattan",RNASEQ_TYPE,"chrX.new.padj.tiff",sep="."))
	tiff(filename=file.name,width=10, height=8,units="in",res=300, compression = 'lzw') #A4 size 
	par(mar=c(5,5,3,1))
	manhattan(dt.man, p="RANK", logp=FALSE, ylab="-log10(P-value)", genomewideline=FALSE, suggestiveline=FALSE, chrlabs=c("X"), font.lab=2, ylim=my.ylim, cex=2, cex.lab=2, cex.axis=1.3, xlab="Chromosome X position (Mb)", yaxt='n')
	abline(h=my.abline, col=c(cbPalette[2], "black",cbPalette[2]), lty=c(2,1,2),lwd=c(2,1,2))
	axis(2, at=seq(round(my.ylim[1],-1),round(my.ylim[2],-2),by=20),labels=abs(seq(round(my.ylim[1],-1),round(my.ylim[2],-2),by=20)), cex.axis=1.5)
	#mtext("B", side=3, adj=0, line=1.2, cex=1.5, font=2) # font=2 to bold face
	dev.off()
}
