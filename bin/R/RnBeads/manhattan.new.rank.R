## Rank2
		my.diff.met<-rnb.diff.met[[i]]

		if(i=="sites"){select<-!is.na(my.diff.met$combinedRank2)}else{select<-!is.na(my.diff.met$combinedRank2) & my.diff.met$num.sites >= median(my.diff.met$num.sites)}
		upper<-log(max(my.diff.met[select,]$combinedRank2))

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
					RANK=upper-log(combinedRank2),
					NUM.SITE=1
				))

		my.filename=file.path(run.dir,paste0(my.project,'.',i,".manhattan.rank.pvalue.abs.meth"))
		tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
		if(grepl("Boy.Girl",my.project)){
			manhattan(dummy, p="RANK", logp=FALSE, ylab=paste0(round(upper,2),"-log(combinedRank)"), genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X")) 
		}else{
			manhattan(dummy, p="RANK", logp=FALSE, ylab=paste0(round(upper,2),"-log(combinedRank)"), genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X","Y")) 
		}
		dev.off()


## Rank3
		my.diff.met<-rnb.diff.met[[i]]

		if(i=="sites"){select<-!is.na(my.diff.met$combinedRank3)}else{select<-!is.na(my.diff.met$combinedRank3) & my.diff.met$num.sites >= median(my.diff.met$num.sites)}
		upper<-log(max(my.diff.met[select,]$combinedRank3))

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
					RANK=upper-log(combinedRank3),
					NUM.SITE=1
				))

		my.filename=file.path(run.dir,paste0(my.project,'.',i,".manhattan.rank.pvalue.rel.meth"))
		tiff(filename=paste0(my.filename,".tiff"),width=16,height=9,units="in",res=300, compression = 'lzw')
		if(grepl("Boy.Girl",my.project)){
			manhattan(dummy, p="RANK", logp=FALSE, ylab=paste0(round(upper,2),"-log(combinedRank)"), genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X")) 
		}else{
			manhattan(dummy, p="RANK", logp=FALSE, ylab=paste0(round(upper,2),"-log(combinedRank)"), genomewideline=FALSE, suggestiveline=FALSE, main=main.title, highlight=top.diff.met, chrlabs=c(1:22, "X","Y")) 
		}
		dev.off()
