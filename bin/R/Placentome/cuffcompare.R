#!/usr/bin/Rscript --vanilla
# to parse '.tracking' file from the cuffcompare
# Sung Gong <sung@bio.cc>

minFpkm=0.1
TR_PREFIX='GRCh38' # GRCh37|GRCh38
ENS_VER=82 # re-construction is based on GRCh38.82 which is used for StringTie
mySource="Placentome" # "FG"
myCohort="POPS" # "CTLM", "CTLF", "POPS", "CTL", "PET", "SGA"
myProject=paste(mySource,myCohort,sep=".") # e.g. Placentome.CTL
top.dir=file.path("~/results/RNA-Seq",mySource) #top.dir="~/results/RNA-Seq/Placentome"

source("~/Pipelines/config/Annotation.R")
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/config/cufflink.R")
source("~/Pipelines/bin/R/Placentome/local.R")

#myTarget=ifelse(myCohort=="CTL", "CSMD1.XLOC_482237",ifelse(myCohort=="CTLF","CSMD1.XLOC_395153","CSMD1.XLOC_370353")) # CTLF|CTLM (GRCh37 only)
#cuffRData=file.path(top.dir, paste0("RData/",myProject,".",TR_PREFIX,".",ENS_VER,".cuffcompare.",myTarget,".RData"))
cuffRData=file.path(top.dir, paste0("RData/",myProject,".",TR_PREFIX,".",ENS_VER,".cuffcompare.RData"))

if(file.exists(cuffRData)){
	cat("loading Cuffcompare RData...\n")
	load (cuffRData) # 
	cat("A list named ", myProject, "loaded\n")
}else{
	# see bin/R/Cuffcompare/local.R
	prep.cuff()
}

######################
## Data Preparation ##
######################
myCuff<-get(myProject) # myCuff<-FG.AGA
eval(parse(text=paste('rm(', myProject, ')'))) # remove the original object from memory

dt.frag<-myCuff$Transfrag   # isa 'data.table' (or data.frame?) (per each transcript (TCONS_))
dt.frag.class<-merge(dt.frag, dt.cuffClass, by="class_code")

#write.csv(data.table(dcast(dt.frag.class, Code + desc + class ~ `Sample Frequency`, value.var="Transcript Number"))[order(class)], file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort, ".tr.count.csv")),row.names=F)

# No of Exon
#fields <- c("ensembl_transcript_id", "ensembl_exon_id")
#dt.iso <- data.table(getBM(attributes =fields, filters = fields[1], values = dt.frag[fpkm.iqr>minFpkm & evi.ratio>=.95 & class_code %in% c("j","c") ,unique(ensembl_transcript_id)], mart = myMart))
#dt.iso <- merge(dt.frag[fpkm.iqr>minFpkm & evi.ratio>=.95 & class_code %in% c("j","c")], dt.iso[,list("enst.exon.cnt"=.N),ensembl_transcript_id], by="ensembl_transcript_id")

#> head(dt.frag)
#TCONS_00000001 TCONS_00000001 XLOC_000005        CICP27       ENST00000442987          = 0.73673925 1.00000000
#TCONS_00000002 TCONS_00000002 XLOC_000004 RP4-669L17.10       ENST00000440038          o 0.25082400 0.01724138

gr.cuff.exons<-myCuff[["Gff"]] # per exon-wise
tcons.gr.list=split(gr.cuff.exons, gr.cuff.exons$transcript_id)

###############################
## Sample & FPKM Information ##
## FRPK from StringTie       ##
###############################
# library and barcode information for corresponding re-constructed transcript (i.e. tcons)
# see  ~/results/RNA-Seq/Placentome/Cuffcompare/POPS/README
#dt.tcon.samples<-fread("zcat ~/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.82.cuffcompare.tcons.samples.csv.gz", header=FALSE, col.names=c("tcon","Library","Barcode","strg"))
dt.tcon.samples<-fread("zcat ~/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.82.cuffcompare.tcons.samples.fpkm.csv.gz", header=FALSE, col.names=c("tcon","Library","Barcode","strg","fpkm","cov","len"))

if(TRUE){
	cat("making PDF file...\n")
	# pdf output filename 
	my.file.name<- file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort, ".ggplot")) # per transcript (TCON_)
	pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=paste0("Transcriptome reconstruction of ",myProject)) # A4 size

	cat("making ggplot of basic stat...\n")
	##################################
	## CSMD1 Re-construction by Sex ##
	##################################
	if(FALSE){
		# 1. 1st exon containing transcript only
		my.exon.start=ifelse(TR_PREFIX=="GRCh37",3040500,31833146)
		my.tcon<-unique(mcols(gr.cuff.exons[start(gr.cuff.exons)>my.exon.start-10,])$transcript_id)
		dt.frag[tcon %in% my.tcon]
		# save CSV 
		my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".",myTarget,".tr.reconstruction.csv"))
		write.csv(dt.frag[tcon %in% my.tcon], file=my.file)
		# write GTF
		my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".",myTarget,".tr.reconstruction.gtf"))
		export.gff(tcons.gr.list[names(tcons.gr.list) %in% my.tcon], my.file)

		myTr=ifelse(myCohort=="CTL", "TCONS_01546312", ifelse(myCohort=="CTLF","TCONS_01105448","TCONS_01100815")) # CTLF|CTLM (GRCh37 only)
		dt.frag[tcon==myTr]
		tcons.gr.list[[myTr]]
		myCuff[["Fpkm"]][myTr,]; myCuff[["Coverage"]][myTr,]; myCuff[["Length"]][myTr,]

		##############################
		##  class_code distribution ##
		##############################
		# 1. a pie chart of class_code by the ratio of cuff.transfrag occurence 
		p <-ggplot(dt.frag.class, aes(x=factor(1), fill=class_code)) + 
			geom_bar(width = 1) + 
			coord_polar(theta = "y") + 
			scale_fill_manual(values=cbPalette) +
			theme_Publication() +
			theme(
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				axis.ticks = element_blank(),
				axis.text.x=element_blank(),
				axis.text.y=element_blank()
				)
		print(p)
	}

	###################
	## General stat. ##
	###################
	if(FALSE){
		# 2. a barchart of class_code (same as above)
		#levels(dt.frag$class_code)=unname(unlist(cuffClass)[unlist(cuffClass) %in% levels(dt.frag$class_code)]) # does not work
		#setattr(dt.frag$class_code,"levels",unname(unlist(cuffClass)[unlist(cuffClass) %in% levels(dt.frag$class_code)])) # does not work
		p <-ggplot(dt.frag.class, aes(class_code)) + 
			geom_bar(aes(fill=class)) + 
			scale_fill_manual(name="Class",values=cbPalette) +
			scale_x_discrete(name='Transcript Types', limits=dt.cuffClass[class_code %in% unique(dt.frag$class_code),class_code], labels=dt.cuffClass[class_code %in% unique(dt.frag$class_code),desc]) +
			ggtitle(paste("Total No. of transcript=",nrow(dt.frag))) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1))
		print(p)

		p <-ggplot(dt.frag.class[fpkm.iqr>=minFpkm], aes(class_code)) + 
			geom_bar(aes(fill=class)) + 
			scale_fill_manual(name="Class",values=cbPalette) +
			scale_x_discrete(name='Transcript Types', limits=dt.cuffClass[class_code %in% unique(dt.frag$class_code),class_code], labels=dt.cuffClass[class_code %in% unique(dt.frag$class_code),desc]) +
			ggtitle(paste("Total No. of transcript=",nrow(dt.frag[fpkm.iqr>=minFpkm]),"(FPKM>=",minFpkm,")")) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1))
		print(p)

		# 3. histogram of log(fpkm.iqr)
		#p <-ggplot(dt.frag, aes(x=log(fpkm.iqr+0.01))) +  geom_histogram(colour="black", fill="white") + ggtitle("histogram of log(fpkm.iqr)") + theme_Publication()
		#print(p)

		# 4. histogram of evi.ratio
		#p <-ggplot(dt.frag, aes(x=evi.ratio)) +  geom_histogram(colour="black", fill="white") + theme_Publication()
		#print(p)
		
		# 5. class_code count 
		#dummy<-sapply(names(sampleFreq), function(i) table(dt.frag[fpkm.iqr >= minFpkm & evi.ratio >=sampleFreq[i],class_code])) # isa 'matrix'
		#dummy[unlist(cuffClass)[unlist(cuffClass) %in% rownames(dummy)],]

		###########################################
		## NO. of class_code by Sample Frequency ##
		## either 'Novel', 'Known' or 'Novel'
		## count class_code by sampleFreq
		###########################################
		cat("making class_code ...\n")
		dummy<-sapply(names(sampleFreq), function(i) table(dt.frag[fpkm.iqr >= minFpkm & evi.ratio >=sampleFreq[i],class_code])) # isa 'matrix'
		dummy<-reshape2::melt(dummy) # isa 'data.frame'
		colnames(dummy)<-c("Code","Sample Frequency","Transcript Number")
		#   Var1  Var2 value
		#1     = >=90% 31503
		#2     . >=90%     0
		#3     c >=90%   301
		#4     e >=90%     0
		#5     i >=90%   751
		dummy<-data.table(merge(dummy, dt.cuffClass, by.x="Code",by.y="class_code"))
		#dummy<-dummy[order(`Sample Frequency`,class,-`Transcript Number`)]
		dummy2<-dummy[,list(`Transcript Number`=sum(`Transcript Number`)),list(`Sample Frequency`,class)][order(`Sample Frequency`,class)]

		p <-ggplot(dummy2, aes(x=`Sample Frequency`, y=`Transcript Number`, fill=class)) + 
			geom_bar(position="dodge", stat="identity", width=0.5) +  
			scale_fill_manual(name="Type",values=cbPalette) +
			ggtitle(paste0("No. of reconstructed transcript from ",myCohort,"(FPKM>=",minFpkm,")")) +
			theme_Publication()	
		print(p)

		p <-ggplot(dt.frag.class[fpkm.iqr>=minFpkm & evi.ratio>=0.25], aes(class_code)) + 
			geom_bar(aes(fill=class)) + 
			scale_fill_manual(name="Class",values=cbPalette) +
			scale_x_discrete(name='Transcript Types', limits=dt.cuffClass[class_code %in% unique(dt.frag$class_code),class_code], labels=dt.cuffClass[class_code %in% unique(dt.frag$class_code),desc]) +
			ggtitle(paste("Total No. of transcript=",nrow(dt.frag[fpkm.iqr>=minFpkm & evi.ratio>=0.25]),"(FPKM>=",minFpkm,"& Frequency>=25%)")) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1))
		print(p)

		p <-ggplot(dummy, aes(x=`Sample Frequency`, y=`Transcript Number`, fill=desc)) + 
			geom_bar(stat="identity", width=0.5) +  
			scale_fill_brewer(name="Class",palette="Set1") +
			ggtitle(paste0("Transcrotome reconstruction from ",myProject)) +
			theme_Publication()	
		print(p)

		#####################################
		## Frequency of transcript_biotype ##
		## For Complete Match              ##
		#####################################
		# myBiotype defined within Annotation.R 
		cat("making transcript_biotype frequency...\n")
		#select<-dt.frag$fpkm.iqr >= minFpkm & dt.frag$class_code=="="
		#dummy<-sapply(names(sampleFreq), function(i) table(merge(dt.frag[select & dt.frag$evi.ratio >= sampleFreq[i],c("tcon","transcript_biotype")], myBiotype, by.x="transcript_biotype", by.y="name")$biotype_group) ) # isa 'matrix'
		dummy<-sapply(names(sampleFreq), function(i) table(merge(as.data.frame(dt.frag[fpkm.iqr >= minFpkm & class_code=="=" & evi.ratio >= sampleFreq[i],.(tcon,transcript_biotype)]), myBiotype, by.x="transcript_biotype", by.y="name")$biotype_group) ) # isa 'matrix'
		foo<-reshape2::melt(dummy[order(rowSums(dummy),decreasing=T),]) # if dummy isa 'matrix'
		colnames(foo)<-c("biotype_group","Sample Frequency","Transcript Number")

		p <-ggplot(foo[foo$`Transcript Number`>0,], aes(x=`Sample Frequency`, y=`Transcript Number`, fill=biotype_group)) + 
			geom_bar(stat="identity", width=0.5) +
			scale_fill_manual(values=cbPalette) +
			ggtitle(paste0("Biotypes of 'Complete match' transcripts from  ",myProject)) +
			theme_Publication()	
		print(p)

		#######################
		## Novel Transcripts ##
		## by Sample Freq    ##
		## j, i, u only      ##
		#######################
		cat("counting novel transcript...\n")
		#dummy<-sapply(rev(seq(0, 1, by=.05)), function(i) table(dt.frag[ fpkm.iqr>=minFpkm & (class_code=="j" | class_code=="i"| class_code=="u") & evi.ratio>=i,class_code]) )[c("j","i","u"),]
		dt.frag[,class_code:=as.factor(class_code)]
		dummy<-sapply(rev(seq(0, 1, by=.05)), function(i) table(dt.frag[ fpkm.iqr>=minFpkm & evi.ratio>=i,class_code]) )
		colnames(dummy)<-paste0(">=",rev(seq(0, 1, by=.05))*100,"%") # isa matrix
		dummy<-reshape2::melt(dummy)
		colnames(dummy)<-c("Code","Sample Frequency","Transcript Number")
		dummy<-data.table(merge(dummy, dt.cuffClass, by.x="Code",by.y="class_code"))
		#dt.novel<-dummy[class=="Novel"]
		dt.novel<-dummy[Code %in% c("j","i","u")]
		dt.novel$Code<-droplevels(dt.novel$Code)
		#dt.novel$desc<-factor(dt.novel$desc, levels=names(unlist(cuffClass$Novel))) #character to factor
		dt.novel$desc<-factor(dt.novel$desc, levels=dt.cuffClass[class_code %in% c("j","i","u"),desc]) #character to factor

		# all bar-chart
		p<-ggplot(dt.novel, aes(x=`Sample Frequency`, y=`Transcript Number`, fill=desc)) +  
			geom_bar(stat="identity", position="dodge") +
			#scale_fill_manual(name="Class", values=cbPalette) +
			scale_fill_brewer(name="Class", palette="Set1") +
			ggtitle(paste0("Novel Transcript of ",mySource)) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1),legend.position=c(.2,.8))
			print(p)

		# 100-5 line
		p<-ggplot(dt.novel, aes(x=`Sample Frequency`, y=`Transcript Number`, colour=desc, group=desc)) +  
			geom_point(aes(shape=desc),size=4) + 
			geom_line(linetype="dashed") + 
			scale_color_brewer(name='Type', palette="Set1") +
			scale_shape_discrete(name='Type') +
			ggtitle(paste0("Novel Transcript of ",myProject)) +
			coord_cartesian(ylim = c(0, dt.novel[`Sample Frequency`==">=5%",max(`Transcript Number`)]+100)) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1),legend.position=c(.2,.8))
		print(p)

		# 100-50 line
		p<-ggplot(dt.novel[`Sample Frequency` %in% c(paste0(">=",rev(seq(0.5, 1, by=.05))*100,"%"))], aes(x=`Sample Frequency`, y=`Transcript Number`, colour=desc, group=desc)) + 
			geom_point(aes(shape=desc),size=4) + 
			geom_line(linetype="dashed") + 
			scale_color_brewer(name='Type', palette="Set1") +
			scale_shape_discrete(name='Type') +
			ggtitle(paste0("Novel Transcript of ",myProject)) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1))
		print(p)

		############################################
		## Novel Transcripts at leat 10% of cohort #
		############################################
		dt.novel1<-dt.frag.class[fpkm.iqr>=minFpkm & evi.ratio>=0.1 & class=="Novel",.N,desc]
		dt.novel1[,`Type`:="Transcript"]

		dt.novel2<-dt.frag.class[fpkm.iqr>=minFpkm & evi.ratio>=0.1 & class=="Novel",.(unique(hgnc_symbol)),desc][,.N,desc]
		dt.novel2[,`Type`:="Gene"]

		dt.novel<-rbind(dt.novel1,dt.novel2)
		dt.novel$desc<-factor(dt.novel$desc, levels=dt.novel1$desc) #character to factor

		p1<-ggplot(dt.novel, aes(desc, N)) +
			geom_bar(aes(fill=Type),stat="identity",position="dodge") +  
			ggtitle(paste(dt.novel1[,sum(N)], "novel transcripts (",dt.novel2[,sum(N)],"known genes)\npresent at least 10% of cohort",sep=" ")) +
			labs(x="Class",y="Number") +
			#scale_fill_brewer(name='Type', palette="Set4") +
			scale_fill_manual(values=cbPalette) +
			theme_Publication() +
			theme(axis.text.x=element_text(angle=45, hjust=1))

		# tiff output filename 
		my.file.name<- file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort, ".novel.10.tiff")) # per transcript (TCON_)
		tiff(filename=my.file.name,width=12,height=9,units="in",res=300, compression = 'lzw')
		print(p1)
		dev.off()
	}

	####################
	# Gene of Interest #
    # e.g. LEP, FLT1   #
	####################
	if(TRUE){
        myGene="LEP" # or LEP
		myDir=file.path(top.dir, paste("Cuffcompare",myCohort,myGene,sep="/"))
        system(paste0("mkdir -p ",myDir))

		my.tcon=dt.frag.class[hgnc_symbol==myGene & fpkm.iqr>=minFpkm & evi.ratio>=0.1 & class!="Artefact",tcon] # transcirpts of this gene

        # transcripts of this gene with minFPKM & minRatio
		dt.my.gene<-merge(
                          merge(dt.frag.class[tcon %in% my.tcon], dt.tcon.samples, by.x="tcon", by.y="tcon"), 
                          dt.pops, 
                          by.x=c("Library","Barcode"), by.y=c("Library","BarCode"))
        # transcripts of this gene without any filter
		#dt.my.gene<-merge(merge(dt.frag.class[hgnc_symbol==myGene],                                         dt.tcon.samples, by.x="tcon", by.y="tcon"), dt.pops, by.x=c("Library","Barcode"), by.y=c("Library","BarCode"))

		dcast(dt.my.gene, CRN+Cohort ~ tcon)
        unique(dt.my.gene[,.(tcon,class_code,CRN,Cohort)])[,.N,"tcon,Cohort"]
        dcast(unique(dt.my.gene[,.(tcon,class_code,CRN,Cohort)])[,.N,"tcon,Cohort"], tcon~Cohort)

		# dt.pops from bin/R/Placentome/local.R
        if(myGene=="LEP"){
            # samples without re-constructed LEP transcript (N=21), why? (probably filtered out by StringTie)
            # also check this out: 
            # for i in `cat ~/results/RNA-Seq/Placentome/Meta/POPS.stringtie.GRCh38.82.gtf.list`; do  printf "$i\t"; bedtools intersect -a $i -b data/LEP/LEP.GRCh38.bed | wc -l; done | awk '$2==0{print $0}' | wc -l
			# for Ulla
			dt.foo<-dcast(dt.my.gene, CRN+Cohort ~ tcon, value.var="fpkm", fun.aggregate=length, fill=0)[order(Cohort)]
			my.file=file.path(myDir,paste(myCohort,myGene,TR_PREFIX,ENS_VER,"cuffcompare.CRN.csv",sep="."))
			write.csv(dt.foo, file=my.file)
            # dt.foo[TCONS_00323031>1]
        }

        # Write transcripts in CSV
		my.file=file.path(myDir,paste(myCohort,myGene,TR_PREFIX,ENS_VER,"reconstructed.tr.csv",sep="."))
		write.csv(dt.frag.class[hgnc_symbol==myGene & fpkm.iqr>=minFpkm & evi.ratio>=0.1], file=my.file,quote=F,row.names=F) 
        ## Write GTF
		my.file=file.path(myDir,paste(myCohort,myGene,TR_PREFIX,ENS_VER,"reconstructed.tr.gtf",sep="."))
		export.gff(unlist(tcons.gr.list[names(tcons.gr.list) %in% my.tcon]), my.file)
		my.file=file.path(myDir,paste(myCohort,myGene,TR_PREFIX,ENS_VER,"reconstructed.tr.gff3",sep="."))
		export.gff(tcons.gr.list[names(tcons.gr.list) %in% my.tcon], my.file)

		if(FALSE){
			# Avg. of FPKM (vertical, row: samples, col: transcripts)
			dt.fpkm<-dcast(dt.my.gene, SampleName+Cohort ~ tcon, value.var="fpkm", fun.aggregate=mean, fill=0)[order(Cohort)] # sample (row) order by cohort
			dt.my.gene.fpkm<-as.matrix(dt.fpkm[,3:eval(length(names(dt.fpkm)))]); rownames(dt.my.gene.fpkm)=dt.fpkm$SampleName
			dt.my.gene.fpkm<-dt.my.gene.fpkm[,names(sort(colSums(dt.my.gene.fpkm!=0), decreasing=T))] # transcript (column) order by the sample freq.

			ann_samples<-data.frame(`Type`=dt.fpkm$Cohort,row.names=rownames(dt.my.gene.fpkm))
			ann_colors<-list(`Type`=c(`CTL`=cbPalette2[1], `SGA`=cbPalette2[2], `PET`=cbPalette2[3]))
            if(capabilities()[["tiff"]]){
			    my.file=file.path(myDir,paste(myCohort,myGene,TR_PREFIX,ENS_VER,"cuffcompare.fpkm.tiff",sep="."))
            }else{
			    my.file=file.path(myDir,paste(myCohort,myGene,TR_PREFIX,ENS_VER,"cuffcompare.fpkm.jpeg",sep="."))
            }
			pheatmap::pheatmap(log(dt.my.gene.fpkm*100+1), cluster_rows=F, cluster_cols=F, annotation_row=ann_samples, annotation_colors=ann_colors, fontsize_row=5, fontsize_col=13)
			pheatmap::pheatmap(log(dt.my.gene.fpkm*100+1), cluster_rows=F, cluster_cols=F, annotation_row=ann_samples, annotation_colors=ann_colors, fontsize_row=5, fontsize_col=13, filename=my.file)

			my.file=file.path(myDir,paste(myCohort,myGene,TR_PREFIX,ENS_VER,"cuffcompare.fpkm.csv",sep="."))
			write.csv(dt.my.gene.fpkm, file=my.file)
		}else{
			# Avg. of FPKM (horizontal, row: transcripts, col: samples))
		}

	}
	####################################################
	## Quantification based on reference vs. novel.10 ##
	## Pick up a gene(s) for Steve                    ##
	####################################################
	if(TRUE){
		dt.known.novel.10<-dt.frag.class[fpkm.iqr>=minFpkm & evi.ratio>=0.1 & class %in% c("Novel","Known")]
		dt.rpt<-dt.frag.class[fpkm.iqr>=minFpkm & evi.ratio>=0.1 & class %in% c("Novel")]

		dt.comp<-list()
		dt.comp[["SLX-11368"]]<-fread("~/results/SLX-11368.Homo_sapiens.v1/Salmon/SLX-11368.NoIndex.quant.POPS.GRCh38.82.ensg.ref.vs.novel.10.txt", col.names=c("ensembl_gene_id","ref","rpt"))
		dt.comp[["SLX-11369"]]<-fread("~/results/SLX-11369.Homo_sapiens.v1/Salmon/SLX-11369.NoIndex.quant.POPS.GRCh38.82.ensg.ref.vs.novel.10.txt", col.names=c("ensembl_gene_id","ref","rpt"))
		#dt.comp[["SLX-11368"]]<-fread("~/results/SLX-11368.Homo_sapiens.v1/Salmon/SLX-11368.NoIndex.quant.POPS.GRCh38.82.ensg.ref.vs.known.novel.10.txt", col.names=c("ensembl_gene_id","ref","rpt"))
		#dt.comp[["SLX-11369"]]<-fread("~/results/SLX-11369.Homo_sapiens.v1/Salmon/SLX-11369.NoIndex.quant.POPS.GRCh38.82.ensg.ref.vs.known.novel.10.txt", col.names=c("ensembl_gene_id","ref","rpt"))

		min.read=100; min.fc=5
		lapply(dt.comp, function(i) i[ref>min.read & rpt>min.read ,fc:=rpt/ref])
		dt.foo<-merge(dt.comp[["SLX-11368"]][ref>min.read & rpt>min.read & fc>min.fc], dt.comp[["SLX-11369"]][ref>min.read & rpt>min.read & fc>min.fc], by="ensembl_gene_id")
		dt.bar<-merge(merge(dt.foo, dt.ensg[,.(ensembl_gene_id, hgnc_symbol, gene_biotype)]), dt.rpt, by="hgnc_symbol")
		dt.bar[exon.cnt>1 & transcript_biotype=="protein_coding"][order((fc.x+fc.y)/2)]

		# xloc overlapping with only 1 known reference gene
		merge(dt.bar[exon.cnt>1 & transcript_biotype=="protein_coding"],  dt.frag.class[,.(`gene_cnt`=length(unique(hgnc_symbol))),xloc][gene_cnt==1], by="xloc")[order(fc.x+fc.y)]
	}
	##############################
	## NO. of isoforms per gene ##
	## of exact protein coding  ##
	##############################
	if(FALSE){
		cat("counting isoform per gene...\n")
		library(plyr)
		#select<-dt.frag$fpkm.iqr >= minFpkm & dt.frag$class_code=="=" & dt.frag$transcript_biotype=="protein_coding"
		dummy<-sapply(names(sampleFreq), function(i){
											foo<-as.data.frame(dt.frag[fpkm.iqr >= minFpkm & class_code=="=" & transcript_biotype=="protein_coding" & evi.ratio >= sampleFreq[i]]) # isa 'data.frame'
											#for(j in colnames(foo)[1:5])(foo[[j]]<-droplevels(foo[[j]])) # droplevel 
											cbind(
												`Sample Frequency`=i, 
												`Transcript Number`=plyr::daply(foo, .(xloc), function(df) length(unique(df$tcon))),
												`ENST Number`=plyr::daply(foo, .(xloc), function(df) length(unique(df$ensembl_transcript_id))),
												`HGNC Number`=plyr::daply(foo, .(xloc), function(df) length(unique(df$hgnc_symbol))),
												`Avg FPKM`=plyr::daply(foo, .(xloc), function(df) mean(df$fpkm.iqr,na.rm=T))
											) # count distinct tcon by xloc
										} 
				) # isa 'list'
		bar<-as.data.frame(do.call("rbind",dummy)) 

		bar[,1]<-factor(bar[,1], levels=names(sampleFreq)) # reorder factor
		for(i in 2:5)(bar[[i]]<-as.numeric(as.character(bar[[i]]))) # a factor to numeric must be proceeded by factor to character, then numeric

		p<-ggplot(as.data.frame(bar), aes(x=`Sample Frequency`, y=`Transcript Number`, fill=`Sample Frequency`)) + geom_boxplot() + ggtitle("Avg. number per 'protein_coding' gene of 'Complete match'") + theme_Publication()
		print(p)
		p<-ggplot(as.data.frame(bar), aes(x=`Sample Frequency`, y=`ENST Number`, fill=`Sample Frequency`)) + geom_boxplot() + ggtitle("Avg. number per 'protein_coding' gene of 'Complete match'") + theme_Publication()
		print(p)
		p<-ggplot(as.data.frame(bar), aes(x=`Sample Frequency`, y=`HGNC Number`, fill=`Sample Frequency`)) + geom_boxplot() + ggtitle("Avg. number per 'protein_coding' gene of 'Complete match'") + theme_Publication()
		print(p)
	}


	#############
	# Write GFF #
	#############
	# already done
	if(FALSE){
		my.freq=10
		my.dt.frag<-dt.frag.class[fpkm.iqr>=minFpkm & evi.ratio>=my.freq*0.01 & class %in% c("Novel","Known")]
		my.dt.frag<-dt.frag.class[fpkm.iqr>=minFpkm & evi.ratio>=my.freq*0.01 & class %in% c("Novel")]
		my.tcon<-my.dt.frag[,tcon]
		# save CSV 
		my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".known.novel.",my.freq,".tr.reconstruction.csv"))
		write.csv(dt.frag.class[tcon %in% my.tcon], file=my.file)

		# write GTF
		my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".known.novel.",my.freq,".tr.reconstruction.gtf"))
		export.gff(tcons.gr.list[names(tcons.gr.list) %in% my.tcon], my.file)

		# to use 'salmon quant -g'
		# 1. TCON XLOC
		dt.tcon.xloc<-my.dt.frag[,.(tcon,xloc)]
		my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".known.novel.",my.freq,".tr.reconstruction.xloc.txt"))
		write.table(dt.tcon.xloc, sep="\t", my.file, row.names=F, col.names=F, quote=F)

		# 2. TCON hgnc_symbol 
		# 2.1 novel.x
		dt.tcon.hgnc_symbol<-unique(my.dt.frag[,.(tcon,hgnc_symbol)])
		my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".known.novel.",my.freq,".tr.reconstruction.hgnc_symbol.txt"))
		write.table(dt.tcon.hgnc_symbol, sep="\t", my.file, row.names=F, col.names=F, quote=F)
		# 2.2 all novel
		#dt.tcon.hgnc_symbol<-unique(dt.frag.class[,.(tcon,hgnc_symbol)])
		#my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".tr.reconstruction.hgnc_symbol.txt"))
		#write.table(dt.tcon.hgnc_symbol, sep="\t", my.file, row.names=F, col.names=F, quote=F)

		# 3. TCON ensembl_gene_id
		# 3.1. novel.x
		dt.tcon.ensg<-unique(merge(my.dt.frag[,.(tcon,ensembl_transcript_id)], dt.enst[,.(ensembl_transcript_id, ensembl_gene_id)])[,.(tcon,ensembl_gene_id)])
		my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".known.novel.",my.freq,".tr.reconstruction.ensg.txt"))
		write.table(dt.tcon.ensg, sep="\t", my.file, row.names=F, col.names=F, quote=F)
		# 3.2. all novel
		#dt.tcon.ensg<-unique(merge(dt.frag.class[,.(tcon,ensembl_transcript_id)], dt.enst[,.(ensembl_transcript_id, ensembl_gene_id)])[,.(tcon,ensembl_gene_id)])
		#my.file=file.path(top.dir, paste0("Cuffcompare/",myCohort,"/",myCohort,".",TR_PREFIX,".tr.reconstruction.ensg.txt"))
		#write.table(dt.tcon.ensg, sep="\t", my.file, row.names=F, col.names=F, quote=F)
	}

	####################################
	## Transcript Frequency by Cohort ## 
	####################################
	if(FALSE){
		dt.frag[hgnc_symbol=="CSMD1" & evi.ratio*ncol(myCuff$Fpkm)>=5][dt.frag[hgnc_symbol=="CSMD1" & evi.ratio*ncol(myCuff$Fpkm)>=5, order(evi.ratio,decreasing=T)]]
		dt.sample.cohort[Sample %in% names(myCuff$Fpkm["TCONS_00352614",][!is.na(myCuff$Fpkm["TCONS_00352614",])]),][,.N,"Cohort"]

		chisq.test(
			dt.sample.cohort[Sample %in% names(myCuff$Fpkm["TCONS_00352614",][!is.na(myCuff$Fpkm["TCONS_00352614",])]),][,.N/table(!is.na(myCuff$Fpkm["TCONS_00352614",]))["TRUE"],"Cohort"][,V1],
			p=dt.cnt.cohort[,.N, by="Cohort,Library,BarCode"][,list(freq=.N/ncol(myCuff$Fpkm)),"Cohort"][,freq]
		)

		chisq.test(
			dt.sample.cohort[Sample %in% names(myCuff$Fpkm["TCONS_00352614",][!is.na(myCuff$Fpkm["TCONS_00352614",])]),][,.N,"Cohort"][,N],
			p=dt.cnt.cohort[,.N, by="Cohort,Library,BarCode"][,list(cnt=.N),"Cohort"][,cnt],
			rescale.p = TRUE
		)

		prop.test(
			dt.sample.cohort[Sample %in% names(myCuff$Fpkm["TCONS_00352614",][!is.na(myCuff$Fpkm["TCONS_00352614",])]),][,.N,"Cohort"][,N],
			dt.cnt.cohort[,.N, by="Cohort,Library,BarCode"][,list(cnt=.N),"Cohort"][,cnt]
		)
	}

	################
	## RepeatMask ##
	################
	if(FALSE){
		gr.repeat<-import.bed("~/data/Annotation/UCSC/Homo_sapiens/GRCh38/rmsk.bed.gz",genome="hg38")
		select<-seqnames(gr.repeat) %in% extractSeqlevels("Homo_sapiens","UCSC") # canonical chromosome only (chr1, chr2.. chrY, chrM)
		gr.repeat<-gr.repeat[select,] # filter non-canonical chr
		seqlevels(gr.repeat)<-seqlevelsInUse(gr.repeat) # drop un-used seqlevels
		seqlevels(gr.repeat)=mapSeqlevels( seqlevels(gr.repeat),style="NCBI" ) # chrX=>X
		#seqlevels(gr.repeat)=sub("chr","",seqlevels(gr.repeat)) # same as above

		gr.novel<-gr.cuff.exons[gr.cuff.exons$transcript_id %in% as.character(dt.frag[class_code=="i" & evi.ratio>=0.9,tcon]),]
		grl.novel<-tcons.gr.list[as.character(dt.frag[class_code=="i" & evi.ratio>=0.9,tcon])]

		intersect(gr.novel, gr.repeat)
		lapply(grl.novel, function(i) intersect(i, gr.repeat))
	}

	###########################
	## GO Analysis via goseq ##
	###########################
	if(FALSE){
		do.go.cuff()
	}
	dev.off()
}

cat("All is done\n")
