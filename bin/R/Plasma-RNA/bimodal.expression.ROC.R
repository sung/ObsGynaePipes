#!/usr/bin/Rscript --vanilla
library(DESeq2)
TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/GTEx/local.R")  # 'getGTExFpkm' 

load("~/results/RNA-Seq/GTEx/DESeq2.1.18.1/GTEx.PT.DESeq2.RData")   # loads dl.gtex.deseq & dl.Fpkm
                                                                    # see bin/R/GTEx/local.R

##############################################
## Bimodal Distribution of Gene Expression  ##
## Observed from Placenta                   ##
## Same in Plasma?                          ##
##############################################
if(FALSE){
    # Bimodal Distribution Pattern in Plasma & Placenta
    if(TRUE){
        # see bin/R/Placentome/placenta.specific.expression.R
        my.Tau=0.99; my.pt.fpkm=0.1; my.pl.fpkm=0.1; my.pt.ratio=100; my.pl.ratio=100
        my.file=paste0("~/results/RNA-Seq/Placentome/CSV/placenta.specific.protein.coding.genes.FPKM.Tau.",my.Tau,".PT.FPKM.",my.pt.fpkm,".ratio.",my.pt.ratio,".csv")
        dt.pt.specific<-fread(my.file) 

        my.genes=c("ERAP2","HLA-C","PSG7") # known bimodal expression (originally reported by Paul Kirk?)
        this.col=c("ensembl_gene_id","hgnc_symbol","chromosome_name")
        dt.target<-rbind(
                         dt.pt.specific[Tau.plasma>0.99 & Plasma/meanFpkmGTEx>2][order(-Plasma)][,..this.col], # plasma-specific too (n=14)
                         dt.ensg[hgnc_symbol %in% my.genes, ..this.col] # bimodal expression genes (n=3)
                         )
        #dt.ensg[hgnc_symbol %in% my.genes, ..this.col]
        dt.target[,(this.col) := lapply(.SD,as.character), .SDcols=this.col]

        ##########################
        # 1. Plasma RNA-Seq Data #
        ##########################
        load("~/results/RNA-Seq/Plasma.2017/DESeq2.1.18.1/ALL.Salmon/deseq.ALL.Salmon.RData") # in-house Plasma (SMRT-Seq CloneTech)
        #load("~/results/RNA-Seq/Plasma.ILM.2017/DESeq2.1.18.1/ALL.Salmon/deseq.ALL.Salmon.RData") # illumina Plamsa (Swift)
        dds.pl<-dds
        dt.pl.samples<-data.table(as.data.frame(colData(dds.pl)), SampleID=rownames(colData(dds.pl)))

        # FPKM
        dt.pl.Fpkm<-merge(data.table(`ensembl_gene_id`=rownames(fpkm(dds.pl)), fpkm(dds.pl)), dt.target)
        dt.melt.pl.Fpkm<-melt.data.table(dt.pl.Fpkm, id.vars=this.col, variable.name="SampleID", value.name="FPKM")

        # FPKM of each samples
        dt.sample.pl.fpkm<-merge(dt.melt.pl.Fpkm, dt.pl.samples)

        ############################################
        # 2. Placenta Tissue (Breach) RNA-Seq Data #
        ############################################
        ## 2. Illumina Placenta Tissue RNA-Seq N=19 (F=9, M=10)
        #load("/home/ssg29/results/RNA-Seq/Boy.Girl.ILM.GRCh38/DESeq2.1.18.1/BR.Salmon/deseq.BR.Salmon.RData")
        ## 2. In-house (JD-BR) N=19 (F=9, M=10) ## sample saples from ILM placneta tissues (n=19)
        load("/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/DESeq2.1.18.1/BR.Salmon/deseq.BR.Salmon.RData")
        dds.br.pt<-dds
        dt.br.pt.samples<-data.table(as.data.frame(colData(dds.br.pt)), SampleID=rownames(colData(dds.br.pt)))

        # FPKM
        dt.br.pt.Fpkm<-merge(data.table(`ensembl_gene_id`=rownames(fpkm(dds.br.pt)), fpkm(dds.br.pt)), dt.target)
        dt.melt.br.pt.Fpkm<-melt.data.table(dt.br.pt.Fpkm, id.vars=this.col, variable.name="SampleID", value.name="FPKM")
        
        # FPKM of each samples
        dt.sample.br.pt.fpkm<-merge(dt.melt.br.pt.Fpkm, dt.br.pt.samples)

        #########################################
        # 3. Placenta Tissue (All) RNA-Seq Data #
        #########################################
        ## 3. In-house (all samples)
        load("/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/DESeq2.1.18.1/ALL/deseq.ALL.RData")
        dds.pt<-dds
        dt.pt.samples<-data.table(as.data.frame(colData(dds.pt)), SampleID=rownames(colData(dds.pt)))

        # FPKM
        dt.pt.Fpkm<-merge(data.table(`ensembl_gene_id`=rownames(fpkm(dds.pt)), fpkm(dds.pt)), dt.target)
        dt.melt.pt.Fpkm<-melt.data.table(dt.pt.Fpkm, id.vars=this.col, variable.name="SampleID", value.name="FPKM")
        
        # FPKM of each samples
        dt.sample.pt.fpkm<-merge(dt.melt.pt.Fpkm, dt.pt.samples)

        ############################
        ## Merge Plasma & Placenta #
        ############################
        dt.foo<-rbind(
            dt.sample.pl.fpkm[,.(ensembl_gene_id,hgnc_symbol,chromosome_name,CRN,Sex,Source,SampleID,GA,FPKM)],
            dt.sample.br.pt.fpkm[,.(ensembl_gene_id,hgnc_symbol,chromosome_name,CRN,Sex,Source="Placenta",SampleID,GA=36,FPKM)]
            )
        dt.pl.pt<-dcast.data.table(dt.foo, ensembl_gene_id+hgnc_symbol+chromosome_name+CRN+Sex~Source+GA,value.var="FPKM")[order(hgnc_symbol,-Placenta_36)] # long (row-based) to wide (column-based)

        #####################
        ### FPKM from GTEX ##
        #####################
        #dl.Fpkm from ~/results/RNA-Seq/GTEx/DESeq2.1.18.1/GTEx.PT.DESeq2.RData
        dl.Fpkm[["Placenta"]] <- data.table(`Tissue`="Placenta", `ensembl_gene_id`=rownames(fpkm(dds.pt)), fpkm(dds.pt)) # replace by all placenta tissue samples (previously AGA only)
        dl.Fpkm[["Plasma"]] <- data.table(`Tissue`="Plasma", `ensembl_gene_id`=rownames(fpkm(dds.pl)), fpkm(dds.pl))
        dt.fpkm.tissue<-getGTExFpkm(dt.target$ensembl_gene_id) # 'getGTExFpkm' defined bin/R/GTEx/local.R 
    }

    if(TRUE){
        ###############
        ## ROC Curve ##
        ###############
        library(pROC)
        #my.file.name<- "/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/Cluster/ALL.Salmon/ROC.by.GA.ALL.Salmon"
        my.file.name<- "/home/ssg29/results/RNA-Seq/Plasma.2017/Cluster/ALL.Salmon/ROC.by.GA.ALL.Salmon"
        pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=8.5, height=8.3, title="ROC of bimally expressed genes") # A4 size

        dt.median<-dt.melt.br.pt.Fpkm[,.(`median`=median(FPKM)),"hgnc_symbol"] # from Illumina Placenta Tissue RNA-Seq N=19 (F=9, M=10)
        dt.median[hgnc_symbol=="PSG7",median:=100] # manual threshold
        dt.median[hgnc_symbol=="ERAP2",median:=2] # manual threshold
        dt.median[hgnc_symbol=="HLA-C",median:=10] # manual threshold

        dl.test<-list()
        # set 'outcome'
        for(my.gene in dt.median$hgnc_symbol){
	        dl.test[[my.gene]]<-dt.pl.pt[hgnc_symbol==my.gene,.(hgnc_symbol,`outcome`=factor(ifelse(Placenta_36>dt.median[hgnc_symbol==my.gene]$median,"1","0")), CRN,Sex,`Placenta`=Placenta_36, `12WK`=Plasma_12, `20WK`=Plasma_20, `28WK`=Plasma_28, `36WK`=Plasma_36)]
        }
        dcast.data.table(rbindlist(dl.test)[!is.na(outcome),.N,"hgnc_symbol,outcome"], hgnc_symbol~outcome)
        # plot ROC
        my.gene="PSG7"
        for(my.gene in dt.median$hgnc_symbol){
            #my.roc.12<-plot.roc(dl.test[[my.gene]]$outcome, dl.test[[my.gene]]$`12WK`, main=paste("AUC of ",my.gene,"by GA (Swift)"), legacy.axes=TRUE, ci=T)
            my.roc.12<-plot.roc(dl.test[[my.gene]]$outcome, dl.test[[my.gene]]$`12WK`, main=paste("AUC of ",my.gene,"by GA (SMART-Seq, CloneTech)"), legacy.axes=TRUE, ci=T)
            my.line.20<-lines.roc(dl.test[[my.gene]]$outcome, dl.test[[my.gene]]$`20WK`, col=cbPalette[2])
            my.line.28<-lines.roc(dl.test[[my.gene]]$outcome, dl.test[[my.gene]]$`28WK`, col=cbPalette[3])
            my.line.36<-lines.roc(dl.test[[my.gene]]$outcome, dl.test[[my.gene]]$`36WK`, col=cbPalette[4])

            if(my.gene=="PSG7"){
                # Illumina setting
                #text(.9, .52, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.roc.12$auc)), 2)))
                #text(.78, .85, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.20$auc)), 2)), col=cbPalette[2])
                #text(.6, .62, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.28$auc)), 2)), col=cbPalette[3])
                #text(.6, .98, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.36$auc)), 2)), col=cbPalette[4])
                # SMART-Seq setting
                text(.9, .58, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.roc.12$auc)), 2)))
                text(.78, .37, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.20$auc)), 2)), col=cbPalette[2])
                text(.6, .62, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.28$auc)), 2)), col=cbPalette[3])
                text(.6, .3, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.36$auc)), 2)), col=cbPalette[4])
            }else{
                text(.1, .2, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.roc.12$auc)), 2)))
                text(.1, .23, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.20$auc)), 2)), col=cbPalette[2])
                text(.1, .26, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.28$auc)), 2)), col=cbPalette[3])
                text(.1, .29, labels=paste("AUC =", round(as.numeric(gsub("Area under the curve:","", my.line.36$auc)), 2)), col=cbPalette[4])
            }
            legend("bottomright", legend=c("12Wk","20Wk","28Wk","36Wk"), col=cbPalette2, lwd=2)
        }
        dev.off()
        # save outcome & FPKM
        #my.file.name<- "/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/Cluster/ALL.Salmon/ROC.by.GA.ALL.Salmon.bimodal.genes" # Illumina
        my.file.name<- "/home/ssg29/results/RNA-Seq/Plasma.2017/Cluster/ALL.Salmon/ROC.by.GA.ALL.Salmon.bimodal.genes" # Clonetech
        write.csv(rbindlist(dl.test), file=paste(my.file.name,format(Sys.time(),'%Y-%m-%d_%I%p'),'csv',sep='.'),row.names=F,quote=F)

        #my.file.name<- "/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/Cluster/ALL.Salmon/ROC.by.GA.ALL.Salmon.pt.specific" # Illumina
        #my.file.name<- "/home/ssg29/results/RNA-Seq/Plasma.2017/Cluster/ALL.Salmon/ROC.by.GA.ALL.Salmon.pt.specific" # CloneTech
        #write.csv(rbindlist(dl.pl.pt), file=paste(my.file.name,format(Sys.time(),'%Y-%m-%d_%I%p'),'csv',sep='.'),row.names=F,quote=F)
    }

    ############################################
    ## FPKM Density Plot & Tissue Comparision ##
    ## for placenta & plasma specific genes
    ############################################
    if(TRUE){
        my.file.name <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/Cluster/ALL/placenta.specific.gene.FPKM.density"
        pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title="FPKM density of placenta specific genes") # A4 size

        for(my.gene in dt.target$hgnc_symbol){
            p<-ggplot(dt.melt.pt.Fpkm[hgnc_symbol==my.gene], aes(log2(FPKM))) + geom_histogram(col="black",fill="grey") + facet_grid(~hgnc_symbol) + theme_Publication()
            print(p)
        }

        dt.foo<-merge(dt.fpkm.tissue[,`Source`:=ifelse(Tissue=="Placenta","Placenta",ifelse(Tissue=="Plasma","Plasma","GTEx"))],dt.target)
        dt.foo<-getGTExFpkm(dt.ensg[hgnc_symbol=="SMS"]$ensembl_gene_id)
        dt.foo[,`Source`:=ifelse(Tissue=="Placenta","Placenta",ifelse(Tissue=="Plasma","Plasma","GTEx"))]
        dt.foo<-merge(dt.foo,dt.ensg)
        for(my.gene in dt.target$hgnc_symbol){
            p<-ggplot(dt.foo[hgnc_symbol==my.gene], aes(Tissue,log2(FPKM))) +
                geom_boxplot(aes(fill=`Source`),outlier.shape=NA,alpha=.7) +
                facet_grid(~hgnc_symbol) +
                scale_fill_manual(values=c(`GTEx`=cbPalette[1],`Placenta`=cbPalette[2],`Plasma`=cbPalette[3])) + 
                theme_Publication() +
                theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none")
            print(p)
        }
        dev.off()
        ## SMS
        dt.foo<-getGTExFpkm(dt.ensg[hgnc_symbol=="SMS"]$ensembl_gene_id)
        dt.foo[,`Source`:=ifelse(Tissue=="Placenta","Placenta",ifelse(Tissue=="Plasma","Plasma","GTEx"))]
        dt.foo<-merge(dt.foo,dt.ensg)
        p<-ggplot(dt.foo[hgnc_symbol==my.gene], aes(Tissue,log2(FPKM))) +
            geom_boxplot(aes(fill=`Source`),outlier.shape=NA,alpha=.7) +
            facet_grid(~hgnc_symbol) +
            scale_fill_manual(values=c(`GTEx`=cbPalette[1],`Placenta`=cbPalette[2],`Plasma`=cbPalette[3])) + 
            theme_Publication() +
            theme(axis.text.x=element_text(angle=45, hjust=1),legend.position="none")
        print(p)
    }

    ###########
    ## GDF15 ##
    ###########
    ## GDF15 for Gordon
    if(FALSE){
        #my.ensg=c(`GDF15`="ENSG00000130513",`PSG7`="ENSG00000221878")
        dt.ase<-fread("~/results/RNA-Seq/Placentome/CSV/GDF15.JD.allele.count.csv") #N=72 (71 unique samples)
        #merge(dt.ase,dt.samples,by.x="SampleName",by.y="SampleID")
        my.target.sample<-merge(dt.ase,dt.samples,by.x="SampleName",by.y="SampleID")$SampleName# N=61

        dt.my.fpkm<-merge(dt.melt.Fpkm[hgnc_symbol=="GDF15"], dt.samples[Source!="FG"]) # JD only (n=184)
        dt.my.fpkm[,`rs1058587`:=ifelse(SampleName %in% my.target.sample,'present','absent')]
        dt.my.fpkm[,.N,rs1058587]

        #ggplot(dt.my.fpkm, aes(log2(FPKM))) + geom_density() + facet_grid(~rs1058587) + theme_Publication()
        my.pval=wilcox.test(log2(dt.my.fpkm[rs1058587=="absent"]$FPKM), log2(dt.my.fpkm[rs1058587=="present"]$FPKM))$p.value

        my.file.name <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/Cluster/ALL/ALL.sample.FPKM.by.SNP"
        pdf(file=paste(my.filename,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title=paste0(myProject,":",sampleType)) # A4 size
        ggplot(dt.my.fpkm, aes(log2(FPKM),fill=rs1058587)) + geom_density(alpha=.6) + ggtitle(paste("Density plot of GDF15 expression (p-value=",my.pval,")")) + theme_Publication()
        ggplot(dt.my.fpkm, aes(rs1058587,log2(FPKM),fill=rs1058587)) + geom_boxplot(alpha=.6) + ggtitle(paste("Boxplot of GDF15 expression (p-value=",my.pval,")")) + scale_x_discrete(labels=c("absent (n=123)","present (n=61)")) + theme_Publication()
        dev.off()
    }
}
