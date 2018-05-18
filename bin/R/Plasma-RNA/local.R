# dl.pl.deseq : $Salmon$FG, $Salmon$ILM, $FtCount$FG, $FtCount$ILM
# dl.pl.ga.deseq : $GA12$FG, $GA12$ILM ... $GA36$FG, $GA36$ILM
# dl.pt.deseq : $Salmon$FG, $Salmon$ILM, $FtCount$FG, $FtCount$ILM
# dl.cnt : (readcount) $ILM, $POPS
prep.dt.deg<-function(deseq.RData){
    if(file.exists(deseq.RData)){
        cat("loading DESeq2 RData...\n")
        cat(paste("reading from ", deseq.RData,"\n"))
        load(deseq.RData)
        cat("dds, res, loaded\n")
    }else{
		stop(paste(deseq.RData, " not found"))
    }
    if(priorInfo(res)$type=="none"){
        cat("Applying lfcShrink()...\n")
        print(system.time(res <- lfcShrink(dds,contrast=c(my.contrast,levels(colData(dds)[[my.contrast]])[2],levels(colData(dds)[[my.contrast]])[1]),res=res,parallel=TRUE)))
        contrast=c(my.contrast,levels(colData(dds)[[my.contrast]])[2],levels(colData(dds)[[my.contrast]])[1])
    }
    # create RData for the first time 
    dt.res<-data.table(data.frame(res)); dt.res[[my.filter]]<-rownames(res);
    rn=rownames(colData(dds)) #as.character(samples[["SampleName"]]) # samples$SampleName (or samples[,c("SampleName")])
    rn.control=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1==0,])
    rn.case=rownames(colData(dds)[as.numeric(colData(dds)[[my.contrast]])-1!=0,])
    ddsFpkm <-fpkm(dds) # isa 'matrix'
    ddsFpm <-fpm(dds) # isa 'matrix'

    foo<-data.table(meanFpkm=rowMeans(ddsFpkm,na.rm=T),
                    meanFpkm.case=rowMeans(ddsFpkm[,rn.case],na.rm=T),
                    meanFpkm.control=rowMeans(ddsFpkm[,rn.control],na.rm=T),
                    meanFpm=rowMeans(ddsFpm,na.rm=T),
                    meanFpm.case=rowMeans(ddsFpm[,rn.case],na.rm=T),
                    meanFpm.control=rowMeans(ddsFpm[,rn.control],na.rm=T))
    foo[[my.filter]]<-rownames(ddsFpkm)
    deseq.anno=merge(dt.res, foo, all.x=TRUE)
    deseq.anno<-merge(deseq.anno, dt.ensg[,.(ensembl_gene_id,chromosome_name,hgnc_symbol,gene_biotype)], all.x=TRUE)[order(pvalue)]
    return(deseq.anno)
}#end of prep.dt.deg

###############################################
## illumina Plasma RNA-Seq Samples           ##
## regardless of GA                          ##
## quantification by Salmon and featureCount ##
## 'dl.pl.deseq'                             ##
###############################################
if(TRUE){
    my.RData.file="/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/DESeq2.1.18.1/dl.pl.deseq.FG.vs.ILM.RData"
    if(file.exists(my.RData.file)){
        load(my.RData.file) # load 'dl.pl.deseq'
    }else{
        dl.pl.deseq<-list() 
        ######################
        ## 1. RNA-Seq by FG ##
        ######################
        TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
        ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
        source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
        ## 1-a. quantification by featureCounts
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Plasma.2017/DESeq2.1.18.1/ALL/deseq.ALL.RData"
        dl.pl.deseq[["FtCount"]][["FG"]]<-prep.dt.deg(deseq.RData)
        cat("FtCount-FG loaded\n")
        ## 1-b. quantification by Salmon 
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Plasma.2017/DESeq2.1.18.1/ALL.Salmon/deseq.ALL.Salmon.RData"
        dl.pl.deseq[["Salmon"]][["FG"]]<-prep.dt.deg(deseq.RData)
        cat("Salmon-FG loaded\n")

        ############################
        ## 2. RNA-Seq by Illumina ##
        ############################
        ## 2-a. quantification by Salmon 
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/DESeq2.1.18.1/ALL.Salmon/deseq.ALL.Salmon.RData"
        dl.pl.deseq[["Salmon"]][["ILM"]]<-prep.dt.deg(deseq.RData)
        cat("Salmon-ILM loaded\n")

        ## 1-b. quantification by featureCounts
        TR_PREFIX='GRCh37' # GRCh37 | GRCh38 
        ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
        source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/DESeq2.1.18.1/ALL/deseq.ALL.RData"
        dl.pl.deseq[["FtCount"]][["ILM"]]<-prep.dt.deg(deseq.RData)
        cat("FtCount-ILM loaded\n")

        ##
        save(dl.pl.deseq, file=my.RData.file)
        cat("dl.pl.deseq saved\n")
        ##
    }

    my.Quant="FtCount"
    my.Quant="Salmon"
    lapply(dl.pl.deseq[[my.Quant]], nrow)
    lapply(dl.pl.deseq[[my.Quant]], function(i) nrow(i[padj<0.05]))
    lapply(dl.pl.deseq[[my.Quant]], function(i) nrow(i[padj<0.05 & meanFpm>=0.25])) # 1 in 4 million reads
    lapply(dl.pl.deseq[[my.Quant]], function(i) nrow(i[padj<0.05 & abs(log2FoldChange)>=1])) # at least 2 fold-change
    lapply(dl.pl.deseq[[my.Quant]], function(i) nrow(i[log2FoldChange>0 & padj<0.05])) # Female > Male
    lapply(dl.pl.deseq[[my.Quant]], function(i) nrow(i[log2FoldChange<=0 & padj<0.05])) # Female <= Male

    lapply(dl.pl.deseq[[my.Quant]], function(i) i[chromosome_name=="Y" & padj<0.05,.(hgnc_symbol,chromosome_name,padj,log2FoldChange,baseMean,meanFpkm,meanFpm)])
    lapply(dl.pl.deseq[[my.Quant]], function(i) i[chromosome_name=="X" & padj<0.05,.(hgnc_symbol,chromosome_name,padj,log2FoldChange,baseMean,meanFpkm,meanFpm)])
}

#######################################
## illumina Plasma RNA-Seq Samples   ##
## controlled by GA                  ##
## quantification by Salmon          ##
## 'dl.pl.ga.deseq'                  ##
#######################################
if(TRUE){
    my.RData.file="/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/DESeq2.1.18.1/dl.pl.ga.deseq.FG.vs.ILM.RData"
    if(file.exists(my.RData.file)){
        load(my.RData.file) # load 'dl.pl.ga.deseq'
    }else{
        TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
        ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
        source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg

        dl.pl.ga.deseq<-list()
        ######################
        ## 1. RNA-Seq by FG ##
        ######################
        my.source=list(`Plasma.2017`="FG",`Plasma.ILM.2017`="ILM")
        for(j in paste("GA",c(12, 20, 28, 36),sep=".")){
            for(i in c("Plasma.2017", "Plasma.ILM.2017")){
                deseq.RData <- file.path("/home/ssg29/results/RNA-Seq",paste0(i,"/DESeq2.1.18.1/",j, "/deseq.",j,".RData"))
                dl.pl.ga.deseq[[j]][[my.source[[i]]]]<-prep.dt.deg(deseq.RData)
            }
        }
        ##
        save(dl.pl.ga.deseq, file=my.RData.file)
        cat("dl.pl.ga.deseq saved\n")
        ##
    }
}


#####################################
## Placenta tissue RNA-Seq Samples ##
## 'dl.pt.deseq'
#####################################
if(TRUE){
    my.RData.file="/home/ssg29/results/RNA-Seq/Boy.Girl.ILM.GRCh37/DESeq2.1.18.1/dl.pt.deseq.FG.vs.ILM.RData"
    if(file.exists(my.RData.file)){
        load(my.RData.file) # load 'dl.pt.deseq'
    }else{
        dl.pt.deseq<-list() 
        ######################
        ## 1. RNA-Seq by FG ##
        ######################
        TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
        ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
        source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
        ## 1-a. quantification by featureCounts
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/DESeq2.1.18.1/BR/deseq.BR.RData"
        dl.pt.deseq[["FtCount"]][["FG"]]<-prep.dt.deg(deseq.RData)
        cat("FtCount-FG loaded\n")

        ## 1-b. quantification by Salmon 
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/DESeq2.1.18.1/BR.Salmon/deseq.BR.Salmon.RData"
        dl.pt.deseq[["Salmon"]][["FG"]]<-prep.dt.deg(deseq.RData)
        cat("Salmon-FG loaded\n")

        ## 1-b. quantification by Salmon 
        #cat("Salmon-FG loaded\n")
        ############################
        ## 2. RNA-Seq by Illumina ##
        ############################
        ## 2-a. quantification by Salmon 
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.ILM.GRCh38/DESeq2.1.18.1/BR.Salmon/deseq.BR.Salmon.RData"
        dl.pt.deseq[["Salmon"]][["ILM"]]<-prep.dt.deg(deseq.RData)
        cat("Salmon-ILM loaded\n")

        ## 2-b. quantification by featureCounts
        TR_PREFIX='GRCh37' # GRCh37 | GRCh38 
        ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
        source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
        deseq.RData <- "/home/ssg29/results/RNA-Seq/Boy.Girl.ILM.GRCh37/DESeq2.1.18.1/BR/deseq.BR.RData"
        dl.pt.deseq[["FtCount"]][["ILM"]]<-prep.dt.deg(deseq.RData)
        cat("FtCount-ILM loaded\n")

        ##
        save(dl.pt.deseq, file=my.RData.file)
        cat("dl.pt.deseq saved\n")
        ##
    }
}

##########################
## Load Read Count Data ##
##########################
if(TRUE){
    my.RData.file="/home/ssg29/results/RNA-Seq/Placentome/RData/POPS.ILM.read.cnt.RData"
    if(file.exists(my.RData.file)){
        load(my.RData.file) # load 'dl.cnt
    }else{
        dl.cnt<-list()
        ## Illumina Read Count
        my.cmd<-"awk '/r_1/{split($1,a,\"/\"); split(a[7],b,\".\"); print a[6],b[2],$2}' /home/ssg29/results/SLX-ILM-Plasma2017.Homo_sapiens.v1/SLX-ILM-Plasma2017.fq.read.base.cnt.txt"
        dummy<-system(my.cmd,intern=T)
        dt.foo<-data.table(t(simplify2array(strsplit(dummy, " "))))  # isa 'matrix' 
        setnames(dt.foo, c("Library","BarCode","ReadCnt"))
        dt.foo$ReadCnt<-as.numeric(dt.foo$ReadCnt)

        dl.cnt[["ILM"]]<-merge(dt.foo,
                        rbind(
                                fread("/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/Meta/meta.ALL.Salmon.csv")[,.(Library,BarCode,SampleName,CRN,Sex,GA,`Source`="ILM",Tissue="Plasma")],
                                fread("/home/ssg29/results/RNA-Seq/Boy.Girl.ILM.GRCh38/Meta/meta.ALL.Salmon.csv")[,.(Library,BarCode,SampleName,CRN,Sex,`GA`=NA,`Source`="ILM",Tissue="Placenta")]
                        )
            )

        # POPS Placenta and Plasma total-RNA-Seq Read Count
        FILES=c("/home/ssg29/results/SLX-9168.Homo_sapiens.SE125.v2/SLX-9168.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-9169.Homo_sapiens.SE125.v2/SLX-9169.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-10281.Homo_sapiens.v1/SLX-10281.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-10402.Homo_sapiens.v1/SLX-10402.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-9792.Homo_sapiens.v1/SLX-9792.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-10284.Homo_sapiens.v1/SLX-10284.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-10285.Homo_sapiens.v1/SLX-10285.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-10283.Homo_sapiens.v1/SLX-10283.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-10287.Homo_sapiens.v1/SLX-10287.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-9342.Homo_sapiens.PE75.v2/SLX-9342.fq.read.base.cnt.txt",
                "/home/ssg29/results/SLX-9345.Homo_sapiens.PE75.v2/SLX-9345.fq.read.base.cnt.txt"
                )
        dl.readCount<-list()
        dl.readCount<-lapply(FILES, function(i){
                                        if(file.exists(i)){
                                                my.cmd<-paste0("awk '/r_1/{split($1,a,\"/\"); split(a[7],b,\".\"); print a[6],b[2],$2}' ",i)
                                                dummy<-system(my.cmd,intern=T)
                                                dt.foo<-data.table(t(simplify2array(strsplit(dummy, " ")))) # isa 'matrix' 
                                                setnames(dt.foo, c("Library","BarCode","ReadCnt"))
                                                dt.foo$ReadCnt<-as.numeric(dt.foo$ReadCnt)
                                                dt.foo<-dt.foo[,.(`ReadCnt`=sum(ReadCnt)),"Library,BarCode"]
                                                #dl.readCount[[cnt]]=dt.foo
                                                #cnt<-cnt+1
                                        }else{
                                                stop(paste(i,"Not Found\n"))
                                        }
                                    }
                        )
            
        dl.cnt[["POPS"]]<-merge(rbindlist(dl.readCount),
                        # meta info of POPS placenta and plasma tissue
                        rbind(
                            fread("/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh38/Meta/meta.ALL.csv")[,.(Library,BarCode,SampleName,CRN,Sex,`GA`=NA,`Source`="POPS",Tissue="Placenta")],
                            fread("/home/ssg29/results/RNA-Seq/Plasma.2017/Meta/meta.ALL.csv")[,.(Library,BarCode,SampleName,CRN,Sex,GA,`Source`="POPS",Tissue="Plasma")]
                        )
                    )
        save(dl.cnt, file=my.RData.file)
    }
}

##################################
## Comparison between FG and ILM #
## based on Salmon-quant         #
##################################
if(FALSE){
    ##############
    # No. of DEG #
    ##############
    #1-a. Placenta tissue without abudance threshold
    dl.pt.deseq$Salmon$FG[padj<0.05,.N,(`Female>Male`=log2FoldChange>0)]
    dl.pt.deseq$Salmon$FG[padj<0.05,.N,.(`Female<Male`=log2FoldChange<0,`chromosome_name==Y`=chromosome_name=="Y")]
    dl.pt.deseq$Salmon$FG[padj<0.05,.N,.(`Female>=Male`=log2FoldChange>=0,`chromosome_name==X`=chromosome_name=="X")]

    dl.pt.deseq$Salmon$ILM[padj<0.05,.N,(`Female>Male`=log2FoldChange>0)]
    dl.pt.deseq$Salmon$ILM[padj<0.05,.N,.(`Female<Male`=log2FoldChange<0,`chromosome_name==Y`=chromosome_name=="Y")]
    dl.pt.deseq$Salmon$ILM[padj<0.05,.N,.(`Female>=Male`=log2FoldChange>=0,`chromosome_name==X`=chromosome_name=="X")]

    merge(dl.pt.deseq$Salmon$FG[padj<0.05 & log2FoldChange>0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], dl.pt.deseq$Salmon$ILM[padj<0.05 & log2FoldChange>0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], by=c("hgnc_symbol","chromosome_name"))[,.N,chromosome_name]
    merge(dl.pt.deseq$Salmon$FG[padj<0.05 & log2FoldChange<0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], dl.pt.deseq$Salmon$ILM[padj<0.05 & log2FoldChange<0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], by=c("hgnc_symbol","chromosome_name"))[,.N,chromosome_name]

    #1-b. Placenta tissue with abudance threshold
    dl.pt.deseq$Salmon$FG[padj<0.05 & meanFpm>1 & abs(log2FoldChange)>1 & meanFpm>=1,.N,(`Female>Male`=log2FoldChange>0)]
    dl.pt.deseq$Salmon$FG[padj<0.05,.N,.(`Female<Male`=log2FoldChange<0,`chromosome_name==Y`=chromosome_name=="Y")]
    dl.pt.deseq$Salmon$FG[padj<0.05,.N,.(`Female>=Male`=log2FoldChange>=0,`chromosome_name==X`=chromosome_name=="X")]

    dl.pt.deseq$Salmon$FG[chromosome_name=="Y" & padj<0.05 & abs(log2FoldChange)>1 & meanFpm>=1,.(hgnc_symbol,baseMean,log2FoldChange,padj,meanFpm,chromosome_name)][order(-meanFpm)]
    dl.pt.deseq$Salmon$ILM[chromosome_name=="Y" & padj<0.05 & abs(log2FoldChange)>1 & meanFpm>=1,.(hgnc_symbol,baseMean,log2FoldChange,padj,meanFpm,chromosome_name)][order(-meanFpm)]

    dl.pt.deseq$Salmon$FG[hgnc_symbol=="NLGN4Y"]
    dl.pt.deseq$Salmon$ILM[hgnc_symbol=="NLGN4Y"]
    dl.pt.deseq$Salmon$ILM[hgnc_symbol=="DDX3Y"]
    dl.pt.deseq$Salmon$ILM[hgnc_symbol=="RPS4Y1"]

    merge(dl.pt.deseq$Salmon$FG[padj<0.05 &  meanFpm>1 & abs(log2FoldChange)>1 & meanFpm>=1,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], dl.pt.deseq$Salmon$ILM[padj<0.05 &  meanFpm>1 & abs(log2FoldChange)>1 & meanFpm>=1,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], by=c("hgnc_symbol","chromosome_name"))

    #2. Plasma
    dl.pl.deseq$Salmon$FG[padj<0.05,.N,(`Female>Male`=log2FoldChange>0)]
    dl.pl.deseq$Salmon$FG[padj<0.05,.N,.(`Female<Male`=log2FoldChange<0,`chromosome_name==Y`=chromosome_name=="Y")]
    dl.pl.deseq$Salmon$FG[padj<0.05,.N,.(`Female>=Male`=log2FoldChange>=0,`chromosome_name==X`=chromosome_name=="X")]

    dl.pl.deseq$Salmon$ILM[padj<0.05,.N,(`Female>Male`=log2FoldChange>0)]
    dl.pl.deseq$Salmon$ILM[padj<0.05,.N,.(`Female<Male`=log2FoldChange<0,`chromosome_name==Y`=chromosome_name=="Y")]
    dl.pl.deseq$Salmon$ILM[padj<0.05,.N,.(`Female>=Male`=log2FoldChange>=0,`chromosome_name==X`=chromosome_name=="X")]

    dl.pl.deseq$Salmon$FG[chromosome_name=="Y" & padj<0.05 & abs(log2FoldChange)>1 & meanFpm>=1,.(hgnc_symbol,baseMean,log2FoldChange,padj,meanFpm,chromosome_name)][order(-meanFpm)]
    dl.pl.deseq$Salmon$ILM[chromosome_name=="Y" & padj<0.05 & abs(log2FoldChange)>1 & meanFpm>=1,.(hgnc_symbol,baseMean,log2FoldChange,padj,meanFpm,chromosome_name)][order(-meanFpm)]

    dl.pl.deseq$Salmon$ILM[hgnc_symbol=="KDM5D"]

    merge(dl.pl.deseq$Salmon$FG[padj<0.05 & log2FoldChange>0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], dl.pl.deseq$Salmon$ILM[padj<0.05 & log2FoldChange>0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], by=c("hgnc_symbol","chromosome_name"))[,.N,chromosome_name]
    merge(dl.pl.deseq$Salmon$FG[padj<0.05 & log2FoldChange<0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], dl.pl.deseq$Salmon$ILM[padj<0.05 & log2FoldChange<0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], by=c("hgnc_symbol","chromosome_name"))[,.N,chromosome_name]

    #3. DEG of Plasma (Salmon) by GA (P<0.05)
    for(j in paste("GA",c(12, 20, 28, 36),sep=".")){
        # FG
        dl.pl.ga.deseq[[j]]$FG[padj<0.05,.N,(`Female>Male`=log2FoldChange>0)]
        dl.pl.ga.deseq[[j]]$FG[padj<0.05,.N,.(`Female<Male`=log2FoldChange<0,`chromosome_name==Y`=chromosome_name=="Y")]
        dl.pl.ga.deseq[[j]]$FG[padj<0.05,.N,.(`Female>=Male`=log2FoldChange>=0,`chromosome_name==X`=chromosome_name=="X")]
        # ILM
        dl.pl.ga.deseq[[j]]$ILM[padj<0.05,.N,(`Female>Male`=log2FoldChange>0)]
        dl.pl.ga.deseq[[j]]$ILM[padj<0.05,.N,.(`Female<Male`=log2FoldChange<0,`chromosome_name==Y`=chromosome_name=="Y")]
        dl.pl.ga.deseq[[j]]$ILM[padj<0.05,.N,.(`Female>=Male`=log2FoldChange>=0,`chromosome_name==X`=chromosome_name=="X")]

        # common F>M
        merge(dl.pl.ga.deseq[[j]]$FG[padj<0.05 & log2FoldChange>0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], 
              dl.pl.ga.deseq[[j]]$ILM[padj<0.05 & log2FoldChange>0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], 
              by=c("hgnc_symbol","chromosome_name"))[,.N,chromosome_name]
        # common F<M
        merge(dl.pl.ga.deseq[[j]]$FG[padj<0.05 & log2FoldChange<0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], 
              dl.pl.ga.deseq[[j]]$ILM[padj<0.05 & log2FoldChange<0,.(hgnc_symbol,chromosome_name,log2FoldChange,padj)], 
              by=c("hgnc_symbol","chromosome_name"))[,.N,chromosome_name]
    }
}


get_plasma_by_gene<-function(my.target.genes){
        #2. top 4 plasma chrY genes by GA: PCDH11Y, DDX3Y, ZFY, RPS4Y1
        #my.target.genes=c('PCDH11Y','DDX3Y','ZFY','RPS4Y1')
        dl.dummy<-list(); cnt<-1
        for(i in c("FG", "ILM")){
            for(j in paste("GA",c(12, 20, 28, 36),sep=".")){
                dl.dummy[[cnt]]<-cbind(`Source`=i, `GA`=limma::strsplit2(j,"GA.")[1,2], dl.pl.ga.deseq[[j]][[i]][hgnc_symbol %in% my.target.genes])
                cnt<-cnt+1
            }
        }
        dt.dummy<-rbindlist(dl.dummy)
        dt.dummy$hgnc_symbol<-droplevels(dt.dummy$hgnc_symbol)
        dt.dummy$hgnc_symbol<-factor(dt.dummy$hgnc_symbol, levels=my.target.genes)
        return(dt.dummy)
}

get_plasma_by_pval<-function(my.pval){
    dl.dummy<-list(); cnt<-1
    for(i in c("FG", "ILM")){
        for(j in paste("GA",c(12, 20, 28, 36),sep=".")){
            if(nrow(dl.pl.ga.deseq[[j]][[i]][padj< my.pval])>=1){
                dl.dummy[[cnt]]<-cbind(`Source`=i, `GA`=limma::strsplit2(j,"GA.")[1,2], dl.pl.ga.deseq[[j]][[i]][padj< my.pval])
                cnt<-cnt+1
            }
        }
    }
    dt.dummy<-rbindlist(dl.dummy)
    dt.dummy$hgnc_symbol<-droplevels(dt.dummy$hgnc_symbol)
    return(dt.dummy)
}

