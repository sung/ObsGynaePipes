#!/usr/bin/Rscript --vanilla
#library(ggpubr)
library(scales) # for 'trans_breaks'
source("~/Pipelines/config/graphic.R")
TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
ENS_VER=82 # 75 (for igenome), 82 (for manual downlaod) 90 (for Cholestatis project)
source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
source("~/Pipelines/bin/R/Plasma-RNA/local.R") # load necessary R data files for DATA VIZ
# dl.pl.deseq : $Salmon$FG, $Salmon$ILM, $FtCount$FG, $FtCount$ILM
# dl.pl.ga.deseq : $GA12$FG, $GA12$ILM ... $GA36$FG, $GA36$ILM
# dl.pt.deseq : $Salmon$FG, $Salmon$ILM, $FtCount$FG, $FtCount$ILM
# dl.cnt : (readcount) $ILM, $POPS

#############
## DATA VIZ #
#############
if(TRUE){
    # prep data
    if(TRUE){
        # raw read-count data
        dt.cnt<-rbind(
                    # POPS Placenta and illumina Placenta
                    merge(dl.cnt[["ILM"]][grepl("^RNA",BarCode)], dl.cnt[["POPS"]][Tissue=="Placenta" & grepl("^BR",SampleName) & !SampleName %in% c("BR08","BR17")], by=c("CRN","GA","Sex","Tissue"), all=T), # n=22 (ILM=19 & POPS=22)
                    # POPS Plasma and illumina Plasma 
                    merge(dl.cnt[["ILM"]][grepl("^PL",BarCode)], dl.cnt[["POPS"]][Tissue=="Plasma"], by=c("CRN","GA","Sex","Tissue"), all=T) # n=87 (ILM=83, POPS=86)
                )

        #########################
        ## FPM of target genes ##
        #########################
        #1. top 11 plasma chrY genes regardless of GA
        my.target.chrY=c('PRKY','NLGN4Y','TXLNGY','EIF1AY','PCDH11Y','KDM5D','ZFY','UTY','USP9Y','RPS4Y1','DDX3Y')
        dt.top.chrY<-list()
        dt.top.chrY[["plasma"]]<-rbind(
                                    cbind(`Source`="FG",dl.pl.deseq[["Salmon"]][["FG"]][hgnc_symbol %in% my.target.chrY]),
                                    cbind(`Source`="ILM",dl.pl.deseq[["Salmon"]][["ILM"]][hgnc_symbol %in% my.target.chrY])
                                )
        dt.top.chrY[["plasma"]]$hgnc_symbol<-droplevels(dt.top.chrY[["plasma"]]$hgnc_symbol)
        dt.top.chrY[["plasma"]]$hgnc_symbol<-factor(dt.top.chrY[["plasma"]]$hgnc_symbol, levels=my.target.chrY)

        dt.top.chrY[["placenta"]]<-rbind(
                                    cbind(`Source`="FG",dl.pt.deseq[["Salmon"]][["FG"]][hgnc_symbol %in% my.target.chrY]),
                                    cbind(`Source`="ILM",dl.pt.deseq[["Salmon"]][["ILM"]][hgnc_symbol %in% my.target.chrY])
                                )
        dt.top.chrY[["placenta"]]$hgnc_symbol<-droplevels(dt.top.chrY[["placenta"]]$hgnc_symbol)
        dt.top.chrY[["placenta"]]$hgnc_symbol<-factor(dt.top.chrY[["placenta"]]$hgnc_symbol, levels=my.target.chrY)

        #2. top 4 plasma chrY genes by GA: PCDH11Y, DDX3Y, ZFY, RPS4Y1
        my.target.genes=c('PCDH11Y','DDX3Y','ZFY','RPS4Y1')
        dt.top.chrY.GA<-get_plasma_by_gene(my.target.genes)
        #dl.dummy<-list(); cnt<-1
        #for(i in c("FG", "ILM")){
        #    for(j in paste("GA",c(12, 20, 28, 36),sep=".")){
        #        dl.dummy[[cnt]]<-cbind(`Source`=i, `GA`=limma::strsplit2(j,"GA.")[1,2], dl.pl.ga.deseq[[j]][[i]][hgnc_symbol %in% my.target.genes])
        #        cnt<-cnt+1
        #    }
        #}
        #dt.top.chrY.GA<-rbindlist(dl.dummy)
        #dt.top.chrY.GA$hgnc_symbol<-droplevels(dt.top.chrY.GA$hgnc_symbol)
        #dt.top.chrY.GA$hgnc_symbol<-factor(dt.top.chrY.GA$hgnc_symbol, levels=my.target.genes)

        #3. top 10 plasma chrX genes regardless of GA
        my.target.chrX=c('TSIX','HDHD1','EIF1AX','ZFX','KDM6A','SMC1A','KDM5C','DDX3X','STS','XIST')
        dt.top.chrX<-rbind(
                        cbind(`Source`="FG",dl.pl.deseq[["Salmon"]][["FG"]][hgnc_symbol %in% my.target.chrX]),
                        cbind(`Source`="ILM",dl.pl.deseq[["Salmon"]][["ILM"]][hgnc_symbol %in% my.target.chrX])
                        )
        dt.top.chrX$hgnc_symbol<-droplevels(dt.top.chrX$hgnc_symbol)
        dt.top.chrX$hgnc_symbol<-factor(dt.top.chrX$hgnc_symbol, levels=my.target.chrX)

        #4. FG-ILM Plasma FPM
        dt.pl.fg.ilm<-merge(dl.pl.deseq[["Salmon"]][["FG"]][,.(hgnc_symbol,chromosome_name,log2FoldChange,padj,meanFpm)],
                            dl.pl.deseq[["Salmon"]][["ILM"]][,.(hgnc_symbol,chromosome_name,log2FoldChange,padj,meanFpm)],
                            by=c("hgnc_symbol","chromosome_name"))
    }

    if(TRUE){
        my.file.name<-"/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/Cluster/top.chrY.genes.plasma.FG.vs.ILM"
        pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title="RNA-Seq In-house vs. ILM") # A4 size

        #
        p0<-ggplot(dt.cnt, aes(ReadCnt.x/10^6, ReadCnt.y/10^6, color=Tissue)) +
            geom_point(size=5,alpha=.7) +  
            coord_cartesian(ylim=c(0,170),xlim=c(0,170)) +
            geom_abline(slope=1,linetype="dotted")+
            xlab("illumina (PE51)") + ylab("In-house (PE75|SE125)") + ggtitle("No. of Raw Read Count (per Million)") +
            theme_Publication() +
            ggpubr::color_palette("jco")
        print(p0)
        # by ggpubr
        #ggscatter(dt.cnt[!is.na(Source)], x="ReadCnt.x", y="ReadCnt.y", color="Source", size=5, alpha=.7, palette="jco")
        if(FALSE){
            library(cowplot)
            # Marginal boxplot along x axis
            p0.ybox <- ggplot(dt.cnt, aes(Tissue,ReadCnt.y/10^6, color=Tissue)) +
                        geom_boxplot(alpha=0.5, outlier.shape=NA, width=.3) +  
                        theme_Publication() +
                        ggpubr::color_palette("jco")
            p0.y <- axis_canvas(p0, axis = "y")

            # Marginal boxplot along y axis
            p0.ybox<- axis_canvas(p0, axis = "y", coord_flip=TRUE) +
                        ggplot(dt.cnt, aes(Tissue,ReadCnt.y/10^6, color=Tissue)) +
                        geom_boxplot(alpha=0.5, outlier.shape=NA, width=.3) +  
                        theme_Publication() +
                        ggpubr::color_palette("jco") +
                        coord_flip()
            p0_insert_xaxis <- grob(p0, p0.xbox, grid::unit(.2, "null"), position = "top")
            p0_insert_yaxis <- grob(p0_insert_xaxis, p0.ybox, grid::unit(.2, "null"), position = "right")
            ggdraw(p0_insert_yaxis)
        }

        p1<-ggplot(dt.top.chrY[["placenta"]], aes(hgnc_symbol, meanFpm)) +
            geom_bar(stat="identity", position="dodge") +
            xlab("Gene") + ylab("FPM (Fragment Per Million)") +
            ggtitle("Placenta Tissue") +
            facet_wrap(~Source,ncol=2) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p1)

        dt.foo<-melt(dt.top.chrY[["placenta"]][,.(Source,hgnc_symbol,`Female`=meanFpm.case,`Male`=meanFpm.control)],id=c("Source","hgnc_symbol"),measure=c("Female","Male"),variable.name="Sex",value.name="FPM")  # col to row 
        p1.a<-ggplot(dt.foo, aes(hgnc_symbol, FPM, fill=Sex)) +
            geom_bar(stat="identity", position="dodge") +
            xlab("Gene") + ylab("FPM (Fragment Per Million)") +
            ggtitle("Placenta Tissue") +
            facet_wrap(~Source,ncol=2) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p1.a)

        p1.b<-ggplot(dt.top.chrY[["plasma"]], aes(hgnc_symbol, meanFpm)) +
            geom_bar(stat="identity", position="dodge") +
            xlab("Gene") + ylab("FPM (Fragment Per Million)") +
            ggtitle("Plasma") +
            facet_wrap(~Source,ncol=2) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p1.b)

        dt.foo<-melt(dt.top.chrY[["plasma"]][,.(Source,hgnc_symbol,`Female`=meanFpm.case,`Male`=meanFpm.control)],id=c("Source","hgnc_symbol"),measure=c("Female","Male"),variable.name="Sex",value.name="FPM")  # col to row 
        p1.c<-ggplot(dt.foo, aes(hgnc_symbol, FPM, fill=Sex)) +
            geom_bar(stat="identity", position="dodge") +
            xlab("Gene") + ylab("FPM (Fragment Per Million)") +
            ggtitle("Plasma") +
            facet_wrap(~Source,ncol=2) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p1.c)

        p2<-ggplot(dt.top.chrY.GA, aes(hgnc_symbol, meanFpm, fill=GA)) +
            geom_bar(stat="identity", position="dodge") +
            xlab("Gene") + ylab("FPM (Fragment Per Million)") +
            facet_wrap(~Source,ncol=2) +
            theme_Publication()
        print(p2)

        p3<-ggplot(dt.top.chrY.GA, aes(hgnc_symbol, baseMean,fill=GA)) +
            geom_bar(stat="identity", position="dodge") +
            xlab("Gene") + ylab("Avg. Count") +
            facet_wrap(~Source,ncol=2) +
            theme_Publication()
        print(p3)

        p4<-ggplot(dt.top.chrX, aes(hgnc_symbol, meanFpm)) +
            geom_bar(stat="identity", position="dodge") +
            xlab("Gene") + ylab("FPM (Fragment Per Million)") +
            facet_wrap(~Source,ncol=2) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p4)

        ##############################
        # R^2 for all genes of FPM>1 #
        ##############################
        # FPM>1
        if(TRUE){
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1], aes(meanFpm.x, meanFpm.y,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(data=dt.pl.fg.ilm[meanFpm.x>40000 | meanFpm.y>50000],check_overlap=TRUE,vjust=-.9,size=4) + 
            geom_abline(slope=1,linetype="dotted")+
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            ggtitle(paste0("All genes of FPM>1 (R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1, .(meanFpm.x, meanFpm.y)]))$r.squared,2),", n=",nrow(dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1]),")")) +
            theme_Publication()
        print(p5)
        }
        # R^2 for all genes of FPM>1 (zoom-in)
        if(TRUE){
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1], aes(meanFpm.x, meanFpm.y,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(data=dt.pl.fg.ilm[meanFpm.x>5000 | meanFpm.y>7500],check_overlap=T,vjust=-1,size=4) + 
            coord_cartesian(ylim=c(0,25000),xlim=c(0,25000)) +
            geom_abline(slope=1,linetype="dotted")+
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            ggtitle(paste0("R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.x<25000 & meanFpm.y>1 & meanFpm.y<25000, .(meanFpm.x, meanFpm.y)]))$r.squared,2))) +
            ggtitle(paste0("All genes of FPM<25K (R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1, .(meanFpm.x, meanFpm.y)]))$r.squared,2),", n=",nrow(dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1 & meanFpm.x<25000 & meanFpm.y<25000]),")")) +
            theme_Publication()
        print(p5)
        }
        # R^2 for all genes of FPM>1 (zoom-in log-scale)
        if(TRUE){
        ## log2(FPM)
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1], aes(meanFpm.x, meanFpm.y)) +
            geom_point(size=5,alpha=.6) +
            #coord_trans(y="log2", x="log2") +
            scale_x_continuous(trans = log2_trans(),
                               breaks = trans_breaks("log2", function(x) 2^x),
                               labels = trans_format("log2", math_format(2^.x))) +
            scale_y_continuous(trans = log2_trans(),
                               breaks = trans_breaks("log2", function(x) 2^x),
                               labels = trans_format("log2", math_format(2^.x))) +
            geom_abline(slope=1,linetype="dotted")+
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            ggtitle(paste0("All genes of FPM>1 (R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1, .(meanFpm.x, meanFpm.y)]))$r.squared,2),", n=",nrow(dt.pl.fg.ilm[meanFpm.x>1 & meanFpm.y>1]),")")) +
            theme_Publication()
        print(p5)
        }

        if(TRUE){
        # chrY
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y >0 & chromosome_name=="Y"], aes(meanFpm.x, meanFpm.y,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(data=dt.pl.fg.ilm[meanFpm.x+meanFpm.y>1000 & chromosome_name=="Y"],hjust=1.1,size=4) + 
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            geom_abline(slope=1,linetype="dotted")+
            ggtitle(paste0("chrY genes of FPM>0 (n=",nrow(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y >0 & chromosome_name=="Y"]),")")) +
            theme_Publication()
        print(p5)
        }

        if(TRUE){
        # chrY zoom-in 1
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<1000 & chromosome_name=="Y"], aes(meanFpm.x, meanFpm.y,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & chromosome_name=="Y" & meanFpm.x<1000 & (meanFpm.x>.5 | meanFpm.y>5)],check_overlap=TRUE,vjust=-1,size=4) + 
            geom_abline(slope=1,linetype="dotted")+
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            ggtitle(paste0("chrY genes of FPM>0 & FPM<1K (R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<1000 & chromosome_name=="Y", .(meanFpm.x, meanFpm.y)]))$r.squared,2),", n=",nrow(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y >0 & meanFpm.x<1000 & chromosome_name=="Y"]),")")) +
            theme_Publication()
        print(p5)
        }

        if(TRUE){
        # chrY zoom-in 2
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<1000 & chromosome_name=="Y"], aes(meanFpm.x, meanFpm.y,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & chromosome_name=="Y" & meanFpm.x<1000 & (meanFpm.x>.1 | meanFpm.y>5)],check_overlap=TRUE,vjust=-1,size=4) + 
            coord_cartesian(xlim=c(0,0.5)) +
            geom_abline(slope=1,linetype="dotted")+
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            ggtitle(paste0("chrY genes of FPM>0 & FPM(In-house)<0.5 (n=",nrow(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y >0 & meanFpm.x<0.5 & chromosome_name=="Y"]),")")) +
            ggtitle(paste0("chrY genes of FPM>0 & FPM(In-house)<0.5 (R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<0.5 & chromosome_name=="Y", .(meanFpm.x, meanFpm.y)]))$r.squared,2),", n=",nrow(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y >0 & meanFpm.x<0.5 & chromosome_name=="Y"]),")")) +
            theme_Publication()
        print(p5)
        }

        if(TRUE){
        # chrY zoom-in 1 (log2 scale)
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<1000 & chromosome_name=="Y"], aes(meanFpm.x, meanFpm.y,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & chromosome_name=="Y" & meanFpm.x<1000 & (meanFpm.x>.1 | meanFpm.y>5)],check_overlap=TRUE,vjust=-1,size=4) + 
            coord_trans(y="log2", x="log2") +
            geom_abline(slope=1,linetype="dotted")+
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            ggtitle(paste0("chrY genes of FPM>0 & FPM<1K (R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<1000 & chromosome_name=="Y", .(meanFpm.x, meanFpm.y)]))$r.squared,2),", n=",nrow(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y >0 & meanFpm.x<1000 & chromosome_name=="Y"]),")")) +
            theme_Publication()
        print(p5)
        }

        if(TRUE){
        # chrY zoom-in 1 (log2 scale)
        p5<-ggplot(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<1000 & chromosome_name=="Y"], aes(meanFpm.x, meanFpm.y,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & chromosome_name=="Y" & meanFpm.x<1000 & (meanFpm.x>.1 | meanFpm.y>5)],check_overlap=TRUE,vjust=-1,size=4) + 
            scale_x_continuous(trans = log2_trans(),
                               breaks = trans_breaks("log2", function(x) 2^x),
                               labels = trans_format("log2", math_format(2^.x))) +
            scale_y_continuous(trans = log2_trans(),
                               breaks = trans_breaks("log2", function(x) 2^x),
                               labels = trans_format("log2", math_format(2^.x))) +
            geom_abline(slope=1,linetype="dotted")+
            xlab("FPM (In-House)") + ylab("FPM (illumina)") +
            ggtitle(paste0("chrY genes of FPM>0 & FPM<1K (R^2=",round(summary(lm(data=dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y>0 & meanFpm.x<1000 & chromosome_name=="Y", .(meanFpm.x, meanFpm.y)]))$r.squared,2),", n=",nrow(dt.pl.fg.ilm[meanFpm.x>0 & meanFpm.y >0 & meanFpm.x<1000 & chromosome_name=="Y"]),")")) +
            theme_Publication()
        print(p5)
        }
        dev.off()
    }

    ###########################
    # Placenta Specific Genes #
    ###########################
    if(TRUE){
        #5. 106 placenta specific genes (see ~/Pipelines/bin/R/Placentome/placenta.specific.expression.R)
        pt.specific.gene<-fread("~/results/RNA-Seq/Placentome/CSV/placenta.specific.protein.coding.genes.FPKM.Tau.0.99.PT.FPKM.0.1.ratio.100.csv")[,ensembl_gene_id]
        dt.pl.top.pt.gene<-dcast.data.table(
                                rbindlist(lapply(names(dl.pl.deseq$Salmon), function(i) dl.pl.deseq$Salmon[[i]][ensembl_gene_id %in% pt.specific.gene,.(`Type`=i,ensembl_gene_id,meanFpkm)])) ,
                                ensembl_gene_id ~ Type, 
                                value.var="meanFpkm") # long (row-based) to wide (column-based)
        my.file.name<-"/home/ssg29/results/RNA-Seq/Plasma.ILM.2017/Cluster/placenta.specific.genes.plasma.FPKM.FG.vs.ILM"
        pdf(file=paste(my.file.name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'), width=11.7, height=8.3, title="RNA-Seq In-house vs. ILM") # A4 size

        if(TRUE){
        print(dt.pl.top.pt.gene[,.N,.(`In-house < illumina`=FG<ILM)])
        p6<-ggplot(merge(dt.pl.top.pt.gene,dt.ensg), aes(FG,ILM,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(check_overlap=TRUE,hjust=1.2,size=4) + 
            scale_x_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) +
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) +
            geom_abline(slope=1,linetype="dotted")+
            ggtitle(paste(nrow(dt.pl.top.pt.gene), "Placenta Specific Genes")) +
            xlab("FPKM (In-House)") + ylab("FPKM (illumina)") +
            theme_Publication()
        print(p6)
        }
        if(TRUE){
        my.cutoff<-0.1
        print(dt.pl.top.pt.gene[,.N,.(`In-house < illumina`=FG<ILM)])
        dt.foo<-merge(dt.pl.top.pt.gene,dt.ensg)
        dt.foo[hgnc_symbol=="CSH1"]
        p6<-ggplot(dt.foo, aes(FG,ILM,label=hgnc_symbol)) +
            geom_point(data=dt.foo[ILM>my.cutoff],size=5,alpha=.6) +
            geom_point(data=dt.foo[hgnc_symbol=="PSG7"],size=5,alpha=.6,color="blue") +
            geom_point(data=dt.foo[ILM<=my.cutoff],size=1,alpha=.6) +
            geom_text(data=dt.foo[ILM>my.cutoff],check_overlap=TRUE,hjust=1.2,size=4) + 
            scale_x_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) +
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) +
            geom_abline(slope=1,linetype="dotted")+
            geom_hline(yintercept=my.cutoff, linetype="dashed", color = "red", size=1) +
            ggtitle(paste(nrow(dt.pl.top.pt.gene[ILM>my.cutoff]), "Placenta Specific Genes (FPKM (illumina) > ",my.cutoff,")")) +
            xlab("FPKM (In-House)") + ylab("FPKM (illumina)") +
            theme_Publication()
        print(p6)
        }
        if(TRUE){
        print(dt.pl.top.pt.gene[,.N,.(`In-house < illumina`=FG<ILM)])
        dt.pl.top.pt.gene[,.N,.(`In-house < illumina`=FG<ILM)]
        p6<-ggplot(merge(dt.pl.top.pt.gene,dt.ensg), aes(FG,ILM,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
            coord_cartesian(xlim=c(0,70),ylim=c(0,70)) +
            geom_abline(slope=1,linetype="dotted")+
            ggtitle(paste(nrow(dt.pl.top.pt.gene), "Placenta Specific Genes")) +
            xlab("FPKM (In-House)") + ylab("FPKM (illumina)") +
            theme_Publication()
        print(p6)
        }
        if(FALSE){
        print(dt.pl.top.pt.gene[FG<10 & ILM<10,.N,.(`In-house < illumina`=FG<ILM)])
        p6<-ggplot(merge(dt.pl.top.pt.gene,dt.ensg), aes(FG,ILM,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
            coord_cartesian(xlim=c(0,10),ylim=c(0,10)) +
            geom_abline(slope=1,linetype="dotted")+
            ggtitle(paste(nrow(dt.pl.top.pt.gene[FG<10 & ILM<10]), "Placenta Specific Genes (FPKM<10)")) +
            xlab("FPKM (In-House)") + ylab("FPKM (illumina)") +
            theme_Publication()
        print(p6)
        }
        if(FALSE){
        print(dt.pl.top.pt.gene[FG<2 & ILM<2,.N,.(`In-house < illumina`=FG<ILM)])
        p6<-ggplot(merge(dt.pl.top.pt.gene,dt.ensg), aes(FG,ILM,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
            #coord_trans(y="log2", x="log2") +
            coord_cartesian(xlim=c(0,2),ylim=c(0,2)) +
            geom_abline(slope=1,linetype="dotted")+
            ggtitle(paste(nrow(dt.pl.top.pt.gene[FG<2 & ILM<2]), "Placenta Specific Genes (FPKM<2)")) +
            xlab("FPKM (In-House)") + ylab("FPKM (illumina)") +
            theme_Publication()
        print(p6)
        }
        if(FALSE){
        print(dt.pl.top.pt.gene[FG<.2 & ILM<.2,.N,.(`In-house < illumina`=FG<ILM)])
        p6<-ggplot(merge(dt.pl.top.pt.gene,dt.ensg), aes(FG,ILM,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
            coord_cartesian(xlim=c(0,.2),ylim=c(0,.2)) +
            geom_abline(slope=1,linetype="dotted")+
            ggtitle(paste(nrow(dt.pl.top.pt.gene[FG<.2 & ILM<.2]), "Placenta Specific Genes (FPKM<.2)")) +
            xlab("FPKM (In-House)") + ylab("FPKM (illumina)") +
            theme_Publication()
        print(p6)
        }
        if(FALSE){
        print(dt.pl.top.pt.gene[FG<.02 & ILM<.02,.N,.(`In-house < illumina`=FG<ILM)])
        p6<-ggplot(merge(dt.pl.top.pt.gene,dt.ensg), aes(FG,ILM,label=hgnc_symbol)) +
            geom_point(size=5,alpha=.6) +
            geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
            coord_cartesian(xlim=c(0,.02),ylim=c(0,.02)) +
            geom_abline(slope=1,linetype="dotted")+
            ggtitle(paste(nrow(dt.pl.top.pt.gene[FG<0.02 & ILM<0.02]), "Placenta Specific Genes (FPKM<.02)")) +
            xlab("FPKM (In-House)") + ylab("FPKM (illumina)") +
            theme_Publication()
        print(p6)
        }
        if(FALSE){
        dt.pl.top.pt.gene[,`Type`:=ifelse(FG<ILM,"FPKM(in-house) < FPKM(illumina)","FPKM(in-house) >= FPKM(illumina)")]
        p7<-ggplot(melt.data.table(dt.pl.top.pt.gene, id.vars=c("ensembl_gene_id","Type"), variable.name="Source", value.name="FPKM"),aes(log2(FPKM)))  + 
            geom_density(aes(fill=Source),alpha=0.5) +
            facet_grid(.~Type) +
            theme_Publication()
        print(p7)
        }
    }
    dev.off()
    ## end of pt-specific

}

if(FALSE){
    my.target.genes=c("SMS")
    dt.top.GA<-get_plasma_by_gene(my.target.genes)
    dt.top.GA<-get_plasma_by_pval(0.05)
    dt.top.GA[,.N,"Source,GA"]
}
