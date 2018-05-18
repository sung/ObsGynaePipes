TR_PREFIX='GRCh38' # GRCh37 | GRCh38 
ENS_VER=90 # GTEx and Placenta RNA-Seq quantification was made GRCh37.82, but just to borrow recent gene_biotype
source("~/Pipelines/config/Annotation.R") # load dt.ensg gr.ensg
source("~/Pipelines/bin/R/Placentome/local.R") # load fTau()
source("~/Pipelines/config/graphic.R") # for theme_Publication()
source("~/Pipelines/bin/R/Plasma-RNA/local.R") # load 'dl.pl.deseq'
load("~/results/RNA-Seq/GTEx/DESeq2.1.18.1/GTEx.PT.DESeq2.RData")   # loads dl.gtex.deseq & dl.Fpkm
                                                                    # GRCh37.82
                                                                    # placenta from healthy babies 
                                                                    # see bin/R/GTEx/local.R
#########################################
# Placenta & Plasma vs 20 GTEx-tissue ###
#########################################
if(TRUE){
    dt.gtex.deseq<-rbind(
            rbindlist(dl.gtex.deseq),  # see above (Plasma from Boy.Girl.FG.JD.GRCh37 AGA samples)
			#/home/ssg29/results/RNA-Seq/Boy.Girl.FG.JD.GRCh37/DESeq2.1.18.1/AGA/filtered_toptags_deseq.all.gene.AGA.csv.gz
            #cbind(dl.pl.deseq[["Salmon"]][["ILM"]], `Tissue`="Plasma")[,-c("meanFpm","meanFpm.case","meanFpm.control")] # GRCh38.82
            cbind(dl.pl.deseq[["FtCount"]][["ILM"]], `Tissue`="Plasma")[,-c("meanFpm","meanFpm.case","meanFpm.control")] # GRCh37.82
            )
    dt.gtex.deseq[,.N,Tissue]

    dt.fpkm<-dcast.data.table(dt.gtex.deseq, ensembl_gene_id~Tissue, value.var="meanFpkm") # long (row-based) to wide (column-based)
    dt.fpkm<-merge(dt.fpkm,dt.gtex.deseq[!Tissue %in% c("Placenta","Plasma"),.(`meanFpkmGTEx`=mean(meanFpkm)),ensembl_gene_id])

    # placenta-specific
    m.fpkm1<-as.matrix(dt.fpkm[,-c("ensembl_gene_id","Plasma","meanFpkmGTEx")]) # a matrix of FPKM
    rownames(m.fpkm1)<-dt.fpkm$ensembl_gene_id; head(m.fpkm1)
    Tau.placenta<-apply(m.fpkm1, 1, fTau)

    # plasma-specific
    m.fpkm2<-as.matrix(dt.fpkm[,-c("ensembl_gene_id","Placenta","meanFpkmGTEx")]) # a matrix of FPKM
    rownames(m.fpkm2)<-dt.fpkm$ensembl_gene_id; head(m.fpkm2)
    Tau.plasma<-apply(m.fpkm2, 1, fTau)

    dt.gtex.fpkm<-cbind(dt.fpkm,Tau.placenta,Tau.plasma)

    my.Tau=0.99; my.pt.fpkm=0.1; my.pl.fpkm=0.1; my.pt.ratio=100; my.pl.ratio=100
    # placenta-specific (any) genes (n=1504)
    dt.pt.specific.any<-merge(dt.gtex.fpkm[Tau.placenta>my.Tau & Placenta>my.pt.fpkm & Placenta/meanFpkmGTEx>=my.pt.ratio], dt.ensg[,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Placenta)]
    # placenta-specific protein coding genes (n=106)
    dt.pt.specific<-merge(dt.gtex.fpkm[Tau.placenta>my.Tau & Placenta>my.pt.fpkm & Placenta/meanFpkmGTEx>=my.pt.ratio], dt.ensg[gene_biotype=="protein_coding",.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Placenta)]

    # plasma-specific (any) genes (n=400)
    dt.pl.specific.any<-merge(dt.gtex.fpkm[Tau.plasma>my.Tau & Plasma>my.pl.fpkm & Plasma/meanFpkmGTEx>=my.pl.ratio], dt.ensg[,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Plasma)]
    # plasma-specific genes (n=123)
    dt.pl.specific<-merge(dt.gtex.fpkm[Tau.plasma>my.Tau & Plasma>my.pl.fpkm & Plasma/meanFpkmGTEx>=my.pl.ratio], dt.ensg[gene_biotype=="protein_coding",.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Plasma)]

    # common betwen the two any genes (n=37)
    dt.pt.pl.specific.any<-merge(dt.gtex.fpkm[Tau.plasma>my.Tau & Tau.placenta>my.Tau & Placenta>my.pt.fpkm & Plasma>my.pl.fpkm & Plasma/meanFpkmGTEx>=my.pt.ratio & Placenta/meanFpkmGTEx>=my.pl.ratio], dt.ensg[,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Plasma)]
    # common betwen the two protein coding (n=13)
    dt.pt.pl.specific<-merge(dt.gtex.fpkm[Tau.plasma>my.Tau & Tau.placenta>my.Tau & Placenta>my.pt.fpkm & Plasma>my.pl.fpkm & Plasma/meanFpkmGTEx>=my.pt.ratio & Placenta/meanFpkmGTEx>=my.pl.ratio], dt.ensg[gene_biotype=="protein_coding",.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Plasma)]

    dt.pt.specific[grepl("HIST",hgnc_symbol)]
    dt.pt.specific[grepl("PSG",hgnc_symbol)]

    my.file=paste0("~/results/RNA-Seq/Placentome/CSV/placenta.specific.protein.coding.genes.FPKM.Tau.",my.Tau,".PT.FPKM.",my.pt.fpkm,".ratio.",my.pt.ratio,".csv")
    write.csv(dt.pt.specific, file=my.file, row.names=F, quote=F)

    my.file=paste0("~/results/RNA-Seq/Placentome/CSV/plasma.specific.protein.coding.genes.FPKM.Tau.",my.Tau,".PL.FPKM.",my.pl.fpkm,".ratio.",my.pl.ratio,".csv")
    write.csv(dt.pl.specific, file=my.file, row.names=F, quote=F)

    dl.pt.deseq[["Salmon"]]
    , `Tissue`="Plasma")[,-c("meanFpm","meanFpm.case","meanFpm.control")] # GRCh38.82
    cbind(dl.pl.deseq[["FtCount"]][["ILM"]], `Tissue`="Plasma")[,-c("meanFpm","meanFpm.case","meanFpm.control")] # GRCh37.82

    ## X-Y plot of FPKM for the placena specific genes
    if(TRUE){
        summary(lm(data=dt.pt.specific[,.(Placenta,Plasma)]))$r.squared
        dt.pt.specific[,.(hgnc_symbol,Placenta,Plasma)]
        p1<-ggplot(dt.pt.specific, aes(Placenta,Plasma,label=hgnc_symbol)) +
                geom_point(size=5,alpha=.6) +
                geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
                #coord_trans(y="log2", x="log2") +
                #coord_cartesian(xlim=c(0,.2),ylim=c(0,.2)) +
                ggtitle(paste(nrow(dt.pt.specific), "Placenta Specific Protein-coding Genes")) +
                xlab("FPKM(Placenta)") + ylab("FPKM (Plasma - Illumina)") +
                theme_Publication()
        print(p1)
    }
    if(TRUE){
        library(scales) # for 'trans_breaks'
        summary(lm(data=dt.pt.specific[Plasma!=0,.(log2(Placenta),log2(Plasma))]))$r.squared
        p2<-ggplot(dt.pt.specific, aes(Placenta,Plasma,label=hgnc_symbol)) +
                geom_point(size=5,alpha=.6) +
                geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
                #coord_trans(y="log2", x="log2") +
                scale_x_continuous(trans = log2_trans(),
                                breaks = trans_breaks("log2", function(x) 2^x),
                                labels = trans_format("log2", math_format(2^.x))) +
                scale_y_continuous(trans = log2_trans(),
                                breaks = trans_breaks("log2", function(x) 2^x),
                                labels = trans_format("log2", math_format(2^.x))) +
                ggtitle(paste(nrow(dt.pt.specific), "Placenta Specific Protein-coding Genes")) +
                xlab("FPKM(Placenta)") + ylab("FPKM (Plasma - Illumina)") +
                theme_Publication()
        print(p2)
    }

    if(TRUE){
        library(scales) # for 'trans_breaks'
        p3<-ggplot(dt.pt.pl.specific, aes(Placenta,Plasma,label=hgnc_symbol)) +
                geom_point(size=5,alpha=.6) +
                geom_text(check_overlap=TRUE,vjust=-1,size=4) + 
                #coord_trans(y="log2", x="log2") +
                scale_x_continuous(trans = log2_trans(),
                                breaks = trans_breaks("log2", function(x) 2^x),
                                labels = trans_format("log2", math_format(2^.x))) +
                scale_y_continuous(trans = log2_trans(),
                                breaks = trans_breaks("log2", function(x) 2^x),
                                labels = trans_format("log2", math_format(2^.x))) +
                ggtitle(paste(nrow(dt.pt.pl.specific), "Placenta & Plasma Specific Protein-coding Genes")) +
                xlab("FPKM(Placenta)") + ylab("FPKM (Plasma - Illumina)") +
                theme_Publication()
        print(p3)
    }

    ###############################
    ## Gondon's rank-based method #
    ## PT or PL top 1 or top2     #
    ###############################
    dt.gtex.deseq[hgnc_symbol=="PSG7",.(Tissue,meanFpkm,rank=rank(meanFpkm)),"ensembl_gene_id,hgnc_symbol,gene_biotype"][order(hgnc_symbol,rank)]
    #
    my.ensg<-dt.gtex.deseq[,.(Tissue,meanFpkm,rank=rank(meanFpkm)),"ensembl_gene_id,hgnc_symbol,gene_biotype"][rank==21 | rank==22][Tissue %in% c("Placenta","Plasma"),.N,"ensembl_gene_id,hgnc_symbol"][N==2,as.character(ensembl_gene_id)]
    #dt.pt.pl.top.rank<-merge(dt.gtex.fpkm[ensembl_gene_id %in% my.ensg], dt.ensg[gene_biotype=="protein_coding",.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Placenta)]
    dt.pt.pl.top.rank<-merge(dt.gtex.fpkm[ensembl_gene_id %in% my.ensg], dt.ensg[,.(chromosome_name,ensembl_gene_id,hgnc_symbol,gene_biotype)])[order(-Placenta)]
    table(my.ensg %in% as.character(dt.pt.specific$ensembl_gene_id))
    table(my.ensg %in% as.character(dt.pt.specific.any$ensembl_gene_id))
    
    dt.pt.specific
}
