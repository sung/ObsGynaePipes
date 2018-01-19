#!/USr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# First created 18/Mar/2016
# Last modified 18/Mar/2016

TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R") # load global config
source("~/Pipelines/config/graphic.R")
source("~/Pipelines/bin/R/RoadMap/local.R")

###################################################
# Overlap between promoter and cpgi by chromsome ##
###################################################
# http://stackoverflow.com/questions/27574775/is-it-possible-to-use-the-r-data-table-funcion-foverlaps-to-find-the-intersectio

# 0. genes
gr.genes<-get.region("Genes")
# 1. cpg islands
gr.cpgi<-get.region("CPGi")

# promoter 15.15
gr.pro.15.15<-get.region("Promo_15.15")
# promoter 15.10
gr.pro.15.10<-get.region("Promo_15.10")
# promoter 15.05
gr.pro.15.05<-get.region("Promo_15.05")
# promoter 15.00
gr.pro.15.00<-get.region("Promo_15.00")

# promoter 10.15
gr.pro.10.15<-get.region("Promo_10.15")
# promoter 10.10
gr.pro.10.10<-get.region("Promo_10.10")
# promoter 10.05
gr.pro.10.05<-get.region("Promo_10.05")
# promoter 10.00
gr.pro.10.00<-get.region("Promo_10.00")

# cpg islands shores
#gr.cgis<-get.region("CPGi-shores")
# enhancer (FANTOM5)
#gr.enhancer<-get.region("Enhancer")
#gr.enhancer.ts<-get.region("Enhancer-ts")

file.name<-file.path("~/results/RoadMap",paste("promoter.cpgi.overlap",time.stamp,"pdf",sep="."))
pdf(file=file.name, width=11.7, height=8.3, title="The size of overlap amongst gene, promoters, CPGi, CPGi-shroes")

df.overlap=plot.intersect(gr.pro.15.15, gr.cpgi, c("Promo_15.15","CPGi"), is.ignore.strand=TRUE)
Promo_15.15=data.frame(Promoter="Promo_15.15", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2]) 

df.overlap=plot.intersect(gr.pro.15.10, gr.cpgi, c("Promo_15.10","CPGi"), is.ignore.strand=TRUE)
Promo_15.10=data.frame(Promoter="Promo_15.10", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2])

df.overlap=plot.intersect(gr.pro.15.05, gr.cpgi, c("Promo_15.05","CPGi"), is.ignore.strand=TRUE)
Promo_15.05=data.frame(Promoter="Promo_15.05", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2])

df.overlap=plot.intersect(gr.pro.15.00, gr.cpgi, c("Promo_15.00","CPGi"), is.ignore.strand=TRUE)
Promo_15.00=data.frame(Promoter="Promo_15.00", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2])

df.overlap=plot.intersect(gr.pro.10.15, gr.cpgi, c("Promo_10.15","CPGi"), is.ignore.strand=TRUE)
Promo_10.15=data.frame(Promoter="Promo_10.15", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2])

df.overlap=plot.intersect(gr.pro.10.10, gr.cpgi, c("Promo_10.10","CPGi"), is.ignore.strand=TRUE)
Promo_10.10=data.frame(Promoter="Promo_10.10", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2])

df.overlap=plot.intersect(gr.pro.10.05, gr.cpgi, c("Promo_10.05","CPGi"), is.ignore.strand=TRUE)
Promo_10.05=data.frame(Promoter="Promo_10.05", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2])

df.overlap=plot.intersect(gr.pro.10.00, gr.cpgi, c("Promo_10.00","CPGi"), is.ignore.strand=TRUE)
Promo_10.00=data.frame(Promoter="Promo_10.00", Overlap.CPGi=colSums(df.overlap)[1] /  colSums(df.overlap)[3], Overlap.promoter=colSums(df.overlap)[1] /  colSums(df.overlap)[2])

foo<-rbind(Promo_15.15, Promo_15.10, Promo_15.05, Promo_15.00, Promo_10.15, Promo_10.10, Promo_10.05, Promo_10.00)

print(ggplot(foo, aes(Promoter,Overlap.CPGi*100)) + geom_bar(stat = "identity")  + ggtitle("% Overlap of CPGi with the promoters"))
print(ggplot(foo, aes(Promoter,Overlap.promoter*100)) + geom_bar(stat = "identity")  + ggtitle("% Overlap of Promoters with the CPGi"))

#plot.intersect(gr.pro.10.10, gr.enhancer, c("Promo_10.10","Enhancer"), is.ignore.strand=TRUE)
#plot.intersect(gr.genes, gr.enhancer, c("Genes","Enhancer"), is.ignore.strand=TRUE)
#plot.intersect(gr.cpgi, gr.enhancer, c("CPGi","Enhancer"), is.ignore.strand=TRUE)
#plot.intersect(gr.cgis, gr.enhancer, c("CPGi-shores","Enhancer"), is.ignore.strand=TRUE)

#plot.intersect(gr.enhancer, gr.enhancer.ts, c("Enhancer","Enhancer-ts"), is.ignore.strand=TRUE)

dev.off()
