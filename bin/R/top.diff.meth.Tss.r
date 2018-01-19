#for each individual and for each Tss region (n=24) pull out the top 1%, 0.1%, 0.001% and 0.001% diff. methylated by % meth.diff
#################################################################################
#load R objects consiting of the methybase files from p2.regions2.r line 26
load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair2.diffMeth.Tss.RData")
rm(myDiff1,myDiff2,myDiff3)
quantile(abs(myDiffp1$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t5.a<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.b<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.c<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.d<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.a$region<-"Tss500"
p2t5.a$percent<-"1"
p2t5.b$region<-"Tss500"
p2t5.b$percent<-"0.1"
p2t5.c$region<-"Tss500"
p2t5.c$percent<-"0.001"
p2t5.d$region<-"Tss500"
p2t5.d$percent<-"0.0001"
Tss500<-rbind(p2t5.a,p2t5.b,p2t5.c,p2t5.d)
quantile(abs(myDiffp2$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t1.a<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.b<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.c<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.d<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.a$region<-"Tss1000"
p2t1.a$percent<-"1"
p2t1.b$region<-"Tss1000"
p2t1.b$percent<-"0.1"
p2t1.c$region<-"Tss1000"
p2t1.c$percent<-"0.001"
p2t1.d$region<-"Tss1000"
p2t1.d$percent<-"0.0001"
Tss1000<-rbind(p2t1.a,p2t1.b,p2t1.c,p2t1.d)
quantile(abs(myDiffp3$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t2.a<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.b<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.c<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.d<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.a$region<-"Tss2000"
p2t2.a$percent<-"1"
p2t2.b$region<-"Tss2000"
p2t2.b$percent<-"0.1"
p2t2.c$region<-"Tss2000"
p2t2.c$percent<-"0.001"
p2t2.d$region<-"Tss2000"
p2t2.d$percent<-"0.0001"
Tss2000<-rbind(p2t2.a,p2t2.b,p2t2.c,p2t2.d)
pair2<-rbind(Tss500,Tss1000,Tss2000)
save(pair2,file="top.diff.meth.Tss.RData")
rm(list=ls())
#######
#PAIR3#
#######
load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair3.diffMeth.Tss.RData")
rm(myDiff1,myDiff2,myDiff3)
quantile(abs(myDiffp1$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t5.a<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.b<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.c<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.d<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.a$region<-"Tss500"
p2t5.a$percent<-"1"
p2t5.b$region<-"Tss500"
p2t5.b$percent<-"0.1"
p2t5.c$region<-"Tss500"
p2t5.c$percent<-"0.001"
p2t5.d$region<-"Tss500"
p2t5.d$percent<-"0.0001"
Tss500<-rbind(p2t5.a,p2t5.b,p2t5.c,p2t5.d)
quantile(abs(myDiffp2$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t1.a<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.b<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.c<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.d<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.a$region<-"Tss1000"
p2t1.a$percent<-"1"
p2t1.b$region<-"Tss1000"
p2t1.b$percent<-"0.1"
p2t1.c$region<-"Tss1000"
p2t1.c$percent<-"0.001"
p2t1.d$region<-"Tss1000"
p2t1.d$percent<-"0.0001"
Tss1000<-rbind(p2t1.a,p2t1.b,p2t1.c,p2t1.d)
quantile(abs(myDiffp3$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t2.a<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.b<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.c<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.d<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.a$region<-"Tss2000"
p2t2.a$percent<-"1"
p2t2.b$region<-"Tss2000"
p2t2.b$percent<-"0.1"
p2t2.c$region<-"Tss2000"
p2t2.c$percent<-"0.001"
p2t2.d$region<-"Tss2000"
p2t2.d$percent<-"0.0001"
Tss2000<-rbind(p2t2.a,p2t2.b,p2t2.c,p2t2.d)
pair3<-rbind(Tss500,Tss1000,Tss2000)
load("top.diff.meth.Tss.RData")
save(pair2,pair3,file="top.diff.meth.Tss.RData")
rm(list=ls())
#######
#PAIR5#
######
load("top.diff.meth.Tss.RData")
load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair5.diffMeth.Tss.RData")
rm(myDiff1,myDiff2,myDiff3)
quantile(abs(myDiffp1$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t5.a<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.b<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.c<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.d<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.a$region<-"Tss500"
p2t5.a$percent<-"1"
p2t5.b$region<-"Tss500"
p2t5.b$percent<-"0.1"
p2t5.c$region<-"Tss500"
p2t5.c$percent<-"0.001"
p2t5.d$region<-"Tss500"
p2t5.d$percent<-"0.0001"
Tss500<-rbind(p2t5.a,p2t5.b,p2t5.c,p2t5.d)
quantile(abs(myDiffp2$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t1.a<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.b<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.c<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.d<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.a$region<-"Tss1000"
p2t1.a$percent<-"1"
p2t1.b$region<-"Tss1000"
p2t1.b$percent<-"0.1"
p2t1.c$region<-"Tss1000"
p2t1.c$percent<-"0.001"
p2t1.d$region<-"Tss1000"
p2t1.d$percent<-"0.0001"
Tss1000<-rbind(p2t1.a,p2t1.b,p2t1.c,p2t1.d)
quantile(abs(myDiffp3$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t2.a<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.b<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.c<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.d<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.a$region<-"Tss2000"
p2t2.a$percent<-"1"
p2t2.b$region<-"Tss2000"
p2t2.b$percent<-"0.1"
p2t2.c$region<-"Tss2000"
p2t2.c$percent<-"0.001"
p2t2.d$region<-"Tss2000"
p2t2.d$percent<-"0.0001"
Tss2000<-rbind(p2t2.a,p2t2.b,p2t2.c,p2t2.d)
pair5<-rbind(Tss500,Tss1000,Tss2000)
save(pair2,pair3,pair5,file="top.diff.meth.Tss.RData")
rm(list=ls())
#######
#PAIR8#
#######
load("top.diff.meth.Tss.RData")
load("/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/pair8.diffMeth.Tss.RData")
rm(myDiff1,myDiff2,myDiff3)
quantile(abs(myDiffp1$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t5.a<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.b<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.c<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.d<-subset(getData(myDiffp1),abs(myDiffp1$meth.diff) >= quantile(abs(myDiffp1$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t5.a$region<-"Tss500"
p2t5.a$percent<-"1"
p2t5.b$region<-"Tss500"
p2t5.b$percent<-"0.1"
p2t5.c$region<-"Tss500"
p2t5.c$percent<-"0.001"
p2t5.d$region<-"Tss500"
p2t5.d$percent<-"0.0001"
Tss500<-rbind(p2t5.a,p2t5.b,p2t5.c,p2t5.d)
quantile(abs(myDiffp2$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t1.a<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.b<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.c<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.d<-subset(getData(myDiffp2),abs(myDiffp2$meth.diff) >= quantile(abs(myDiffp2$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t1.a$region<-"Tss1000"
p2t1.a$percent<-"1"
p2t1.b$region<-"Tss1000"
p2t1.b$percent<-"0.1"
p2t1.c$region<-"Tss1000"
p2t1.c$percent<-"0.001"
p2t1.d$region<-"Tss1000"
p2t1.d$percent<-"0.0001"
Tss1000<-rbind(p2t1.a,p2t1.b,p2t1.c,p2t1.d)
quantile(abs(myDiffp3$meth.diff),c(0.99,0.999,0.9999,0.99999))
p2t2.a<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.b<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.c<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.9999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.d<-subset(getData(myDiffp3),abs(myDiffp3$meth.diff) >= quantile(abs(myDiffp3$meth.diff),0.99999),select=c("chr","start","end","strand","qvalue","meth.diff"))
p2t2.a$region<-"Tss2000"
p2t2.a$percent<-"1"
p2t2.b$region<-"Tss2000"
p2t2.b$percent<-"0.1"
p2t2.c$region<-"Tss2000"
p2t2.c$percent<-"0.001"
p2t2.d$region<-"Tss2000"
p2t2.d$percent<-"0.0001"
Tss2000<-rbind(p2t2.a,p2t2.b,p2t2.c,p2t2.d)
pair8<-rbind(Tss500,Tss1000,Tss2000)
save(pair2,pair3,pair5,pair8,file="top.diff.meth.Tss.RData")
rm(list=ls())
#############################################################
#Add Ensembl transcript ID and Gene ID to all the four files#
#############################################################
#Import dataframe R object from p2.regions2.r line 96 take Tss500 and name it uniquely
load("pair2.annotate.Tss.regions1.RData")
rm(f1,f2,f3)
p2.5<-subset(g1,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA2F,SGA2F))
rm(g1)
p2.5$region<-"Tss500"
p2.1<-subset(g2,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA2F,SGA2F))
rm(g2)
p2.1$region<-"Tss1000"
p2.2<-subset(g3,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA2F,SGA2F))
rm(g3)
p2.2$region<-"Tss2000"
a2<-rbind(p2.5,p2.1,p2.2)
p2.a<-merge(pair2,a2,by=c("chr","start","end","strand","qvalue","meth.diff","region"))
load("pair3.annotate.Tss.regions1.RData")
rm(f1,f2,f3)
p3.5<-subset(g1,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA3F,SGA3F))
rm(g1)
p3.5$region<-"Tss500"
p3.1<-subset(g2,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA3F,SGA3F))
rm(g2)
p3.1$region<-"Tss1000"
p3.2<-subset(g3,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA3F,SGA3F))
rm(g3)
p3.2$region<-"Tss2000"
a3<-rbind(p3.5,p3.1,p3.2)
p3.a<-merge(pair3,a3,by=c("chr","start","end","strand","qvalue","meth.diff","region"))
load("pair5.annotate.Tss.regions1.RData")
rm(f1,f2,f3)
p5.5<-subset(g1,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA5M,SGA5M))
rm(g1)
p5.5$region<-"Tss500"
p5.1<-subset(g2,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA5M,SGA5M))
rm(g2)
p5.1$region<-"Tss1000"
p5.2<-subset(g3,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA5M,SGA5M))
rm(g3)
p5.2$region<-"Tss2000"
a5<-rbind(p5.5,p5.1,p5.2)
p5.a<-merge(pair5,a5,by=c("chr","start","end","strand","qvalue","meth.diff","region"))
load("pair8.annotate.Tss.regions1.RData")
rm(f1,f2,f3)
p8.5<-subset(g1,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA8M,SGA8M))
rm(g1)
p8.5$region<-"Tss500"
p8.1<-subset(g2,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA8M,SGA8M))
rm(g2)
p8.1$region<-"Tss1000"
p8.2<-subset(g3,select=c(chr,start,end,strand,qvalue,meth.diff,feature.name,AGA8M,SGA8M))
rm(g3)
p8.2$region<-"Tss2000"
a8<-rbind(p8.5,p8.1,p8.2)
p8.a<-merge(pair8,a8,by=c("chr","start","end","strand","qvalue","meth.diff","region"))
save(pair2,pair3,pair5,pair8,p2.a,p3.a,p5.a,p8.a,file="top.diff.meth.Tss.RData")
