library(data.table)

dt.c_curve<-fread("~/results/SLX-9342.Homo_sapiens.PE75.v1/TopHat/D701_D501/accepted_hits.preseq.c_curve.txt")
dt.lc_extrap<-fread("~/results/SLX-9342.Homo_sapiens.PE75.v1/TopHat/D701_D501/accepted_hits.preseq.lc_extrap.txt")

dt.foo<-rbind(dt.c_curve[,.(`type`="observed",total_reads,distinct_reads)], dt.lc_extrap[,.(`type`="expected",`total_reads`=TOTAL_READS,`distinct_reads`=EXPECTED_DISTINCT)])
ggplot(dt.foo, aes(total_reads/10^6, distinct_reads/10^6)) + geom_point(aes(col=type),alpha=0.7,size=5) + xlim(c(0,100))
