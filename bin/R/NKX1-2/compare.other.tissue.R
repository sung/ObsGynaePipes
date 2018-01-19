#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

nkx.RData<-'/scratch/mdj30/methcall/methylkit/8075_8081/BisSNP/2/regions/nkx1.2.singleCG.locus.4pairs.5mC.RData'
load(nkx.RData) #the R object you should use is “all.3” 

tss.RData<-"~/scratch/results/Methyl-Seq/SGA.AGA/methylKit/dmr.tss.RData" # tss.5k, tss1k, tss2k, dmr.tss (common direction top tss)
load(tss.RData)
mj.sample<-c("AGA2F","SGA2F","AGA3F","SGA3F","AGA5M","SGA5M","AGA8M","SGA8M")

#NKX1-2
chr<-'chr10'
strand<-'-' # reverse strand
myENST<-"ENST00000440536"
if(myENST=="ENST00000451024"){
	start<-126135592;end<-126138753
}else{
	start<-126135998;end<-126138550
}
if(strand=='-'){tss.start=end}else{tss.start=start}
# http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000229544;r=1:10813-10814
# NKX1-2                     Chromosome 10: 126,135,592-126,138,753 reverse strand.
# NKX1-2-001 ENST00000451024 Chromosome 10: 126,135,592-126,138,753 reverse strand. 1580bp 310aa Protein coding CCDS59221 NM_001146340 NP_001139812 
# NKX1-2-201 ENST00000440536 Chromosome 10: 126,135,998-126,138,550 reverse strand. 1037bp 332aa Protein coding - -

public.data.list<-read.table(file='/home/ssg29/Pipelines/bin/R/NKX1-2/meta.txt', header=TRUE, stringsAsFactors=TRUE)

################
## Single CpG ##
################
nkx.cpg<-matrix(nrow=nrow(public.data.list),ncol=nrow(all.3)) # col: 14 CpG, rows: Tissue 
colnames(nkx.cpg)<-all.3$end # CpG locations
rownames(nkx.cpg)<-public.data.list$Tissue # Tissues

nkx.cpg.numC<-matrix(nrow=nrow(public.data.list),ncol=nrow(all.3)) # col: 14 CpG, rows: Tissue 
colnames(nkx.cpg.numC)<-all.3$end # CpG locations
rownames(nkx.cpg.numC)<-public.data.list$Tissue # Tissues

nkx.cpg.numT<-matrix(nrow=nrow(public.data.list),ncol=nrow(all.3)) # col: 14 CpG, rows: Tissue 
colnames(nkx.cpg.numT)<-all.3$end # CpG locations
rownames(nkx.cpg.numT)<-public.data.list$Tissue # Tissues

nkx.cpg.num<-matrix(nrow=nrow(public.data.list),ncol=nrow(all.3)) # col: 14 CpG, rows: Tissue 
colnames(nkx.cpg.num)<-all.3$end # CpG locations
rownames(nkx.cpg.num)<-public.data.list$Tissue # Tissues
# for each CpG
for(i in 1:nrow(all.3)){ 
	target=all.3[i,]$end
	print(target)
	#for(j in 1:nrow(public.data.list[public.data.list$Publication!='CamObsGynae',])){ # for each public tissue-specific CpG dataset
	for(j in 1:nrow(public.data.list)){ # for each public tissue-specific CpG dataset
		print(paste("processing",public.data.list[j,"Publication"], public.data.list[j,"Tissue"]))
		dummy<-paste0("awk '{if($1==\"chr10\"){if($6==\"-\" && $3==",target,")print $0}}' ", public.data.list[j,]$File)
		result<-system(dummy, intern=TRUE) # isa and should be either a character string (one row) or a vector (multiple rows)
		if(length(result)){
			numC<-NA; numT<-NA;
			print(paste("result=",result[1]))
			# PNAS
			# chr10   69163   69164   1.00    0       +       69163   69163   210,27,27
			if(public.data.list[j,"Publication"]=='PNAS-2013'){
				met=as.numeric(unlist(strsplit(result[1],"\t"))[4])*100
			# Ziller et al. 2013
			# chr1    10815   10817   '1/1'   1000    -
			}else if(public.data.list[j,"Publication"]=='Nature-2013'){
				met=unlist(strsplit(result[1],"\t"))[4] # "'1/1'"
				met=unlist(strsplit(met,"'"))[2] # " strip out "'"
				numC=as.numeric(unlist(strsplit(met,"/")))[1]
				numT=as.numeric(unlist(strsplit(met,"/")))[2] - numC
				met=round(numC/(numC+numT)*100,2) #do the calculation
				#met=round(as.numeric(unlist(strsplit(met,"/")))[1]/as.numeric(unlist(strsplit(met,"/")))[2]*100,2) #do the calculation
			# Encode
			#chr start end name score  strand tStart tEnd rgb count percentMeth
			# chr1    10003286        10003287        RRBS_CM01-201-001       20      -       10003286        10003287        0,255,0 20      0
			}else if(public.data.list[j,"Publication"]=='Encode'){
				met=as.numeric(unlist(strsplit(result[1],"\t"))[11])
				cov=as.numeric(unlist(strsplit(result[1],"\t"))[10]) #
				numC=round(met/100*cov); numT=cov-numC
			#MJ
			}else{
				met=as.numeric(unlist(strsplit(result[1],"\t"))[4]) # 20.6897 
				cov=as.numeric(unlist(strsplit(result[1],"\t"))[5]) # 29 
				numC=round(met/100*cov); numT=cov-numC
			}
			print(paste("%met=",met))
			nkx.cpg[j,i]=met
			nkx.cpg.num[j,i]=length(result)
			nkx.cpg.numC[j,i]=numC
			nkx.cpg.numT[j,i]=numT
		}#end of if length
	} # end of j
}# end of i 

dummy<-rbind(all.3$SGA2F,all.3$SGA3F,all.3$SGA5M,all.3$SGA8M)
rownames(dummy)<-c("SGA2F","SGA3F","SGA5M","SGA8M")
nkx.cpg<-rbind(nkx.cpg, dummy)

write.csv(nkx.cpg, "~/scratch/results/Methyl-Seq/SGA.AGA/NKX1-2.top.14.cpg.csv")

#######################################
### Avg % met within Tss & Gene body ##
#######################################
#awk 'BEGIN{cnt=0}{if($1=="chr10" && $6=="-" && $3<=126138550+2000 && $3>=126138550-100){pct+=$4; cnt++}}END{if(cnt!=0)print pct/cnt}' /home/ssg29/scratch/data/fastq/PNAS-2013/Placenta/GSM1083880_placenta2_MethylC-seq_chr10.hg19.BED
region.type<-c("tss2000","tss1000","tss500","body") 
nkx.region<-matrix(nrow=nrow(public.data.list),ncol=length(region.type)) # col: 14 CpG, rows: Tissue 
colnames(nkx.region)<-region.type
rownames(nkx.region)<-public.data.list$Tissue # Tissues

nkx.region.numC<-matrix(nrow=nrow(public.data.list),ncol=length(region.type)) # col: 14 CpG, rows: Tissue 
colnames(nkx.region.numC)<-region.type
rownames(nkx.region.numC)<-public.data.list$Tissue # Tissues

nkx.region.numT<-matrix(nrow=nrow(public.data.list),ncol=length(region.type)) # col: 14 CpG, rows: Tissue 
colnames(nkx.region.numT)<-region.type
rownames(nkx.region.numT)<-public.data.list$Tissue # Tissues

nkx.region.num<-matrix(nrow=nrow(public.data.list),ncol=length(region.type)) # col: 14 CpG, rows: Tissue 
colnames(nkx.region.num)<-region.type
rownames(nkx.region.num)<-public.data.list$Tissue # Tissues

for(i in 1:length(region.type)){ 
	print(paste0("processing ",region.type[i]))
	for(j in 1:nrow(public.data.list)){ # for each public tissue-specific CpG dataset
		print(paste("  processing",public.data.list[j,"Publication"], public.data.list[j,"Tissue"]))
		if(region.type[i]=="tss2000"){
			upper<-tss.start+2000
			lower<-tss.start-100
		}else if(region.type[i]=="tss1000"){
			upper<-tss.start+1000
			lower<-tss.start-100
		}else if(region.type[i]=="tss500"){
			upper<-tss.start+500
			lower<-tss.start-100
		}else{ # gene body
			upper<-end
			lower<-start
		}
		if(public.data.list[j,"Publication"]=='PNAS-2013'){
			dummy<-paste0("awk 'BEGIN{cnt=0}{if($1==\"chr10\" && $6==\"-\" && $3<=",upper," && $3>=",lower,"){pct+=$4*100; cnt++}}END{if(cnt!=0)print pct/cnt}' ", public.data.list[j,]$File)
		}else{
			dummy<-paste0("awk '{if($1==\"chr10\" && $6==\"-\" && $3<=",upper," && $3>=",lower,")print $0}' ", public.data.list[j,]$File)
		}
		result<-system(dummy, intern=TRUE) # isa and should be either a character string (one row) or a vector (multiple rows)
		if(length(result)){
			# chr10   69163   69164   1.00    0       +       69163   69163   210,27,27
			if(public.data.list[j,"Publication"]=='PNAS-2013'){
				avg.met=result[1]
			}else{
				sumC=0;sumT=0
				for(k in 1:length(result)){
					# Ziller et al Nature Letter 2013
					# chr1    10815   10817   '1/1'   1000    -
					if(public.data.list[j,"Publication"]=='Nature-2013'){
						met=unlist(strsplit(result[k],"\t"))[4] # "'1/1'"
						met=unlist(strsplit(met,"'"))[2] # " strip out "'"
						numC=as.numeric(unlist(strsplit(met,"/")))[1] # 
						numT=as.numeric(unlist(strsplit(met,"/")))[2]-numC # 
					# Encode
					#chr start end name score  strand tStart tEnd rgb count percentMeth
					# chr1    10003286        10003287        RRBS_CM01-201-001       20      -       10003286        10003287        0,255,0 20      0
					}else if(public.data.list[j,"Publication"]=='Encode'){
						met=as.numeric(unlist(strsplit(result[k],"\t"))[11])
						cov=as.numeric(unlist(strsplit(result[k],"\t"))[10]) #
						numC=round(met/100*cov); numT=cov-numC
					# MJ BSSNP bed format
					# chr1    10649   10650   20.6897 29      +       10649   10650   60,180,0        0       14
					}else{
						met=as.numeric(unlist(strsplit(result[k],"\t"))[4]) # 20.6897 
						cov=as.numeric(unlist(strsplit(result[k],"\t"))[5]) # 29 
						#if(cov>=10){ # only for depth>=10
						numC=round(met/100*cov); numT=cov-numC
						#}
					}
					sumC=sumC+numC; sumT=sumT+numT
				}
				print(paste0(" sumC=",sumC," sumT=",sumT))
				avg.met=round(sumC/(sumC+sumT)*100,2)
				nkx.region.numC[j,i]=sumC
				nkx.region.numT[j,i]=sumT
			} # end of publi.cdata.list 
			print(paste("%met=",avg.met))
			nkx.region[j,i]=avg.met
			nkx.region.num[j,i]=length(result)
		}#end of if length(result)
	} # end of j
}# end of length(region.type)

# NKX1-2-001 ENST00000451024 Chromosome 10: 126,135,592-126,138,753
# NKX1-2-201 ENST00000440536 Chromosome 10: 126,135,998-126,138,550
dummy<-rbind("tss2000"=dmr.tss$tss2000[dmr.tss$tss2000$feature.name==myENST,mj.sample],
	"tss1000"=dmr.tss$tss1000[dmr.tss$tss1000$feature.name==myENST,mj.sample],
	"tss500"=dmr.tss$tss500[dmr.tss$tss500$feature.name==myENST,mj.sample],
	"body"=NA)
dummy<-t(dummy)
nkx.region<-rbind(nkx.region, dummy)
write.csv(nkx.region, "~/scratch/results/Methyl-Seq/SGA.AGA/NKX1-2.region.comparision.csv")
