#!/usr/bin/Rscript --vanilla
#http://code.google.com/p/methylkit

source ("~/Pipelines/config/methylkit.R") # load config
#############################
## 1. Read methylation files
#############################
#Set2, 68(SGA) ~/scratch/results/SLX-8074.SLX-8080.v1/BisSNP/A003/SLX-8074.SLX-8080.A003.cpg.filtered.CG.methylkit.txt
#Set2, 66(AGA) ~/scratch/results/SLX-8074.SLX-8080.v1/BisSNP/A001/SLX-8074.SLX-8080.A001.cpg.filtered.CG.methylkit.txt
#Set3, 74(SGA) ~/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A007/SLX-8075.SLX-8077.SLX-8081.A007.cpg.filtered.CG.methylkit.txt
#Set3, 72(AGA) ~/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A005/SLX-8075.SLX-8077.SLX-8081.A005.cpg.filtered.CG.methylkit.txt
#Set5, 77(SGA) ~/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A004/SLX-8075.SLX-8077.SLX-8081.A004.cpg.filtered.CG.methylkit.txt
#Set5, 75(AGA) ~/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A002/SLX-8075.SLX-8077.SLX-8081.A002.cpg.filtered.CG.methylkit.txt
#Set8, 65(SGA) ~/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A006/SLX-8075.SLX-8077.SLX-8081.A006.cpg.filtered.CG.methylkit.txt
#Set8, 73(AGA) ~/scratch/results/SLX-8075.SLX-8077.SLX-8081.v1/BisSNP/A008/SLX-8075.SLX-8077.SLX-8081.A008.cpg.filtered.CG.methylkit.txt

# Below for Set2 only
if(FALSE){
	file.list2=list( paste0(bissnp_call_dir1,"/A003/",project1,".A003.cpg.filtered.CG.methylkit.txt"),
					paste0(bissnp_call_dir1,"/A001/",project1,".A001.cpg.filtered.CG.methylkit.txt"))
	myobj2=read( file.list2, sample.id=list("SGA2F","AGA2F"),
				assembly="hg19",treatment=c(1,0),context="CpG") # is a 'methylRawList'

	file.list3=list(paste0(bissnp_call_dir2,"/A007/",project2,".A007.cpg.filtered.CG.methylkit.txt"),
					paste0(bissnp_call_dir2,"/A005/",project2,".A005.cpg.filtered.CG.methylkit.txt"))
	myobj3=read( file.list3, sample.id=list("SGA3F","AGA3F"),
				assembly="hg19",treatment=c(1,0),context="CpG") # is a 'methylRawList'

	file.list5=list(paste0(bissnp_call_dir2,"/A004/",project2,".A004.cpg.filtered.CG.methylkit.txt"),
					paste0(bissnp_call_dir2,"/A002/",project2,".A002.cpg.filtered.CG.methylkit.txt"))
	myobj5=read( file.list5, sample.id=list("SGA5M","AGA5M"),
				assembly="hg19",treatment=c(1,0),context="CpG") # is a 'methylRawList'

	file.list8=list(paste0(bissnp_call_dir2,"/A006/",project2,".A006.cpg.filtered.CG.methylkit.txt"),
					paste0(bissnp_call_dir2,"/A008/",project2,".A008.cpg.filtered.CG.methylkit.txt"))
	myobj8=read( file.list8, sample.id=list("SGA8M","AGA8M"),
				assembly="hg19",treatment=c(1,0),context="CpG") # is a 'methylRawList'
}

# Below for Set2, Set3, Set5 and Set8
file.list=list( paste0(bissnp_call_dir1,"/A003/",project1,".A003.cpg.filtered.CG.methylkit.txt"),
				paste0(bissnp_call_dir1,"/A001/",project1,".A001.cpg.filtered.CG.methylkit.txt"),
				paste0(bissnp_call_dir2,"/A007/",project2,".A007.cpg.filtered.CG.methylkit.txt"),
				paste0(bissnp_call_dir2,"/A005/",project2,".A005.cpg.filtered.CG.methylkit.txt"),
				paste0(bissnp_call_dir2,"/A004/",project2,".A004.cpg.filtered.CG.methylkit.txt"),
				paste0(bissnp_call_dir2,"/A002/",project2,".A002.cpg.filtered.CG.methylkit.txt"),
				paste0(bissnp_call_dir2,"/A006/",project2,".A006.cpg.filtered.CG.methylkit.txt"),
				paste0(bissnp_call_dir2,"/A008/",project2,".A008.cpg.filtered.CG.methylkit.txt"))

#myobj is a class of 'methylRaw' ('methylRawList')
myobj=read( file.list, sample.id=list("SGA2F","AGA2F","SGA3F","AGA3F","SGA5M","AGA5M","SGA8M","AGA8M"),
			assembly="hg19",treatment=c(1,0,1,0,1,0,1,0),context="CpG") # is a 'methylRawList'

myobj2=reorganize(myobj, sample.id=list("SGA2F","AGA2F"),treatment=c(1,0)) # is a 'methylRawList'
myobj3=reorganize(myobj, sample.id=list("SGA3F","AGA3F"),treatment=c(1,0)) # is a 'methylRawList'
myobj5=reorganize(myobj, sample.id=list("SGA5M","AGA5M"),treatment=c(1,0)) # is a 'methylRawList'
myobj8=reorganize(myobj, sample.id=list("SGA8M","AGA8M"),treatment=c(1,0)) # is a 'methylRawList'

save(myobj2,myobj3,myobj5,myobj8,myobj,file=paste0(methylkit_dir,"/",my_prefix,".myobj.RData")) # 8oxBS.myobj.RData
