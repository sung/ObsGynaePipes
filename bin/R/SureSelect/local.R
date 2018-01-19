library(data.table)
# ~/Pipelines/data/SureSelect/MJ.sureselect.meta.txt
# ~/Pipelines/bin/R/RnBeads/local.R (see process.PT.bed)
process.PT.bed<-function(my.cpg.type){
	cat(paste0("Processing PT:", my.cpg.type,"\n"))

	dt.meta<-fread("~/Pipelines/data/SureSelect/MJ.sureselect.meta.txt", stringsAsFactor=FALSE)

	##############################
	## Technical Replicate (TR) ##
	## 1 female (HAL03)
	## 1 male (HAL01)
	##############################
	print(system.time(dl.tech<-lapply(dt.meta[Condition=='AGA' & Context==my.cpg.type & Replicate=="TR",File], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)

	names(dl.tech)<-dt.meta[Condition=='AGA' & Context==my.cpg.type & Replicate=="TR",Barcode]

	#HAL03 (Female) HAL01 (Male)
	cat("\tmerging male and female...\n")
	PT.SS.Tech.CG.dt.merged<-merge(dl.tech[["HAL03"]][V1!='chrY' & V1!='chrM'], dl.tech[["HAL01"]][V1!='chrY' & V1!='chrM'], by=c("V1","V2","V3","V6"))
	# V1: Chr
	# V2: Start (0-based)
	# V3: End (1-based)
	# V4: %met
	# V5: Depth-of-coverage 
	# V6: Strand

	#> PT.SS.Tech.CG.dt.merged
    #       V1        V2        V3 V6      V4.x V5.x    V4.y V5.y
    #  1: chr1     13078     13079  +   0.00000    1   0.000    1
	#  2: chr1     15084     15085  +   0.00000    1 100.000    1

	PT.SS.Tech.CG.dt.merged[,c("V1","V2","V4.x","V4.y"):=list(substr(V1,4,5), NULL, round(V4.x*V5.x/100,2), round(V4.y*V5.y/100,2) )]
	setnames(PT.SS.Tech.CG.dt.merged, c("V1","V2","V3","V5.x","V6.x","V5.y","V6.y"))
	# V1:chr, V2:Pos, V3:Str, V5.x:cnt-met-female, V6.x:depth-cov-female, V5.y:cnt-met-male. V6.y:depth-cov-male
	#> PT.SS.Tech.CG.dt.merged
	#    V1        V2 V3 V5.x V6.x V5.y V6.y
	# 1:  8   3010282  +  0.0 12.0 15.1 22.5
	save(PT.SS.Tech.CG.dt.merged, file="~/data/RoadMap/BS-Seq/Schultz2015/PT.SS.Tech.CG.dt.merged.RData")

	###############################
	## Biological Replicate (BR) ##
	###############################
	print(system.time(dl.bio<-lapply(dt.meta[Condition=='AGA' & Context==my.cpg.type & Replicate=="BR",File], function(i){fread(i,sep="\t")[,.(V1,V2,V3,V4,V5,V6)]}))) # isa 'list' of data.table (by chr)
	names(dl.bio)<-dt.meta[Condition=='AGA' & Context==my.cpg.type & Replicate=="BR",Barcode]

	# aggregate 3 female replicates (HAL40, HAL71, HAL82)
	dt.bio.female<-rbindlist(dl.bio[dt.meta[Condition=='AGA' & Context==my.cpg.type & Replicate=="BR" & Sex=="F",Barcode]])[,list(`num.c`=round(sum(V4/100*V5),2),`depth`=sum(V5), .N),.(V1,V3,V6)]
	setnames(dt.bio.female, c("V1","V2","V3","V5","V6","N"))
	dt.bio.female[,V1:=substr(V1,4,5)] # chrX=>X chrMT=>MT

	# aggregate 3 male replicates
	dt.bio.male<-rbindlist(dl.bio[dt.meta[Condition=='AGA' & Context==my.cpg.type & Replicate=="BR" & Sex=="M",Barcode]])[,list(`num.c`=round(sum(V4/100*V5),2),`depth`=sum(V5), .N),.(V1,V3,V6)]
	setnames(dt.bio.male, c("V1","V2","V3","V5","V6","N"))
	dt.bio.male[,V1:=substr(V1,4,5)] # chrX=>X chrMT=>MT

	# merge females and males
	PT.SS.Bio.CG.dt.merged<-merge(dt.bio.female, dt.bio.male,by=c("V1","V2","V3"))
	if(FALSE){
		#dt.meta[Condition=='AGA' & Context==my.cpg.type & Replicate=="BR" & Sex=="F",Barcode]
		dt.foo<-merge(
					merge(dl.bio[["HAL40"]][V1!='chrY' & V1!='chrM'], dl.bio[["HAL71"]][V1!='chrY' & V1!='chrM'], by=c("V1","V2","V3","V6")), 
					dl.bio[["HAL82"]][V1!='chrY' & V1!='chrM'], 
					by=c("V1","V2","V3","V6")
					)
		#> dt.foo
		#           V1       V2       V3 V6     V4.x V5.x     V4.y V5.y       V4 V5
		#      1: chr1   167995   167996  +   0.0000    1 100.0000    1 100.0000  1
		#	   2: chr1   168016   168017  + 100.0000    1 100.0000    1 100.0000  1
		dt.bio.female<-dt.foo[,list(`num.c`=round(sum(V4.x/100*V5.x,V4.y/100*V5.y,V4/100*V5),2), `depth`=sum(V5.x,V5.y,V5)),"V1,V3,V6"]
		setnames(dt.bio.female, c("V1","V2","V3","V5","V6"))
		#> dt.bio.female
		#           V1        V2 V3 V5  V6
		#      1: chr1    167996  +  2   3

		# 3 Male technical replicate
		# HAL02, HAL26, HAL96
		dt.bar<-merge(
					merge(dl.bio[["HAL02"]][V1!='chrY' & V1!='chrM'], dl.bio[["HAL26"]][V1!='chrY' & V1!='chrM'], by=c("V1","V2","V3","V6")), 
					dl.bio[["HAL96"]][V1!='chrY' & V1!='chrM'], 
					by=c("V1","V2","V3","V6")
					)
		dt.bio.male<-dt.bar[,list(`num.c`=round(sum(V4.x/100*V5.x,V4.y/100*V5.y,V4/100*V5),2), `depth`=sum(V5.x,V5.y,V5)),"V1,V3,V6"]
		setnames(dt.bio.male, c("V1","V2","V3","V5","V6"))
		# merge 3 females and 3 males
		PT.SS.Bio.CG.dt.merged<-merge(dt.bio.female, dt.bio.male, by=c("V1","V2","V3"))
		PT.SS.Bio.CG.dt.merged[,c("V1"):=list(substr(V1,4,5))]
	}
	save(PT.SS.Bio.CG.dt.merged, file="~/data/RoadMap/BS-Seq/Schultz2015/PT.SS.Bio.CG.dt.merged.RData")
}
