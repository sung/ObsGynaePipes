
	## Load RoadMap Tissue
	dt.tissue.locus<-list()
	for(my.target in c("CSMD1.DMR2","CSMD1.DMR","CSMD1.legacy.promo","CSMD1.putative.promo")){
		my.gr.ensg <- reduce(rtracklayer::import.bed(target.region[[my.target]]))
		names(my.gr.ensg)<-my.target

		# RoadMap
		subject.keys=c("seqnames","start","end")
		query.key=c("V1","V2","End")
		subject<-as.data.frame(my.gr.ensg)
		cat("\tRemoving leading 'chr'...\n")
		#subject$seqnames<-simplify2array(strsplit(as.character(subject$seqnames), "chr"))[2,] # chrX = > X
		subject$seqnames<-substr(subject$seqnames, 4, 5) #chrX=>X
		cat("\tconvert to data.table...\n")
		dt.subject<-as.data.table(subject) # isa data.table
		setkeyv(dt.subject,subject.keys) # regions of interests 
		rm(subject)
		cat("Subject loaded\n")

		for(my.tissue in c(avail.tissues[!avail.tissues %in% "PA"])){ # except 'PT': placenta
			## Load this RoadMap tissue
			dt.roadmap=load.my.tissue.dt.merged(my.tissue,my.cpg.type)

			cat("\tFinding CpGs overlap...\n")
			system.time(dt.overlap<-foverlaps(dt.roadmap, dt.subject, by.x=query.key, type="any", nomatch=0L))
			#> dt.overlap
			#   V1    Start      End Strand       V2 V3  V4 V7 V5.x V6.x V5.y V6.y    i.End
			#1:  X  2746554  2748235      *  2746590  - CGA  1    6   17    2   10  2746590
			#2:  X  2835950  2836236      *  2836235  + CGG  1   14   14    9   10  2836235
			if(nrow(dt.overlap)){
				dt.meth<-rbind(dt.overlap[,.(V1,V2,V5.x/V6.x*100,"Female",my.tissue)], dt.overlap[,.(V1,V2,V5.y/V6.y*100,"Male",my.tissue)]) # isa 'data.table'
				setnames(dt.meth,c("chr","start","methylation_level","Gender","Tissue"))
				dt.tissue.locus[[my.target]][[my.tissue]]<-dt.meth
			}
		}# end of my.tissue 
	}
