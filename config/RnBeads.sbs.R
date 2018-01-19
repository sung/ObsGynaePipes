TR_PREFIX="GRCh37"
source("~/Pipelines/config/Annotation.R")
source("~/Pipelines/config/graphic.R")
options(scipen=999) # diable scientific notation
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I_%M_%p') #2015-01-07_01_30PM

##########
# Config #
##########
min.doc <- 10 # depth of coverage
num.cores <-20  # Set Multi-cores
my.cpg.context <- "CpG" # CpG context
#my.cpg.context <- "CHH" # CHH context
#my.cpg.context <- "CHG" # CHG context

##############################################
# WGBS 4 SGA vs. 5 AGA
#my.project <- "SGA.AGA"
#my.project <- "SGA.AGA.Boy"
#my.project <- "SGA.AGA.Girl"
##############################################

##############################################
# WGBS 2 Boys vs. 2 Girls (each SGA/AGA group)
my.project <- "AGA.Boy.Girl"
#my.project <- "SGA.Boy.Girl"
##############################################

##############################################
#my.project <- "BS.oxBS"
#my.project <- "5hmC.SGA.AGA"
##############################################

##############################################
# Fludigm v1 (SLX-8771)
#my.project <- "FLD.v1.8" # FLD v1: the original 8 samples
#my.project <- "FLD.v1.min.50" # FLD v1: low-pappa-sga 
#my.project <- "FLD.v1.SGA.PAPPA" # FLD v1: AGA Low PAPPA vs. AGA Normal PAPPA
#my.project <- "FLD.v1.PET" # FLD v1: PET vs. non-PET 
##############################################

##############################################
# Fludigm v2 (SLX-8772)
##############################################
#my.project <- "FLD.v2.the.eight" # Set1: the original 4 SGA vs. 4 AGA
#my.project <- "FLD.v2.sga.aga.validation" # Set2: additional 14 SGA vs. 14 AGA validation
#my.project <- "FLD.v2.aga.low.vs.normal.pappa" # Set3: 10 AGA low PAPP-A vs. 10 AGA normal PAPP-A
#my.project <- "FLD.v2.pet.control" # Set4: 8 PET vs. 8 control

##############################################
# Fludigm v3 (SLX-8773)
##############################################
#my.project <- "FLD.v3.the.eight" # Set1: the original 4 SGA vs. 4 AGA (pair8 NA for oxBS)
#my.project <- "FLD.v3.sga.aga.validation" # Set2: additional 15 SGA vs. 15 AGA validation
#my.project <- "FLD.v3.aga.low.vs.normal.pappa" # Set3: 10 AGA low PAPP-A vs. 10 AGA normal PAPP-A 
#my.project <- "FLD.v3.pet.control" # Set4: 8 PET vs. 8 control

##############################################
# Fludigm v4 (SLX-10295)
##############################################
#my.project <- "FLD.v4.the.eight" # Set1: the original 4 SGA vs. 4 AGA (pair8 NA for oxBS)
#my.project <- "FLD.v4.sga.aga.validation" # Set2: additional 15 SGA vs. 15 AGA validation
#my.project <- "FLD.v4.aga.low.vs.normal.pappa" # Set3: 10 AGA low PAPP-A vs. 10 AGA normal PAPP-A 
#my.project <- "FLD.v4.pet.control" # Set4: 8 PET vs. 8 control

##########################
# SureSelect (SLX-10409) #
##########################
#my.project <- "SureSelect.Boy.Girl" # either Set1 or Set2


## Common Varaibles
top.dir <- "~/results/RnBeads"
meta.dir <- file.path(top.dir,"Meta") # RnBeads/Meta
temp.dir <- file.path("~/temp",time.stamp); if(!file.exists(temp.dir)){dir.create(temp.dir)} 

# Define this if you have already preprocessed methylation rnb data set 
if(my.project=="SGA.AGA"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-01-21_04_30_PM' # for SGA.AGA.all or Boys.Girls (for the inital run)
	my.meta.file <- file.path(meta.dir,"cpg.sample.annotation.sbs.csv")
}else if(my.project=="SGA.AGA.Boy"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-05-13_04_24_PM'
	my.meta.file <- file.path(meta.dir,"cpg.sample.annotation.boy.csv")
}else if(my.project=="SGA.AGA.Girl"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-05-15_05_09_PM' 
	my.meta.file <- file.path(meta.dir,"cpg.sample.annotation.girl.csv")
}else if(my.project=="5hmC.SGA.AGA"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-03-27_05_04_PM'
	my.meta.file <- file.path(meta.dir,"cpg.sample.annotation.5hmC.csv") # 5hmC.SGA.AGA
}else if(my.project=="AGA.Boy.Girl"){
	my.group.list=c("FetalSex")
	my.run.date<-'reports-2015-03-30_04_47_PM'
	my.meta.file <- file.path(meta.dir,"Boy.Girl",paste0(my.cpg.context,".sample.annotation.sbs.AGA.csv")) # Boys.Girls (AGA CHH)
}else if(my.project=="SGA.Boy.Girl"){
	my.group.list=c("FetalSex")
	my.run.date<-'reports-2015-04-09_06_08_PM' # Boys.Girls (SGA)
	my.meta.file <- file.path(meta.dir,"Boy.Girl","cpg.sample.annotation.sbs.SGA.csv") # Boys.Girls (SGA)
}else if(my.project=="BS.oxBS"){
	my.group.list=c("SeqType")
	my.meta.file <- file.path(meta.dir,"cpg.sample.annotation.sbs.AGA.5hmC.csv") # BS.oxBS AGA2F.BS - AGA2F.oxBS
}else if(my.project=="FLD.v1.8"){
	my.group.list=c("Condition")
	#my.run.date<-'reports-2015-04-21_02_36_PM' # for CpG coverage (QC)
	my.run.date<-'reports-2015-04-22_04_20_PM' # Set3: for diff meth (n=8)
	my.meta.file <- file.path(meta.dir,'FLD.v1',paste0("cpg.sample.annotation.",my.project,".csv")) #cpg.sample.annotation.SGA.AGA.FLD.v1.8.csv
}else if(my.project=="FLD.v1.min.50"){
	my.group.list=c("Condition")
	#my.run.date<-'reports-2015-04-23_05_08_PM' # Set2: low-pappa-sga (n=20) 
	my.run.date<-'reports-2015-04-23_06_42_PM' # Set3: low-pappa-sga 7 normal-pappa-aga (n=14)
	my.meta.file <- file.path(meta.dir,'FLD.v1',paste0(my.project,".csv"))
}else if(my.project=="FLD.v1.SGA.PAPPA"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-04-28_03_03_PM' # Set4: AGA Low PAPPA vs. AGA Normal PAPPA
	my.meta.file <- file.path(meta.dir,'FLD.v1',"FLD.aga.low.normal.pappa.min.45.cpg.csv")
}else if(my.project=="FLD.v1.PET"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-04-28_03_06_PM' # Set5: AGA Low PAPPA vs. AGA Normal PAPPA
	my.meta.file <- file.path(meta.dir,'FLD.v1',"FLD.pet.control.min.50.cpg.csv" )
}else if(my.project=="FLD.v2.the.eight"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-06-12_06_47_PM' 
	my.meta.file <- file.path(meta.dir,'FLD.v2',paste0(my.project,".min.50.cpg.csv"))
}else if(my.project=="FLD.v2.sga.aga.validation"){ 
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-06-15_12_50_PM'
	my.meta.file <- file.path(meta.dir,'FLD.v2',paste0(my.project,".min.50.cpg.csv"))
}else if(my.project=="FLD.v2.aga.low.vs.normal.pappa"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-06-15_12_51_PM'
	my.meta.file <- file.path(meta.dir,'FLD.v2',paste0(my.project,".min.50.cpg.csv"))
}else if(my.project=="FLD.v2.pet.control"){
	my.group.list=c("Condition")
	my.run.date<-'reports-2015-06-15_02_02_PM'
	my.meta.file <- file.path(meta.dir,'FLD.v2',paste0(my.project,".min.50.cpg.csv"))
}else if(my.project=="FLD.v3.the.eight"){
	my.group.list=c("Condition")
	#my.run.date<-'reports-2015-08-10_04_02_PM' # oxBS, SLX-8773.Homo_sapiens.v1 
	#my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.oxBS.csv")) # reports-2015-08-10_04_02_PM
	#my.run.date<-'reports-2015-08-10_02_15_PM' # BS, SLX-8773.Homo_sapiens.v1 
	#my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.BS.csv")) # reports-2015-08-10_02_15_PM
	#my.run.date <- 'reports-2015-08-14_09_41_AM' # oxBS, SLX-8773.Homo_sapiens.v2	
	#my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v2',paste0(my.project,".cpg.oxBS.csv")) # reports-2015-08-14_09_41_AM
	my.run.date <- 'reports-2015-08-14_09_44_AM' # BS, SLX-8773.Homo_sapiens.v1
	my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v2',paste0(my.project,".cpg.BS.csv")) # reports-2015-08-14_09_44_AM
}else if(my.project=="FLD.v3.sga.aga.validation"){
	my.group.list=c("Condition")
	#my.run.date <- 'reports-2015-08-11_12_15_PM' # oxBS, SLX-8773.Homo_sapiens.v1
	#my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.oxBS.csv")) # reports-2015-08-11_12_15_PM 
	my.run.date <- 'reports-2015-08-11_12_20_PM' # BS, SLX-8773.Homo_sapiens.v1
	my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.BS.csv")) # reports-2015-08-11_12_20_PM
}else if(my.project=="FLD.v3.aga.low.vs.normal.pappa"){
	my.group.list=c("Condition")
	my.run.date <- 'reports-2015-08-13_01_36_PM' #oxBS, SLX-8773.Homo_sapiens.v1
	my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.oxBS.csv")) # reports-2015-08-13_01_36_PM
	#my.run.date <- 'reports-2015-08-12_11_01_AM' #BS, SLX-8773.Homo_sapiens.v1
	#my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.BS.csv")) # reports-2015-08-12_11_01_AM 
}else if(my.project=="FLD.v3.pet.control"){
	my.group.list=c("Condition")
	#my.run.date <- 'reports-2015-08-13_01_30_PM' # oxBS, SLX-8773.Homo_sapiens.v1
	#my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.oxBS.csv")) # reports-2015-08-13_01_30_PM
	my.run.date <- 'reports-2015-08-13_01_33_PM' # BS, SLX-8773.Homo_sapiens.v1
	my.meta.file <- file.path(meta.dir,'FLD.v3','SLX-8773.Homo_sapiens.v1',paste0(my.project,".cpg.BS.csv")) # reports-2015-08-13_01_33_PM
}else if(my.project=="FLD.v4.the.eight"){
	my.group.list=c("Condition")
	my.run.date <- 'reports-2015-09-11_03_34_PM' # oxBS, SLX-10295.v2 
	my.meta.file <- file.path(meta.dir,'FLD.v4',paste0(my.project,".cpg.csv"))
}else if(my.project=="FLD.v4.sga.aga.validation"){
	my.group.list=c("Condition")
	my.run.date <- 'reports-2015-09-18_01_38_PM'
	my.meta.file <- file.path(meta.dir,'FLD.v4',paste0(my.project,".cpg.csv"))
}else if(my.project=="FLD.v4.aga.low.vs.normal.pappa"){
	my.group.list=c("Condition")
	my.run.date <- 'reports-2015-09-28_05_22_PM'
	my.meta.file <- file.path(meta.dir,'FLD.v4',paste0(my.project,".cpg.csv"))
}else if(my.project=="SureSelect.Boy.Girl"){
	my.group.list=c("FetalSex")
	#my.run.date<-'reports-2016-02-05_01_44_PM'
	#my.meta.file <- file.path(meta.dir,my.project,"Set1.cpg.sample.annotation.csv")  # reports-2016-02-05_01_44_PM 
	my.run.date<-'reports-2016-02-04_04_21_PM'
	my.meta.file <- file.path(meta.dir,my.project,"Set2.cpg.sample.annotation.csv") # reports-2016-02-04_04_21_PM 
}else{
	stop(paste0(my.project, " not supported. Stopped!"))
}

# Set annotation file
if(!file.exists(my.meta.file)){ # RnBeads/SGA.AGA/Meta/cpg.sample.annotation.sbs.csv
	stop(paste0(my.meta.file, " not found. Stoppped"))
}

my.defined.regions <- c("genomeTiling500")
my.annotation.csv<- c('~/data/Annotation/RnBeads/hg19.genome.tiling500.csv.gz')
my.annotation.rdata<- c('~/data/Annotation/RnBeads/hg19.genome.tiling500.rnb.RData')

my.annotations <- data.frame(
			csv.file=my.annotation.csv,
			rdata.file=my.annotation.rdata,
			row.names=my.defined.regions,
			stringsAsFactors=FALSE
)

my.region <- "all" # enstTiling500 genomeTiling500
#my.region <- "enstTiling500" # genomeTiling500

my.title <- paste(my.project,my.cpg.context,my.region,sep='.') # CpG.enstTiling500

data.source <- list(bed.dir=NULL, sample.sheet=my.meta.file, bed.column=1L)

project.dir <- file.path(top.dir,my.project); if(!file.exists(project.dir)){dir.create(project.dir)} # RnBeads/SGA.AGA
project.dir <- file.path(project.dir,my.cpg.context); if(!file.exists(project.dir)){dir.create(project.dir)} # RnBeads/SGA.AGA/CpG
rdata.dir <- file.path(project.dir,'RData'); if(!file.exists(rdata.dir)){dir.create(rdata.dir)} # RnBeads/SGA.AGA/CpG/RData
run.dir <- file.path(project.dir,my.region); if(!file.exists(run.dir)){dir.create(run.dir)} # RnBeads/5hmC.SGA.AGA/CpG/all
report.dir <- file.path(run.dir, paste0("reports-",time.stamp)) # RnBeads/SGA.AGA/CpG/enstTiling500/reports-2015-01-07_01PM
rnb.save.dir <- file.path(run.dir, paste0(my.cpg.context,'.rnb.set')) # RnBeads/SGA.AGA/CpG/enstTiling500/CpG.rnb.set(.zip)

if(exists("my.run.date")){
	cat(paste0("Run Date:",my.run.date,"\n"))
	if(file.exists(file.path(run.dir, paste0(my.run.date,'/rnbSet_preprocessed.zip')))){
		preprocessed.rnb.set.dir<-file.path(run.dir, paste0(my.run.date,'/rnbSet_preprocessed.zip')) 
		cat("preprocessed.rnb.set.dir loaded\n")
	}
}

# Set log file
log.dir <- file.path(run.dir,"Log"); if(!file.exists(log.dir)){dir.create(log.dir)} # RnBeads/SGA.AGA/CpG/enstTiling500/Log
log.file <- file.path(log.dir,paste0(my.title,'.',time.stamp,".log")) # RnBeads/SGA.AGA/CpG/enstTiling500/Log/2015-01-07_01PM.log

options(fftempdir=temp.dir) # enable 'ff' package to write on disk
unixtools::set.tempdir(temp.dir) # unixtools downloaded from https://rforge.net/doc/packages/unixtools
cat("Temp dir:\n")
print(tempdir())

cat(paste0("Running RnBeads for ",my.title,"\n"))

