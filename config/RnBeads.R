
time.stamp <- format(Sys.time(), '%Y-%m-%d_%I%p') #2015-01-07_01PM
doc <- 10 # depth of coverage
num.cores <-16  # Set Multi-cores
##########
# Config #
##########
my.project <- "SGA.AGA"
my.cpg.context <- "CpG"
#my.region <-'site'
my.region <- "enstTiling500"
#my.region <- "genomeTiling500"
my.title <- paste0(my.cpg.context,".",my.region) # CpG.enstTiling500

my.annotation.file <- "/new-data/ssg29/data/Annotation/hg19.enst.tiling500.csv.gz"
my.annotation.RData <- "/new-data/ssg29/data/Annotation/hg19.enst.tiling500.rnb.RData" # a path 
#my.annotation.file <- "/new-data/ssg29/data/Annotation/hg19.genome.tiling500.csv.gz"
#my.annotation.RData <- "/new-data/ssg29/data/Annotation/hg19.genome.tiling500.rnb.RData" # a path 

top.dir <- "/new-data/ssg29/results/RnBeads"
temp.dir <- "/new-data/ssg29/temp"
config.dir <- file.path(top.dir,"Meta") # RnBeads/Meta
sample.annotation <- file.path(config.dir, "cpg.sample.annotation.sbs.csv") # a table for the input bed files
data.source <- list(bed.dir=NULL, sample.sheet=sample.annotation,bed.column=1L)

project.dir <- file.path(top.dir,my.project); if(!file.exists(project.dir)){dir.create(project.dir)} # RnBeads/SGA.AGA
project.dir <- file.path(project.dir,my.cpg.context); if(!file.exists(project.dir)){dir.create(project.dir)} # RnBeads/SGA.AGA/CpG
run.dir <- file.path(project.dir,my.region); if(!file.exists(run.dir)){dir.create(run.dir)} # RnBeads/SGA.AGA/CpG/enstTiling500
report.dir <- file.path(run.dir, paste0("reports-",time.stamp)) # RnBeads/SGA.AGA/CpG/enstTiling500/reports-2015-01-07_01PM
save.dir <- file.path(project.dir, paste0(my.cpg.context,'.rnb.set')) # RnBeads/SGA.AGA/CpG/CpG.rnb.set.zip 

# Set log file
log.dir <- file.path(run.dir,"Log"); if(!file.exists(log.dir)){dir.create(log.dir)} # RnBeads/SGA.AGA/CpG/enstTiling500/Log
log.file <- file.path(log.dir,paste0(my.title,'.',time.stamp,".log")) # RnBeads/SGA.AGA/CpG/enstTiling500/Log/2015-01-07_01PM.log

options(fftempdir=temp.dir) # enable 'ff' package to write on disk
