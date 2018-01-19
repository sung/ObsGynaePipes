#!/usr/bin/Rscript --vanilla
# http://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf

library("GenomicFeatures")

#############
## Ensembl ##
#############
#TxEnsembl<-makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") # BioMart dataset version: GRCh37.p13 (7/Oct/2014)
																								# defunct since Bioc 3.2
TxEnsembl<-makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                                  dataset="hsapiens_gene_ensembl",
                                  host="www.ensembl.org") # Ensembl ver: GRCh38.p3 (2/Dec/2015)
# Db type: TxDb
# Supporting package: GenomicFeatures
# Data source: BioMart
# Organism: Homo sapiens
# Resource URL: www.ensembl.org:80
# BioMart database: ENSEMBL_MART_ENSEMBL
# BioMart database version: Ensembl Genes 82
# BioMart dataset: hsapiens_gene_ensembl
# BioMart dataset description: Homo sapiens genes (GRCh38.p3)
# BioMart dataset version: GRCh38.p3
# Full dataset: yes
# miRBase build ID: NA
# transcript_nrow: 216133
# exon_nrow: 733504
# cds_nrow: 293948
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2015-12-03 15:18:38 +0000 (Thu, 03 Dec 2015)
# GenomicFeatures version at creation time: 1.20.6
# RSQLite version at creation time: 1.0.0
# DBSCHEMAVERSION: 1.1

enst<-transcripts(TxEnsembl)
ensx<-exons(TxEnsembl)
ensg<-genes(TxEnsembl)

# get GRanges of this enst
enst[mcols(enst)$tx_name=="ENST00000351205"]

# save to local file
my.local.db<-"~/data/Annotation/GRCh38.p3.TxEnsembl.sqlite"
saveDb(TxEnsembl,file=my.local.db)

# load from local file
txdb <- loadDb(my.local.db) # same as 'TxEnsembl'

# user defined annotation
makeTxDbPackageFromBiomart

##########
## UCSC ##
##########
supportedUCSCtables()
hg19KG <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
hg19ensGene <- makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")

# save to local file
my.local.db<-"~/data/Annotation/hg19ensGene.sqlite"
saveDb(hg19ensGene,file=my.local.db)

# for hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
myTxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Db type: TxDb
# Supporting package: GenomicFeatures
# Data source: UCSC
# Genome: hg38
# Organism: Homo sapiens
# UCSC Table: knownGene
# Resource URL: http://genome.ucsc.edu/
# Type of Gene ID: Entrez Gene ID
# Full dataset: yes
# miRBase build ID: NA
# transcript_nrow: 104178
# exon_nrow: 327656
# cds_nrow: 250176
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2015-03-19 13:57:20 -0700 (Thu, 19 Mar 2015)
# GenomicFeatures version at creation time: 1.19.32
# RSQLite version at creation time: 1.0.0
# DBSCHEMAVERSION: 1.1

# user defined annotation
makeTxDbPackageFromUCSC
