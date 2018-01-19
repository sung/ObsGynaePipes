#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>

###################################
## Collapse Technical Replicates ##
###################################
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = samples, design= deseq.design ) # isa 'DESeqDataSet'

cat("collapsing ddsHTSeq by SampleName...\n")
ddsCollapsed <- collapseReplicates(ddsHTSeq, groupby=ddsHTSeq$SampleName, run=ddsHTSeq$LibName) # isa 'DESeqDataSet'

cat("making dds...\n")
dds <- DESeq(ddsCollapsed) # isa 'DESeqDataSet'

cat("making rld...\n")
rld <- rlog(dds) # isa 'SummarizedExperiment'

####################################
# set exon size per gene for FPKM  #
####################################
hg19ensGene <- loadDb(my.local.db) # load annotation from local file
exons.list.per.gene <- exonsBy(hg19ensGene,"gene") # isa "GRangeList"
keep <- rownames(dds) %in% names(exons.list.per.gene) # remove ensg of no exons info
newCounts <- counts(dds)[keep,]
ddsMatrix <- DESeqDataSetFromMatrix(newCounts, colData(ddsCollapsed), deseq.design) # isa 'DESeqDataSet'
ddsExonicGene <- DESeq(ddsMatrix) # isa 'DESeqDataSet'

# 1. by rowData-way 
rowData(ddsExonicGene) <- exons.list.per.gene
# 2. by mcols[["basepairs"]] (method 1 is much faster) (https://www.biostars.org/p/83901/)
#exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))}) # isa 'list'
#mcols(ddsExonicGene)[["basepairs"]]<-unlist(exonic.gene.sizes[rownames(ddsExonicGene)]) # isa column 'vector'
#mcols(ddsExonicGene)$basepairs<-exonic.gene.sizes[rownames(dds),] # isa column 'vector'
ddsFpkm <-fpkm(ddsExonicGene) # isa 'matrix'

cat("saving DESeq2 RData\n")
save(ddsHTSeq,ddsCollapsed,dds,rld,ddsExonicGene,ddsFpkm, file=deseq.RData) #  save 'dds'
