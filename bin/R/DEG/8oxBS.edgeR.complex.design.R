#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# For detial, see http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# complex design 14.B based on the publication above

myCaller="8oxBS"
source ("~/Pipelines/config/DEG.R") # load config

# pdf output filename 
pdf_file_name<- paste0(edgeRHome,'/',sampleType,'/.plot')
pdf (file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

# 1. Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count
rawCounts = readDGE(samples$HtseqFile, header=FALSE, sep="")$counts
colnames(rawCounts) = samples$SampleName
#> dim(rawCounts)
#[1] 63681   112
#> length(grep("^ENSG", rownames(rawCounts)))
#[1] 63676

# 2. Filter weakly expressed and noninformative (e.g., non-aligned) features using a command like:
noint = rownames(rawCounts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
cpms = cpm(rawCounts)
keep = rowSums(cpms >1) >= nrow(samples)/2 & !noint
#> table(keep)
#keep
#FALSE  TRUE 
#50659 13022
counts = rawCounts[keep,]
nrow(counts)

# 3. Visualize and inspect the count table
head( counts[,order(samples$Condition)], 5 ) # print first 5 case (SGA) records
summary(counts)

# 4. Create a DGEList object (edgeR's container for RNA-seq count data)
d = DGEList(counts=counts, group=samples$Condition) # isa 'DGEList'

# 5. Estimate normalization factors
# The calcNormFactors function normalizes for RNA composition by finding a set of scaling
# factors for the library sizes that minimize the log-fold changes between the samples for most genes
# The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples 
d = calcNormFactors(d)
dim(d)

# 6. Inspect the relationships between samples using a multidimensional scaling (MDS) plot
plotMDS(d, labels=samples$SampleName, col=c("darkgreen","blue")[factor(samples$Condition)])

# 7. Create a design matrix to specify the factors that are expected to affect expression levels:
design <- model.matrix(~ Pair + Condition, samples) # isa 'matrix'

# 8. Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood
d2 = estimateGLMCommonDisp(d, design, verbose=TRUE) # is a DGEList
d2 = estimateGLMTrendedDisp(d2, design, verbose=TRUE)
d2 = estimateGLMTagwiseDisp(d2, design)

# d2$counts is same with counts
# cpm(d2) is different from cpm(counts)
# cpm(d2)["ENSG00000174697",]
# cpm(counts)["ENSG00000174697",]

plotMeanVar(d2, show.tagwise.vars=TRUE, NBline=TRUE)
plotBCV(d2)

# 9. Given the design matrix and dispersion estimates, fit a GLM to each feature:
f = glmFit(d2, design) # is a DGEGLM

# 10. Perform a likelihood ratio test, specifying the difference of interest 
# here, case versus control
# glmLRT(glmfit, coef=ncol(glmfit$design), contrast=NULL, test="chisq")
# test: which test (distribution) to use in calculating the p-values. Possible values are ‘"F"’ or ‘"chisq"’.
# If ‘coef’ is used, the null hypothesis is that all the coefficients indicated by ‘coef’ are equal to zero.
# The data frame ‘table’ contains the following columns:
#   logFC: log2-fold change of expression between conditions being tested.
#  logCPM: average log2-counts per million, the average taken over all libraries in ‘y’.
#      LR: likelihood ratio statistics (only for ‘glmLRT’).
#       F: F-statistics (only for ‘glmQFTest’).
#  PValue: p-values.
de = glmLRT(f)  # the column number for Condition, which is ConditionCTR in our case, in design data.frame
				# isa "DGELRT"
				# check out de$table (logFC, logCPM, LR, PValue)
				# check out de$samples
				# check out de$comparison

o = order(de$table$PValue)
cpm(d2)[o[1:10],]

# 11. Use the topTags function to present a tabular summary of the differential expression statistics 
# topTags(object, n=10, adjust.method="BH", sort.by="PValue")
# topTags extracts the top DE tags in a data frame for a given pair of
# groups, ranked by p-value or absolute log-fold change.
# (note that topTags operates on the output of exactTest or glmLRT)
tt = topTags(de, n=nrow(d2))    # default sort by pvalue
								# is a TopTags
								# BH is a default adjust.method
								# check out tt$comparison (should be the same with de$comparison)
#tt2 = topTags(de, n=nrow(d2), adjust.method="bonferroni")
#tt2 = topTags(de, n=nrow(d2), sort.by="logFC", adjust.method="bonferroni")
edgeR.deg=tt$table
head(edgeR.deg) # isa data.frame

# 12. Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:
nc = cpm(d2, normalized.lib.sizes=TRUE) # isa 'matrix'
rn = rownames(edgeR.deg)
head(nc[rn,order(samples$Condition)],5) # print first 5 entries

# 13. Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot, 
# here showing the genes selected as differentially expressed (with a 5% false discovery rate
table(edgeR.deg$logFC<0)
plotSmear(de, col="#00000020", cex=0.5, main=paste0(nrow(edgeR.deg), " genes (", table(edgeR.deg$logFC<0)["TRUE"], "<0, ", table(edgeR.deg$logFC<0)["FALSE"], ">=0)" ))
abline(h=c(0), col="blue")

################################################
# conventional filtering |logFC|>1 & FDR < 0.1 #
################################################
#deg = rn[(edgeR.deg$logFC <= -1 | edgeR.deg$logFC >= 1) & edgeR.deg$FDR<=0.1]
deg = rn[edgeR.deg$FDR <= 0.1]
#edgeR.deg[deg,]
plotSmear(de, de.tags=deg, col="#00000020", cex=0.5, main="conventional (FC>=1 and FDR<=0.1)")
abline(h=c(-1, 1), col="blue")

# P value distribution
hist(edgeR.deg$PValue, breaks=100, xlab='P-value', main=paste0("P-value distribution of ", nrow(edgeR.deg), " genes"))

#############
# by Pvalue #
#############
table(edgeR.deg$PValue <= myPValue)
deg = rn[edgeR.deg$PValue <= myPValue ]
plotSmear(de, de.tags=deg, col="#00000020", cex=0.5, main=paste0(length(deg)," genes (", table(edgeR.deg[deg,]$logFC<0)["TRUE"], "<0, ", table(edgeR.deg[deg,]$logFC<0)["FALSE"], ">=0) having P-value <=", myPValue))
abline(h=c(0), col="blue")

# FDR distribution of genes having P-value <= myPValue
hist(edgeR.deg[deg,]$FDR, breaks=20, xlab='FDR', main=paste0("FDR distribution of ", length(deg), "genes having P-value <=", myPValue))

####################
# by Pvalue && FDR #
####################
table(Pvalue=(edgeR.deg$PValue <= myPValue), FDR=(edgeR.deg$FDR <= myFDR))
deg = rn[edgeR.deg$PValue <= myPValue & edgeR.deg$FDR <= myFDR ]
# get the list of genes having P-value <= myPvalue and FDR <= myFDR
edgeR.deg[order(edgeR.deg[deg,]$logFC),]
# get the min FC from above
minlogFC=min(abs(edgeR.deg[deg,]$logFC))

plotSmear(de, de.tags=deg, col="#00000020", cex=0.5, main=paste0(length(deg), " genes (", table(edgeR.deg[deg,]$logFC<0)["TRUE"], "<0, ", table(edgeR.deg[deg,]$logFC<0)["FALSE"], ">=0) of FDR <= ", myFDR, " (min logFC=", round(minlogFC,3), ")"))
abline(h=c(-minlogFC, minlogFC), col="blue")

# 14. Save the result table as a CSV file (alternative formats are possible) as follows:
#write.csv(edgeR.deg, file=paste0(edgeRHome,"/toptags_all_pair_edgeR.",sampleType,".csv")

################
## Annotation ##
################
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
top.deg=edgeR.deg[deg,]
top.deg$ensg<-rownames(top.deg)

#fields <- c("ensembl_gene_id", "hgnc_symbol","description", "uniprot_swissprot", "gene_biotype", "phenotype_description", "goslim_goa_accession", "goslim_goa_description", "go_id", "name_1006", "definition_1006", "namespace_1003")
fields <- c("ensembl_gene_id", "hgnc_symbol","description", "gene_biotype")
top.deg.ens=getBM(attributes = fields, filters = "ensembl_gene_id", values = top.deg$ensg, mart = ensembl)

fields <- c("ensembl_gene_id", "hgnc_symbol","description", "gene_biotype", "phenotype_description")
top.deg.pheno=getBM(attributes = fields, filters = "ensembl_gene_id", values = top.deg$ensg, mart = ensembl)

fields <- c("ensembl_gene_id", "hgnc_symbol", "go_id", "name_1006", "definition_1006", "namespace_1003")
top.deg.go=getBM(attributes = fields, filters = "ensembl_gene_id", values = top.deg$ensg, mart = ensembl)

dim(top.deg)
dim(top.deg.ens)
dim(top.deg.go)

##
## FG(n=112) and MJ(n=8)
##
fg.top.deg=read.csv(paste0(edgeRHome,"/filtered_toptags_pair_edgeR.csv")
merge(fg.top.deg, top.deg, by.x='ensg', by.y='ensg')

# join two data.frame by a common column
top.deg.anno=merge(top.deg, top.deg.ens, by.x="ensg", by.y="ensembl_gene_id")
length(unique(top.deg.anno$ensg)) # number of unique DEG
write.csv(top.deg.anno[order(top.deg.anno$logFC,decreasing=TRUE),], file=paste0(edgeRHome,'/filtered_toptags_pair_edgeR.',sampleType,'.csv'))

# plot top.deg.anno with the gene name
myCol <- c(`FDR<0.1`="green", `0.1<=FDR<0.2`="orange", `0.2<=FDR<0.3`="purple", `0.3<=FDR<0.4`="red")
hist(top.deg.anno$FDR, breaks=4,col=myCol, xlab='FDR')

# assign color
top.deg.anno$col<-vector(length=nrow(top.deg.anno))
for(i in 1:nrow(top.deg.anno)){
	if(top.deg.anno[i,]$FDR < 0.1){top.deg.anno[i,]$col=myCol[1]}
	if(top.deg.anno[i,]$FDR >= 0.1 & top.deg.anno[i,]$FDR < 0.2 ){top.deg.anno[i,]$col=myCol[2]}
	if(top.deg.anno[i,]$FDR >= 0.2 & top.deg.anno[i,]$FDR < 0.3 ){top.deg.anno[i,]$col=myCol[3]}
	if(top.deg.anno[i,]$FDR >= 0.3 & top.deg.anno[i,]$FDR < 0.4 ){top.deg.anno[i,]$col=myCol[4]}
}

plot(y=top.deg.anno$logFC, x=top.deg.anno$logCPM, cex=0.5, ylab='logFC', xlab='logCPM')
text(y=top.deg.anno$logFC, x=top.deg.anno$logCPM, cex=2, labels=top.deg.anno$hgnc_symbol, col=top.deg.anno$col)
abline(h=c(-minlogFC, minlogFC), col="blue")
legend("topright", legend=names(myCol), fill=myCol, cex=2 )


# scatterplot3D
#scatterplot3d(y=top.deg.anno$FDR, z=top.deg.anno$logFC, x=top.deg.anno$logCPM)

# plot3D
#plot3d(y=top.deg.anno$FDR, z=top.deg.anno$logFC, x=top.deg.anno$logCPM, ylab='FDR', zlab='logFC', xlab='logCPM')
#text3d(y=top.deg.anno$FDR, z=top.deg.anno$logFC, x=top.deg.anno$logCPM, cex=2, top.deg.anno$hgnc_symbol, col=top.deg.anno$col)

# Pheno
top.deg.pheno.anno=merge(top.deg, top.deg.pheno, by.x="ensg", by.y="ensembl_gene_id")
write.csv(top.deg.pheno.anno[order(top.deg.pheno.anno$logFC,decreasing=TRUE),], file=paste0(edgeRHome,'/filtered_pheno_toptags_pair_edgeR.',sampleType,'.csv'))

# GO
top.deg.go.anno=merge(top.deg, top.deg.go, by.x="ensg", by.y="ensembl_gene_id")
write.csv(top.deg.go.anno[order(top.deg.go.anno$logFC,decreasing=TRUE),], file=paste0(edgeRHome,'/filtered_go_toptags_pair_edgeR.',sampleType,'.csv'))

# RData
save(rawCounts,cpms,edgeR.deg,file=paste0(RDataDir, '/edgeR.',sampleType,'.RData')) #  save 'dds'

dev.off()
cat("All is done", "\n")
