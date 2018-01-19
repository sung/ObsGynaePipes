#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# For detial, see http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# complex design 14.B based on the publication above

# pdf output filename 
pdf_file_name<- '/home/ssg29/scratch/results/SLX-8549.v1/edgeR/8oxBS.plot.complex'
pdf (file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

# define case and control
samples = read.csv("/home/ssg29/scratch/results/SLX-8549.v1/HTSeq/meta.oxBS.csv", stringsAsFactors=FALSE)
samples$Condition = 'CTR'
samples$Condition[seq_len(nrow(samples)) %% 2 == 1]='CASE' # odd numbers are cases (case,control)

# 1. Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count
library(edgeR)
counts = readDGE(samples$HtseqFile)$counts

# 2. Filter weakly expressed and noninformative (e.g., non-aligned) features using a command like:
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms >1) >=4 & !noint
counts = counts[keep,]

# 3. Visualize and inspect the count table
colnames(counts) = samples$SampleName
head( counts[,order(samples$Condition)], 5 ) # print first 5 records
summary(counts)

# 4. Create a DGEList object (edgeR's container for RNA-seq count data)
d = DGEList(counts=counts, group=samples$Condition)

# 5. Estimate normalization factors
# The calcNormFactors function normalizes for RNA composition by finding a set of scaling
# factors for the library sizes that minimize the log-fold changes between the samples for most genes
# The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples 
d = calcNormFactors(d)
dim(d)

# 6. Inspect the relationships between samples using a multidimensional scaling (MDS) plot
plotMDS(d, labels=samples$SampleName, col=c("darkgreen","blue")[factor(samples$Condition)])

# 7. Create a design matrix to specify the factors that are expected to affect expression levels:
design <- model.matrix(~ Pair + Condition, samples)
design

# 8. Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood
d2 = estimateGLMCommonDisp(d, design, verbose=TRUE)
d2 = estimateGLMTrendedDisp(d2, design, verbose=T)
d2 = estimateGLMTagwiseDisp(d2, design)

plotMeanVar(d2, show.tagwise.vars=TRUE, NBline=TRUE)
plotBCV(d2)

# 9. Given the design matrix and dispersion estimates, fit a GLM to each feature:
f = glmFit(d2, design)

# 10. Perform a likelihood ratio test, specifying the difference of interest 
# (here, knockdown versus control, which corresponds to the third column of the above design matrix):
de = glmLRT(f) # the column number for Condition, which is ConditionCTR in our case, in design data.frame

o = order(de$table$PValue)
cpm(d2)[o[1:10],]

# 11. Use the topTags function to present a tabular summary of the differential expression statistics 
# (note that topTags operates on the output of exactTest or glmLRT, but only the latter is shown here)
tt = topTags(de, n=nrow(d2)) # default sort by pvalue?
#tt2 = topTags(de, n=nrow(d), sort.by="logFC")
head(tt$table)

# 12. Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:
nc = cpm(d2, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn,order(samples$Condition)],5) # print first 5 entries

# 13. Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot, 
# here showing the genes selected as differentially expressed (with a 5% false discovery rate
deg = rn[tt$table$FDR < .15 ]
plotSmear(de, de.tags=deg)
abline(h=c(-1, 1), col="blue")

deg = rn[tt$table$PValue < 1.0e-02 ]
plotSmear(de, de.tags=deg)
abline(h=c(-1, 1), col="blue")

deg = rn[(tt$table$logFC < -1 | tt$table$logFC > 1) & tt$table$PValue < 1.0e-02 ]
plotSmear(de, de.tags=deg)
abline(h=c(-1, 1), col="blue")

# 14. Save the result table as a CSV file (alternative formats are possible) as follows:
write.csv(tt$table, file="/home/ssg29/scratch/results/SLX-8549.v1/edgeR/8oxBS.toptags_pair_edgeR.csv")


