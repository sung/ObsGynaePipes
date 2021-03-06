#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 1/May/2014
# Last modified: 1/May/2014 
# Run this script after running tophat/cufflink for all samples 

source /whale-data/ssg29/RNA-Seq/config/rna_seq.config # export envrionment variables
source /whale-data/ssg29/lib/sung.sh #defines 'mkdir_unless'

CASECTR_PAIRS=(R3_R4 R5_R6 R7_R8 R9_R10 R11_R12 R13_R14 
				R15_R16 R17_R18 R19_R20 R21_R22 R23_R24 
				R25_R26 R29_R30 R31_R32 R35_R36 R37_R38 R39_R40 
				R41_R42 R43_R44 R45_R46 R49_R50 R51_R52)

mkdir_unless $RESULT_DIR/$PROJECT/edgeR
##############################################
# 1. make meta-sample info. e.g.: 
# SampleName,Pair,HtseqFile
# R3,R3_R4,/whale-data/ssg29/RNA-Seq/results/Pilot1_v8/HTSeq/R3/R3.igenome.HTSeq.count.txt
# R4,R3_R4,/whale-data/ssg29/RNA-Seq/results/Pilot1_v8/HTSeq/R4/R4.igenome.HTSeq.count.txt
##############################################
META_HTSEQ=$RESULT_DIR/$PROJECT/HTSeq/meta.csv
echo SampleName,Pair,HtseqFile > $META_HTSEQ # header

pair_index=1 # to count pair 
for CASECTR in ${CASECTR_PAIRS[*]}
do
	SAMPLES=${CASECTR/_/ } # R3_R4 => R3 R4

	for Sample in $SAMPLES
	do
		MY_HTSEQ=$RESULT_DIR/$PROJECT/HTSeq/$Sample/$Sample.igenome.HTSeq.count.txt
		if [ -s $MY_HTSEQ ]; then
			echo $Sample,$CASECTR,$MY_HTSEQ >> $META_HTSEQ
			let "sample_index++"
		else
			echo -e "\e[031m$MY_HTSEQ not found. htseq_count failed? \e[0m\n"
			#exit
		fi
	done # end of sample
done # end of case-control pair

echo see $META_HTSEQ

echo "#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# For detial, see http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# simple design 14.A based on the publication above

# pdf output filename 
pdf_file_name<- '$RESULT_DIR/$PROJECT/edgeR/plot_simple'
pdf (file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

# define case and control
samples = read.csv(\"$META_HTSEQ\", stringsAsFactors=FALSE)
samples\$Condition = 'CTR'
samples\$Condition[seq_len(nrow(samples)) %% 2 == 1]='CASE' # odd numbers are cases

# 1. Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count
library("edgeR")
counts = readDGE(samples\$HtseqFile)\$counts

# 2. Filter weakly expressed and noninformative (e.g., non-aligned) features using a command like:
noint = rownames(counts) %in% c(\"no_feature\",\"ambiguous\",\"too_low_aQual\",\"not_aligned\",\"alignment_not_unique\")
cpms = cpm(counts)
keep = rowSums(cpms >1) >=3 & !noint
counts = counts[keep,]

# 3. Visualize and inspect the count table
colnames(counts) = samples\$SampleName
head( counts[,order(samples\$Condition)], 5 ) # print first 5 records

# 4. Create a DGEList object (edgeR's container for RNA-seq count data)
d = DGEList(counts=counts, group=samples\$Condition)

# 5. Estimate normalization factors
d = calcNormFactors(d)

# 6. Inspect the relationships between samples using a multidimensional scaling (MDS) plot
plotMDS(d, labels=samples\$SampleName, col=c(\"darkgreen\",\"blue\")[factor(samples\$Condition)])

# 7. Estimate tagwise dispersion (simple design)
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

# 8. Create a visual representation of the mean-variance relationship using the plotMeanVar (Fig. 5a) and plotBCV (Fig. 5b) functions
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
plotBCV(d)

# 9. Test for differential expression ('classic' edgeR)
de = exactTest(d, pair=c(\"CASE\",\"CTR\"))

# 10. Use the topTags function to present a tabular summary of the differential expression statistics 
tt = topTags(de, n=nrow(d))
head(tt\$table)

# 11. Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt\$table)
head(nc[rn,order(samples\$Condition)],5)

# 12. Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot, 
# here showing the genes selected as differentially expressed (with a 5% false discovery rate; Fig. 6):
deg = rn[tt\$table\$FDR < $FDR ]
plotSmear(d, de.tags=deg)

# 13. Save the result table as a CSV file (alternative formats are possible) as follows
write.csv(tt\$table, file=\"$RESULT_DIR/$PROJECT/edgeR/toptags_simple_edgeR.csv\")

" > $RESULT_DIR/$PROJECT/edgeR/edgeR.simple.design.R
echo R CMD BATCH $RESULT_DIR/$PROJECT/edgeR/edgeR.simple.design.R 
time R CMD BATCH $RESULT_DIR/$PROJECT/edgeR/edgeR.simple.design.R & 

echo "#!/usr/bin/Rscript --vanilla
# Sung Gong <sung@bio.cc>
# For detial, see http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html
# complex design 14.B based on the publication above

# pdf output filename 
pdf_file_name<- '$RESULT_DIR/$PROJECT/edgeR/plot_complex'
pdf (file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

# define case and control
samples = read.csv(\"$META_HTSEQ\", stringsAsFactors=FALSE)
samples\$Condition = 'CTR'
samples\$Condition[seq_len(nrow(samples)) %% 2 == 1]='CASE' # old numbers are case

# 1. Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count
library("edgeR")
counts = readDGE(samples\$HtseqFile)\$counts

# 2. Filter weakly expressed and noninformative (e.g., non-aligned) features using a command like:
noint = rownames(counts) %in% c(\"no_feature\",\"ambiguous\",\"too_low_aQual\",\"not_aligned\",\"alignment_not_unique\")
cpms = cpm(counts)
keep = rowSums(cpms >1) >=3 & !noint
counts = counts[keep,]

# 3. Visualize and inspect the count table
colnames(counts) = samples\$SampleName
head( counts[,order(samples\$Condition)], 5 ) # print first 5 records

# 4. Create a DGEList object (edgeR's container for RNA-seq count data)
d = DGEList(counts=counts, group=samples\$Condition)

# 5. Estimate normalization factors
d = calcNormFactors(d)

# 6. Inspect the relationships between samples using a multidimensional scaling (MDS) plot
plotMDS(d, labels=samples\$SampleName, col=c(\"darkgreen\",\"blue\")[factor(samples\$Condition)])

# 7. Create a design matrix to specify the factors that are expected to affect expression levels:
design <- model.matrix(~Pair+Condition, samples)
design

# 8. Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood
d2 = estimateGLMTrendedDisp(d, design)
d2 = estimateGLMTagwiseDisp(d2, design)

# 9. Given the design matrix and dispersion estimates, fit a GLM to each feature:
f = glmFit(d2, design)

# 10. Perform a likelihood ratio test, specifying the difference of interest 
# (here, knockdown versus control, which corresponds to the third column of the above design matrix):
de = glmLRT(f, coef=XX) # column number for Condition in design data.frame (23 in our case, which is ConditionCTR)

# 11. Use the topTags function to present a tabular summary of the differential expression statistics 
# (note that topTags operates on the output of exactTest or glmLRT, but only the latter is shown here)
tt = topTags(de, n=nrow(d))
head(tt\$table)

# 12. Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt\$table)
head(nc[rn,order(samples\$Condition)],5) # print first 5 entries

# 13. Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot, 
# here showing the genes selected as differentially expressed (with a 5% false discovery rate
deg = rn[tt\$table\$FDR < $FDR ]
plotSmear(d, de.tags=deg)

# 14. Save the result table as a CSV file (alternative formats are possible) as follows:
write.csv(tt\$table, file=\"$RESULT_DIR/$PROJECT/edgeR/toptags_pair_edgeR.csv\")

" > $RESULT_DIR/$PROJECT/edgeR/edgeR.complex.design.R
echo R CMD BATCH $RESULT_DIR/$PROJECT/edgeR/edgeR.complex.design.R 
time R CMD BATCH $RESULT_DIR/$PROJECT/edgeR/edgeR.complex.design.R 
