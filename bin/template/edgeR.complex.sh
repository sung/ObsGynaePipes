#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 27/Aug/2014
# Last modified: 27/Aug/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT

mkdir_unless $PROJECT_DIR	
mkdir_unless $RESULT_DIR/$PROJECT/edgeR

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mProject=$PROJECT \e[0m\n"

##############################################
# 1. make meta-sample info. e.g.: 
# SampleName,Pair,HtseqFile
# R3,R3_R4,/whale-data/ssg29/RNA-Seq/results/Pilot1_v8/HTSeq/R3/R3.igenome.HTSeq.count.txt
# R4,R3_R4,/whale-data/ssg29/RNA-Seq/results/Pilot1_v8/HTSeq/R4/R4.igenome.HTSeq.count.txt
##############################################
META_HTSEQ=$RESULT_DIR/$PROJECT/HTSeq/meta.csv
echo SampleName,Pair,HtseqFile > $META_HTSEQ # header

COEF=$((${#CASECTR_PAIRS[@]}+1))
for CASECTR in ${CASECTR_PAIRS[*]}
do
	BARCODES=${CASECTR/,/ } # D709_D501,D709_D502 => D709_D501 D709_D502 

	for Sample in $BARCODES
	do
		MY_HTSEQ=$RESULT_DIR/$PROJECT/HTSeq/$Sample/$SLX.$Sample.igenome.HTSeq.count.txt
		if [ -s $MY_HTSEQ ]; then
			echo $Sample,${CASECTR/,/|},$MY_HTSEQ >> $META_HTSEQ
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
# complex design 14.B based on the publication above

# pdf output filename 
pdf_file_name<- '$RESULT_DIR/$PROJECT/edgeR/plot_complex'
pdf (file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

# define case and control
samples = read.csv(\"$META_HTSEQ\", stringsAsFactors=FALSE)
samples\$Condition = 'CTR'
samples\$Condition[seq_len(nrow(samples)) %% 2 == 0]='CASE' # even numbers are cases

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
de = glmLRT(f, coef=$COEF) # the column number for Condition, which is ConditionCTR in our case, in design data.frame

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

echo -e "Rscript $RESULT_DIR/$PROJECT/edgeR/edgeR.complex.design.R "
time Rscript $RESULT_DIR/$PROJECT/edgeR/edgeR.complex.design.R 
