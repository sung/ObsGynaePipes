#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 11/Apr/2014
# Last modified: 27/Jun/2014 
# assumes 1 case and 1 control separated by ',' (e.g. D709_D501,D709_D502)

source /whale-data/ssg29/lib/sung.sh #defines 'mkdir_unless'

CASECTR="MY_CASECTR" # case and control (e.g. D709_D501,D709_D502)
BARCODES=(MY_BARCODES) # case and control (e.g. D709_D501 D709_D502)
GENE="MY_GENE" # a gene name to display within R script

CUFFMERG_CASECTR_DIR=$PROJECT_DIR/Cuffmerge/$CASECTR # e.g. Cuffmerge/D709_D501,D709_D502
GTF_ASM=$CUFFMERG_CASECTR_DIR/gtf_assemblies.txt # per-sample 'transcripts.gtf' file
										# e.g. /whale-data/ssg29/RNA-Seq/results/Pilot1/Cuffmerge/R3_R4/gtf_assemblies.txt
CUFFDIFF_CASECTR_DIR=$PROJECT_DIR/Cuffdiff/$CASECTR # e.g. Cuffdiff/D709_D501,D709_D502

mkdir_unless $CUFFMERG_CASECTR_DIR

sample_index=1 # to count sample

echo -e `hostname`
echo -e "NO. of thread=$NT"

<<fastq
/whale-data/ssg29/data/SLX-8546/SLX-8546.D709_D501.C4FN0ANXX.s_4.r_1.fq.gz 
/whale-data/ssg29/data/SLX-8546/SLX-8546.D709_D501.C4FN0ANXX.s_5.r_1.fq.gz 
/whale-data/ssg29/data/SLX-8546/SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1.fq.gz 
fastq

# assuming single-end
# [todo] what if paired end?
if [ $IS_SE -ne 1 ]; then
	echo -e "\e[031This script is for signle-end only\e[0m\n"
	exit
fi

# per each barcode (sample) 
# do 1.FastQC 2.trim 3.tophat 4.cufflink 5.genomcov 6.htseq and 7.cufflink
for Barcode in ${BARCODES[*]} # case/ctrl as a pair 
do
	printf "\n\e[32mBarcode=$Barcode\e[0m\n"

	Trimmed_FastQ_array=() #initialise
	FASTQ_FILES=`ls $FASTQ_DIR/$SLX.$Barcode*.fq.gz | grep -v lost`
	for FastQ_file in $FASTQ_FILES	# per fastq file
	do
		if [ ! -s $FastQ_file ]; then
			echo -e "\e[031m$FastQ_file not found\e[0m\n"
			exit
		fi
		##################################################
		# 1. Run the inital fastqc 
		##################################################
		if [ $RUN_FASTQC -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/FastQC/$Barcode
			echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
			time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
		fi
		############################################################
		# 2. Trim 	
		############################################################
		if [ $RUN_TRIM -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/Trim/$Barcode
			echo -e "\ntrim_galore -a $TR_ADAPTOR1 -o $PROJECT_DIR/Trim/$Barcode --fastqc --quality $TR_QUAL --stringency $TR_STRNCY $FastQ_file"
			time trim_galore \
				--fastqc \
				--output_dir $PROJECT_DIR/Trim/$Barcode \
				--adapter $TR_ADAPTOR1 \
				--quality $TR_QUAL \
				--stringency $TR_STRNCY \
				$FastQ_file
		fi

		#/whale-data/ssg29/RNA-Seq/data/SLX-8546/SLX-8546.D709_D501.C4FN0ANXX.s_4.r_1.fq.gz
		#/whale-data/ssg29/RNA-Seq/results/SLX-8546.v1/Trim/D709_D501/SLX-8546.D709_D501.C4FN0ANXX.s_4.r_1_trimmed.fq.gz
		Trimmed_FastQ=${FastQ_file/data\/$SLX/results\/$PROJECT\/Trim\/$Barcode}
		Trimmed_FastQ=${Trimmed_FastQ%.fq.gz}_trimmed.fq.gz #SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1_trimmed.fq.gz 
		if [ -s $Trimmed_FastQ ]; then
			Trimmed_FastQ_array+=($Trimmed_FastQ) #push it
		else
			echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
			exit
		fi
	done # end of per FastQ_file 

	####################################################################
	# 2. tophat: Map the reads for each sample to the reference genome #
	# real    61m11.381s (for R3 with 2 threads) 82m22.024s (for R4)   #
	# output: $TOPHAT_OUT/accepted_hits.bam                            #
	####################################################################
	TOPHAT_OUT=$PROJECT_DIR/TopHat/$Barcode
	if [ $RUN_TOPHAT -eq 1 ]; then
		mkdir_unless $TOPHAT_OUT
		Merged_FastQ=$(printf ",%s" "${Trimmed_FastQ_array[@]}")
		Merged_FastQ=${Merged_FastQ:1} # /whale-data/dscj1/PE_data/R3.3.fq.gz,/whale-data/dscj1/PE_data/R3.4.fq.gz,/whale-data/dscj1/PE_data/R3.5.fq.gz
		echo -e "\ntophat2 -p $NT --library-type $LIB_TYPE --output-dir $TOPHAT_OUT --max-multihits $TH_MH $TH_PREFILTER -transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX $BOWTIE2_INDEX_BASE $Merged_FastQ"

		time tophat2 \
			--num-threads $NT \
			--library-type $LIB_TYPE \
			--output-dir $TOPHAT_OUT \
			--max-multihits $TH_MH \
			$TH_PREFILTER \
			--transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX \
			$BOWTIE2_INDEX_BASE $Merged_FastQ
	fi

	# check out result
	if [ -s $TOPHAT_OUT/accepted_hits.bam ];then #e.g. /whale-data/ssg29/RNA-Seq/results/Pilot1/TopHat/R3/accepted_hits.bam 
		cuffdiff_input+="$TOPHAT_OUT/accepted_hits.bam " #e.g. /whale-data/ssg29/RNA-Seq/results/Pilot1/TopHat/R3/accepted_hits.bam 
	else
		echo -e "\e[031m$TOPHAT_OUT/accepted_hits.bam not found. tophat failed? \e[0m\n"
		exit
	fi

	######################
	# 3. CoverageBed
	# Input: bam from tophat
	# Output: $sample.genomecov.txt
	######################
	if [ $RUN_GENOMECOV -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Coverage/$Barcode
		echo -e "bedtools genomecov $TOPHAT_OUT/accepted_hits.bam";
		time bedtools genomecov -split -ibam $TOPHAT_OUT/accepted_hits.bam > $PROJECT_DIR/Coverage/$Barcode/$Barcode.genomecov.split.txt
	fi

	######################
	# 4. HTSeq
	# For paired-end data, the alignment have to be sorted either by read name or by alignment position.
	# 'sort -T $SCRATCH_OUT -k 1,1' will do a sort by read name, or
	# samtools sort -n $TOPHAT_OUT/accepted_hits.bam TOPHAT_OUT/accepted_hits.nsorted.bam
	# output:
	######################
	HTSEQ_OUT=$PROJECT_DIR/HTSeq/$Barcode
	SCRATCH_OUT=$PROJECT_DIR/scratch/$Barcode
	if [ $RUN_HTSEQ -eq 1 ]; then
		mkdir_unless $HTSEQ_OUT
		mkdir_unless $SCRATCH_OUT

		# count based on the latest Ensembl GTF
		echo -e "\nsamtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count --stranded=$HTSEQ_STRAND --type=$HTSEQ_TYPE --idattr=$HTSEQ_IDATTR - $GTF2 > $HTSEQ_OUT/$Barcode.GRCh37.75.HTSeq.count.txt"
		time samtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count \
			--stranded=$HTSEQ_STRAND \
			--type=$HTSEQ_TYPE \
			--quiet \
			--idattr=$HTSEQ_IDATTR \
			- $GTF2 > $HTSEQ_OUT/$Barcode.GRCh37.75.HTSeq.count.txt

		# count based on illumina iGenome GTF
		echo -e "\nsamtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count --stranded=$HTSEQ_STRAND --type=$HTSEQ_TYPE --quiet --idattr=$HTSEQ_IDATTR - $GTF > $HTSEQ_OUT/$Barcode.igenome.HTSeq.count.txt"
		time samtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count \
			--stranded=$HTSEQ_STRAND \
			--type=$HTSEQ_TYPE \
			--quiet \
			--idattr=$HTSEQ_IDATTR \
			- $GTF > $HTSEQ_OUT/$Barcode.igenome.HTSeq.count.txt
	fi

	#####################################################
	# 5. cufflinks: Assemble transcripts for each sample #
	# The SAM file supplied to Cufflinks must be sorted by reference position. 
	# If you aligned your reads with TopHat, your alignments will be properly sorted already
	# o/w, sort -k 3,3 -k 4,4n hits.sam > hits.sam.sorted
	# computation time? real    352m42.053s (for R4)
	# output: $PROJECT_DIR/Cufflink/$Barcode/transcripts.gtf 
	#####################################################
	CUFFLINK_OUT=$PROJECT_DIR/Cufflink/$Barcode
	# cufflinks will make the output dir if not exist 
	if [ $RUN_CUFFLINK -eq 1 ]; then
		# --GTF: reference transcript only (http://cufflinks.cbcb.umd.edu/manual.html)
		# --GTF-guide: reference transcript + novel (http://cufflinks.cbcb.umd.edu/howitworks#hrga)
		echo "cufflinks -p $NT --library-type $LIB_TYPE --GTF-guide $GTF --output-dir $CUFFLINK_OUT $TOPHAT_OUT/accepted_hits.bam"
		time cufflinks \
			--num-threads $NT \
			--library-type $LIB_TYPE \
			--GTF-guide $GTF \
			--frag-bias-correct $GENOME \
			--multi-read-correct \
			--output-dir $CUFFLINK_OUT \
			$TOPHAT_OUT/accepted_hits.bam
			#--library-norm-method classic-fpkm # default
	fi

	###############################
	# 6. write a gtf assembly file #
	# to run cuffmerge later
	# output: $CUFFMERG_CASECTR_DIR/gtf_assemblies.txt, which is $GTF_ASM
	###############################
	if [ -s $CUFFLINK_OUT/transcripts.gtf ]; then
			if [ $sample_index -eq 1 ]; then
				echo $CUFFLINK_OUT/transcripts.gtf > $GTF_ASM
			else
				echo $CUFFLINK_OUT/transcripts.gtf >> $GTF_ASM
			fi
	else
		echo -e "\e[031m$CUFFLINK_OUT/transcripts.gtf not found. cufflinks failed? \e[0m\n"
		exit
	fi

	echo -e "\e[32m$Barcode done\e[0m" 
	let "sample_index++"
done # end of per sample

#########################################################################################################
# 6. cuffmerge: Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation: #
# per case-ctr
# output: $CUFFMERG_CASECTR_DIR/merged.gtf
# running time: 291m (5h) for R3-R4
# running time: 43m51.700s for R3-R4
#########################################################################################################
if [ $RUN_CUFFMERGE -eq 1 ]; then
	if [ -s $GTF_ASM ]; then
		echo -e "\ncuffmerge -o $CUFFMERG_CASECTR_DIR -g $GTF -s $GENOME -p $NT $GTF_ASM"
		time cuffmerge \
			--num-threads $NT \
			--ref-sequence $GENOME \
			--ref-gtf $GTF \
			-o $CUFFMERG_CASECTR_DIR \
			$GTF_ASM
	else
		echo -e "\e[031m$GTF_ASM not found. cufflinks failed? \e[0m\n"
		exit
	fi
fi

if [ $RUN_CUFFQUANT -eq 1 ]; then
	time cuffquant \

fi

#################################################################################################################################
# 7. cuffdiff: Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate: #
# output: *.diff, *.fpkm_tracking, *.count_tracking, *.read_group_tracking, read_groups.info, run.info and bias_params.info
# Q: which dispersion model to use?
# Warning: No conditions are replicated, switching to 'blind' dispersion method
# http://cufflinks.cbcb.umd.edu/manual.html#library_norm_meth
# running time: 302m (5h) for R3-R4
# running time: 1226m8.763s for R3-R4
################################################################################################################################
if [ $RUN_CUFFDIFF -eq 1 ]; then
	if [ -s $CUFFMERG_CASECTR_DIR/merged.gtf ]; then
		echo -e "cuffdiff -o $CUFFDIFF_CASECTR_DIR --library-type $LIB_TYPE -b $GENOME -p $NT -L $CASECTR -u $CUFFMERG_CASECTR_DIR/merged.gtf $cuffdiff_input"
		time cuffdiff \
			--num-threads $NT \
			--output-dir $CUFFDIFF_CASECTR_DIR \
			--library-type $LIB_TYPE \
			--labels $CASECTR \
			--frag-bias-correct $GENOME \
			--multi-read-correct \
			$CUFFMERG_CASECTR_DIR/merged.gtf \
			$cuffdiff_input &
			#--FDR 0.05 \ # default
			#--library-norm-method geometric \ # default
			#--dispersion-method pooled \ # default
	else
		echo -e "\e[031m$CUFFMERG_CASECTR_DIR/merged.gtf not found. cuffmerge failed? \e[0m\n"
		exit
	fi
fi

#################################
# 8. cuffnorm 
# per case-control
# computing time: real    46m14.464s (e.g. R3_R4)
# output: Cuffnorm/R3_R4/genes.fpkm_table, Cuffnorm/R3_R4/genes.count_table...
#################################
if [ $RUN_CUFFNORM -eq 1 ]; then
	if [ -s $CUFFMERG_CASECTR_DIR/merged.gtf ]; then
		echo -e "\ncuffnorm -o $PROJECT_DIR/Cuffnorm/$CASECTR -p $NT -L $CASECTR $CUFFMERG_CASECTR_DIR/merged.gtf $cuffdiff_input"
		time cuffnorm \
			-p $NT \
			-L $CASECTR \
			-o $PROJECT_DIR/Cuffnorm/$CASECTR \
			$CUFFMERG_CASECTR_DIR/merged.gtf \
			$cuffdiff_input 
	else
		echo -e "\e[031m$CUFFMERG_CASECTR_DIR/merged.gtf not found. cuffmerge failed? \e[0m\n"
		exit
	fi
fi

###############################################################
# 9. cummeRbund R script
# running time: real    17m37.568s
###############################################################
if [ $RUN_CUMMERBUND -eq 1 ]; then
	echo "#!/usr/bin/Rscript --vanilla
	# Sung Gong <sung@bio.cc>
	# see http://www.nature.com/nprot/journal/v7/n3/pdf/nprot.2012.016.pdf
	# install unless there (install failed on butterfly--bio)
	# source('http://www.bioconductor.org/biocLite.R')
	# biocLite('cummeRbund')

	library(cummeRbund)

	# output filename 
	pdf_file_name<- '$CUFFDIFF_CASECTR_DIR/diff'
	pdf (file=paste(pdf_file_name,format(Sys.time(), '%Y-%m-%d_%I%p'), 'pdf', sep ='.'))

	cuff_data <- readCufflinks('$CUFFDIFF_CASECTR_DIR')

	# Plot the distribution of expression levels for each sample
	csDensity(genes(cuff_data))

	# Compare the expression of each gene in two conditions with a scatter plot 
	csScatter(genes(cuff_data), '${BARCODES[0]}', '${BARCODES[1]}')

	# Create a volcano plot to inspect differentially expressed genes
	csVolcano(genes(cuff_data), '${BARCODES[0]}', '${BARCODES[1]}')

	# Plot expression levels for genes of interest with bar plots
	mygene<- getGene(cuff_data,'$GENE')
	expressionBarplot(mygene)

	# Plot individual isoform expression levels of selected genes of interest with bar plots
	expressionBarplot(isoforms(mygene))

	# quickly inspect the number of genes and transcripts that are differentially expressed between two samples
	gene_diff_data <- diffData(genes(cuff_data))
	sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
	write.table(sig_gene_data, '$CUFFDIFF_CASECTR_DIR/diff_genes.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# isoform
	isoform_diff_data <- diffData(isoforms(cuff_data), '${BARCODES[0]}', '${BARCODES[1]}')
	sig_isoform_data <- subset(isoform_diff_data, (significant == 'yes'))
	write.table(sig_isoform_data, '$CUFFDIFF_CASECTR_DIR/diff_isoform.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# tss
	tss_diff_data <- diffData(TSS(cuff_data), '${BARCODES[0]}', '${BARCODES[1]}')
	sig_tss_data <- subset(tss_diff_data, (significant == 'yes'))
	write.table(sig_tss_data, '$CUFFDIFF_CASECTR_DIR/diff_tss.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# cds
	cds_diff_data <- diffData(CDS(cuff_data), '${BARCODES[0]}', '${BARCODES[1]}')
	sig_cds_data <- subset(cds_diff_data, (significant == 'yes'))
	write.table(sig_cds_data, '$CUFFDIFF_CASECTR_DIR/diff_cds.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# promoter
	promoter_diff_data <- distValues(promoters(cuff_data))
	sig_promoter_data <- subset(promoter_diff_data, (significant == 'yes'))
	write.table(sig_promoter_data, '$CUFFDIFF_CASECTR_DIR/diff_promoter.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# splicing
	splicing_diff_data <- distValues(splicing(cuff_data))
	sig_splicing_data <- subset(splicing_diff_data, (significant == 'yes'))
	write.table(sig_splicing_data, '$CUFFDIFF_CASECTR_DIR/diff_splicing.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# relCDS
	relCDS_diff_data <- distValues(relCDS(cuff_data))
	sig_relCDS_data <- subset(relCDS_diff_data, (significant == 'yes'))
	write.table(sig_relCDS_data, '$CUFFDIFF_CASECTR_DIR/diff_relCDS.txt', sep='\t', row.names = F, col.names = T, quote = F)

	" > $CUFFDIFF_CASECTR_DIR/cummeRbund.script.R

	cd $CUFFDIFF_CASECTR_DIR
	echo R CMD BATCH cummeRbund.script.R
	time R CMD BATCH cummeRbund.script.R
fi
