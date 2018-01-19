 #!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 11/Apr/2014
# Last modified: 1/May/2014 
# assumes 1 case and 1 control

source /whale-data/ssg29/lib/sung.sh #defines 'mkdir_unless'

CASECTR="MY_CASECTR" # case and control (e.g. R3_R4)
SAMPLES="MY_SAMPLES" # case and control (e.g. R3 R4)
CUFFDIFF_LABEL=${SAMPLES/ /,} #'R3 R4'=>'R3,R4'
GENE="MY_GENE" # a gene name to display within R script
CUFFMERG_CASECTR_DIR=$PROJECT_DIR/Cuffmerge/$CASECTR # e.g. Cuffmerge/R3_R4
GTF_ASM=$CUFFMERG_CASECTR_DIR/gtf_assemblies.txt # per-sample 'transcripts.gtf' file
										# e.g. /whale-data/ssg29/RNA-Seq/results/Pilot1/Cuffmerge/R3_R4/gtf_assemblies.txt
CUFFDIFF_CASECTR_DIR=$PROJECT_DIR/Cuffdiff/$CASECTR # e.g. Cuffdiff/R3_R4

mkdir_unless $CUFFMERG_CASECTR_DIR # e.g. /whale-data/ssg29/RNA-Seq/results/Pilot1/Cuffmerge/R3_R4

<<case-control
#case 
/whale-data/dscj1/PE_data/R3.3.fq.gz
/whale-data/dscj1/PE_data/R3.4.fq.gz
/whale-data/dscj1/PE_data/R3.5.fq.gz
#control
/whale-data/dscj1/PE_data/R4.3.fq.gz
/whale-data/dscj1/PE_data/R4.4.fq.gz
/whale-data/dscj1/PE_data/R4.5.fq.gz
case-control

Sample_array=() #an array for $SAMPLES
sample_index=1 # to count sample
# for each sample 
# do 1.tophat and 2.cufflink
for Sample in $SAMPLES
do
	printf "\n\e[32mSample=$Sample\e[0m\n"

	Sample_array+=($Sample) #push it

	FastQ_array=() #an array per sample
	for Lane in $LANES
	do
		##################################################
		# 0. Run the inital fastqc 
		# output: results/Pilot1_v1/FastQC/R3.3.fq_fastqc/
		##################################################
		FastQ_file=$FASTQ_DIR/$Sample.$Lane.fq.gz #e.g. /whale-data/dscj1/PE_data/R3.3.fq.gz
		if [ $RUN_FASTQC -eq 1 ]; then
			if [ -s $FastQ_file ]; then
					echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC"
					time fastqc $FastQ_file -o $PROJECT_DIR/FastQC 
			else
				echo -e "\e[031m$FastQ_file not found\e[0m\n"
				exit
			fi
		fi
		
		############################################################
		# 1. Trim 	
		# output: results/Pilot1_v4/Trim/R3.3_trimmed.fq.gz
		############################################################
		if [ $RUN_TRIM -eq 1 ]; then
			echo -e "\ntrim_galore -a $TR_ADAPTOR1 -o $PROJECT_DIR/Trim --fastqc --quality $TR_QUAL --stringency $TR_STRNCY $FastQ_file"
			time trim_galore --fastqc \
			--output_dir $PROJECT_DIR/Trim \
			--adapter $TR_ADAPTOR1 \
			--quality $TR_QUAL \
			--stringency $TR_STRNCY \
			$FastQ_file

			Trim_FastQ_file=$PROJECT_DIR/Trim/$Sample.$Lane\_trimmed.fq.gz #e.g. results/Pilot1_v8/Trim/R3.3_trimmed.fq.gz
			if [ -s $Trim_FastQ_file ]; then
				FastQ_array+=($Trim_FastQ_file) #push it
			else
				echo -e "\e[031m$Trim_FastQ_file not found\e[0m\n"
				exit
			fi
		else
			FastQ_array+=($FastQ_file) #push it
		fi

	done # end of per lane

	####################################################################
	# 2. tophat: Map the reads for each sample to the reference genome #
	# real    61m11.381s (for R3 with 2 threads) 82m22.024s (for R4)
	# output: $TOPHAT_OUT/accepted_hits.bam
	####################################################################
	TOPHAT_OUT=$PROJECT_DIR/TopHat/$Sample
	if [ $RUN_TOPHAT -eq 1 ]; then
		mkdir_unless $TOPHAT_OUT
		Merged_FastQ=$(printf ",%s" "${FastQ_array[@]}")
		Merged_FastQ=${Merged_FastQ:1} # /whale-data/dscj1/PE_data/R3.3.fq.gz,/whale-data/dscj1/PE_data/R3.4.fq.gz,/whale-data/dscj1/PE_data/R3.5.fq.gz
		echo -e "\ntophat2 -p $NT --library-type $LIB_TYPE --output-dir $TOPHAT_OUT --max-multihits $TH_MH $TH_PREFILTER -transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX $BOWTIE2_INDEX_BASE $Merged_FastQ"

		time tophat2 -p $NT \
		--library-type $LIB_TYPE \
		--output-dir $TOPHAT_OUT \
		--max-multihits $TH_MH \
		$TH_PREFILTER \
		--transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX \
		$BOWTIE2_INDEX_BASE $Merged_FastQ

		cuffdiff_input+="$TOPHAT_OUT/accepted_hits.bam " #e.g. /whale-data/ssg29/RNA-Seq/results/Pilot1/TopHat/R3/accepted_hits.bam 
	fi

	######################
	# CoverageBed
	# Input: bam from tophat
	# Output: $sample.genomecov.txt
	######################
	if [ $RUN_GENOMECOV -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Coverage/$Sample
		if [ -s $TOPHAT_OUT/accepted_hits.bam ]; then
			echo -e "bedtools genomecov $TOPHAT_OUT/accepted_hits.bam";
			#time bedtools genomecov -ibam $TOPHAT_OUT/accepted_hits.bam > $PROJECT_DIR/Coverage/$Sample/$Sample.genomecov.txt
			time bedtools genomecov -split -ibam $TOPHAT_OUT/accepted_hits.bam > $PROJECT_DIR/Coverage/$Sample/$Sample.genomecov.split.txt
		else
			echo -e "\e[031m$TOPHAT_OUT/accepted_hits.bam not found. tophat failed? \e[0m\n"
		fi
	fi

	######################
	# 3. HTSeq
	# For paired-end data, the alignment have to be sorted either by read name or by alignment position.
	# 'sort -T $SCRATCH_OUT -k 1,1' will do a sort by read name 
	# computation time?
	# output:
	######################
	HTSEQ_OUT=$PROJECT_DIR/HTSeq/$Sample
	if [ $RUN_HTSEQ -eq 1 ]; then
		SCRATCH_OUT=$PROJECT_DIR/scratch/$Sample
		mkdir_unless $HTSEQ_OUT
		mkdir_unless $SCRATCH_OUT

		# count based on latest Ensembl
		echo -e "\nsamtools view $TOPHAT_OUT/accepted_hits.bam | sort -T $SCRATCH_OUT -k 1,1 | htseq-count --stranded=$HTSEQ_STRAND --type=$HTSEQ_TYPE --idattr=$HTSEQ_IDATTR - $GTF2 > $HTSEQ_OUT/$Sample.GRCh37.75.HTSeq.count.txt"
		time samtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count \
			--stranded=$HTSEQ_STRAND \
			--type=$HTSEQ_TYPE \
			--quiet \
			--idattr=$HTSEQ_IDATTR \
			- $GTF2 > $HTSEQ_OUT/$Sample.GRCh37.75.HTSeq.count.txt

		# count based on illumina iGenome
		echo -e "\nsamtools view $TOPHAT_OUT/accepted_hits.bam | sort -T $SCRATCH_OUT -k 1,1 | htseq-count --stranded=$HTSEQ_STRAND --type=$HTSEQ_TYPE --quiet --idattr=$HTSEQ_IDATTR - $GTF > $HTSEQ_OUT/$Sample.igenome.HTSeq.count.txt"
		time samtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count \
			--stranded=$HTSEQ_STRAND \
			--type=$HTSEQ_TYPE \
			--quiet \
			--idattr=$HTSEQ_IDATTR \
			- $GTF > $HTSEQ_OUT/$Sample.igenome.HTSeq.count.txt
	fi

	#####################################################
	# 4. cufflinks: Assemble transcripts for each sample #
	# The SAM file supplied to Cufflinks must be sorted by reference position. 
	# If you aligned your reads with TopHat, your alignments will be properly sorted already
	# o/w, sort -k 3,3 -k 4,4n hits.sam > hits.sam.sorted
	# computation time? real    352m42.053s (for R4)
	# output: $PROJECT_DIR/Cufflink/$Sample/transcripts.gtf 
	#####################################################
	CUFFLINK_OUT=$PROJECT_DIR/Cufflink/$Sample
	# cufflinks will make the output dir if not exist 
	if [ $RUN_CUFFLINK -eq 1 ]; then
		if [ -s $TOPHAT_OUT/accepted_hits.bam ]; then
				# --GTF: reference transcript only (http://cufflinks.cbcb.umd.edu/manual.html)
				# --GTF-guide: reference transcript + novel (http://cufflinks.cbcb.umd.edu/howitworks#hrga)
				echo "cufflinks -p $NT --library-type $LIB_TYPE --GTF-guide $GTF --output-dir $CUFFLINK_OUT $TOPHAT_OUT/accepted_hits.bam"
				time cufflinks -p $NT --library-type $LIB_TYPE --GTF-guide $GTF --output-dir $CUFFLINK_OUT $TOPHAT_OUT/accepted_hits.bam
		else
			echo -e "\e[031m$TOPHAT_OUT/accepted_hits.bam not found. tophat failed? \e[0m\n"
			exit
		fi
	fi

	###############################
	# 5. write a gtf assembly file #
	# output: $CUFFMERG_CASECTR_DIR/gtf_assemblies.txt, which is $GTF_ASM
	###############################
	if [ $RUN_CUFFMERGE -eq 1 ]; then
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
	fi

	echo -e "\e[32m$Sample done\e[0m" 
	let "sample_index++"
done # end of per sample

#########################################################################################################
# 6. cuffmerge: Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation: #
# output: $CUFFMERG_CASECTR_DIR/merged.gtf
# running time: 291m (5h) for R3-R4
# running time: 43m51.700s for R3-R4
#########################################################################################################
if [ $RUN_CUFFMERGE -eq 1 ]; then
	if [ -s $GTF_ASM ]; then
			echo -e "\ncuffmerge -o $CUFFMERG_CASECTR_DIR -g $GTF -s $GENOME -p $NT $GTF_ASM"
			time cuffmerge -o $CUFFMERG_CASECTR_DIR -g $GTF -s $GENOME -p $NT $GTF_ASM
	else
		echo -e "\e[031m$GTF_ASM not found. cufflinks failed? \e[0m\n"
		exit
	fi
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
		echo "cuffdiff -o $CUFFDIFF_CASECTR_DIR --library-type $LIB_TYPE -b $GENOME -p $NT -L $CUFFDIFF_LABEL -u $CUFFMERG_CASECTR_DIR/merged.gtf $cuffdiff_input"
		time cuffdiff -o $CUFFDIFF_CASECTR_DIR --library-type $LIB_TYPE -b $GENOME -p $NT -L $CUFFDIFF_LABEL -u $CUFFMERG_CASECTR_DIR/merged.gtf $cuffdiff_input &
	else
		echo -e "\e[031m$CUFFMERG_CASECTR_DIR/merged.gtf not found. cuffmerge failed? \e[0m\n"
		exit
	fi
fi

#################################
# 8. cuffnorm per case-control
# computing time: real    46m14.464s (e.g. R3_R4)
# output: Cuffnorm/R3_R4/genes.fpkm_table, Cuffnorm/R3_R4/genes.count_table...
#################################
if [ $RUN_CUFFNORM -eq 1 ]; then
	if [ -s $CUFFMERG_CASECTR_DIR/merged.gtf ]; then
		echo -e "\ncuffnorm -o $PROJECT_DIR/Cuffnorm/$CASECTR -p $NT -L $CUFFDIFF_LABEL $CUFFMERG_CASECTR_DIR/merged.gtf $cuffdiff_input"
		time cuffnorm -o $PROJECT_DIR/Cuffnorm/$CASECTR -p $NT -L $CUFFDIFF_LABEL $CUFFMERG_CASECTR_DIR/merged.gtf $cuffdiff_input 
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
	csScatter(genes(cuff_data), '${Sample_array[0]}', '${Sample_array[1]}')

	# Create a volcano plot to inspect differentially expressed genes
	csVolcano(genes(cuff_data), '${Sample_array[0]}', '${Sample_array[1]}')

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
	isoform_diff_data <- diffData(isoforms(cuff_data), '${Sample_array[0]}', '${Sample_array[1]}')
	sig_isoform_data <- subset(isoform_diff_data, (significant == 'yes'))
	write.table(sig_isoform_data, '$CUFFDIFF_CASECTR_DIR/diff_isoform.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# tss
	tss_diff_data <- diffData(TSS(cuff_data), '${Sample_array[0]}', '${Sample_array[1]}')
	sig_tss_data <- subset(tss_diff_data, (significant == 'yes'))
	write.table(sig_tss_data, '$CUFFDIFF_CASECTR_DIR/diff_tss.txt', sep='\t', row.names = F, col.names = T, quote = F)

	# cds
	cds_diff_data <- diffData(CDS(cuff_data), '${Sample_array[0]}', '${Sample_array[1]}')
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
