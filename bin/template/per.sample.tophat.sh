#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Jul/2014
# Last modified: 28/Mar/2017
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh # defines user defined functions (e.g. make_run_script)
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir -p $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"

<<fastq
SLX-9169.D704_D501.C6A20ANXX.s_1.r_1.fq.gz
SLX-9169.D704_D501.C6GNPANXX.s_5.r_1.fq.gz
SLX-9169.D704_D501.C6GNPANXX.s_6.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_1.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_2.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_3.r_1.fq.gz
fastq

if [ $IS_SE -eq 1 ]; then
# per each barcode (sample) 
# do 1)FastQC 2)trim 3)tophat 4)bedtool, 5)htseq 
	Trimmed_FastQ_array=() #initialise
	FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`
	for FastQ_file in $FASTQ_FILES	# per fastq file (per lane)
	do
		##################################################
		# 0. Run the inital fastqc 
		##################################################
		if [ $RUN_FASTQC -eq 1 ]; then
			mkdir -p $PROJECT_DIR/FastQC/$Barcode
			echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
			time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
		fi
		############################################################
		# 1. Trimming via trim_galore 	
		############################################################
		#/home/ssg29/data/fastq/SLX-9176/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1_trimmed.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz_trimming_report.txt
		Trimmed_FastQ=${FastQ_file/data\/fastq\/$SLX/results\/$PROJECT\/Trim\/$Barcode}
		Trimmed_FastQ=${Trimmed_FastQ%.fq.gz}_trimmed.fq.gz #SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1_trimmed.fq.gz 
		if [ $RUN_TRIM -eq 1 ]; then
			mkdir -p $PROJECT_DIR/Trim/$Barcode
			echo -e "\ntrim_galore -a $TR_ADAPTOR1 -o $PROJECT_DIR/Trim/$Barcode --fastqc --quality $TR_QUAL --stringency $TR_STRNCY $FastQ_file"
			time trim_galore \
			--output_dir $PROJECT_DIR/Trim/$Barcode \
			--adapter $TR_ADAPTOR1 \
			--quality $TR_QUAL \
			--stringency $TR_STRNCY \
			$FastQ_file
			#--fastqc \
		fi # end of run trim

		# Plasma samples prepared by SMARTer kit of clonetech
		if [ $RUN_TRIM_SMARTER -eq 1 ]; then
			mkdir -p $PROJECT_DIR/Trim/$Barcode
			echo -e "\ncutadapt $FastQ_file"
			time cutadapt \
				-q $TR_QUAL \
				-O $TR_STRNCY \
				-m $TR_LEN \
				-a $TR_ADAPTOR1 \
				-u $TR_CUT_LEN \
				-o $Trimmed_FastQ \
				$FastQ_file &> ${Trimmed_FastQ%_trimmed.fq.gz}.fq.gz_trimming_report.txt
			if [ ! -s $Trimmed_FastQ ];then
				echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
				exit
			else
				echo -e "\nfastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode"
				#time fastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode 
			fi
		fi # end of run trim

		# scrap trimmed reads
		if [ -s $Trimmed_FastQ ]; then
			Trimmed_FastQ_array+=($Trimmed_FastQ) #push it
		else
			echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
			exit
		fi
	done # end of per FastQ_file 

	###################################################################
	# 2. TOPHAT 
	# Map the reads for each sample to the reference genome           #
	# real    61m11.381s (for R3 with 2 threads) 82m22.024s (for R4)
	# output: $TOPHAT_OUT/accepted_hits.bam (sorted by coordinate by default)
	###################################################################
	TOPHAT_OUT=$PROJECT_DIR/TopHat/$Barcode
	if [ $RUN_TOPHAT -eq 1 ]; then
		mkdir -p $TOPHAT_OUT
		Merged_FastQ=$(printf ",%s" "${Trimmed_FastQ_array[@]}")
		Merged_FastQ=${Merged_FastQ:1} # remove first comma 
		echo -e "\ntophat2 -p $NT --library-type $LIB_TYPE --output-dir $TOPHAT_OUT --max-multihits $TH_MH $TH_PREFILTER -transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX $BOWTIE2_INDEX_BASE $Merged_FastQ"

		time tophat2 \
			--num-threads $NT \
			--library-type $LIB_TYPE \
			--output-dir $TOPHAT_OUT \
			--max-multihits $TH_MH \
			$TH_PREFILTER \
			--transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX \
			$BOWTIE2_INDEX_BASE $Merged_FastQ
			#--no-coverage-search # for SLX-9615.A016, SLX-9790.L2V18DRBC02

		# check out result
		if [ ! -s $TOPHAT_OUT/accepted_hits.bam ];then 
			echo -e "\e[031m$TOPHAT_OUT/accepted_hits.bam not found. tophat failed? \e[0m\n"
			exit
		fi
	fi # end of RUN_TOPHAT
	##################
	## 3. BAM2FASTQ ##
	##################
	# unmapped.fastq.gz   : un-paired fastq
	if [ $RUN_BAM2FQ -eq 1 ]; then
		# convert bam to compressed fastq
		if [ -s $TOPHAT_OUT/unmapped.bam ];then 
			echo -e "\nbamutils tofastq $TOPHAT_OUT/unmapped.bam | gzip -f -9 > $TOPHAT_OUT/unmapped.fastq.gz"
			time bamutils tofastq $TOPHAT_OUT/unmapped.bam | gzip -f -9 > $TOPHAT_OUT/unmapped.fastq.gz
		fi
	fi
	######################
	# 3. CoverageBed
	# Input: bam from tophat
	# (NB, to use 'bedtools genomecov', BAM *must* be sorted by position, 
	# which is already done by tophat)
	# Output: $sample.genomecov.txt
	######################
	if [ $RUN_GENOMECOV -eq 1 ]; then
		mkdir -p $PROJECT_DIR/Coverage/$Barcode
		echo -e "bedtools genomecov $TOPHAT_OUT/accepted_hits.bam";
		time bedtools genomecov -split -ibam $TOPHAT_OUT/accepted_hits.bam > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.txt
		time bedtools genomecov -bg -split -ibam $TOPHAT_OUT/accepted_hits.bam > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.bedgraph
	fi

	######################
	# 4. HTSeq
	# For paired-end data, the alignment have to be sorted either by read name or by alignment position.
	# 'sort -T $SCRATCH_OUT -k 1,1' will do a sort by read name, or
	# samtools sort -n $TOPHAT_OUT/accepted_hits.bam TOPHAT_OUT/accepted_hits.sorted.bam
	# output:
	######################
	HTSEQ_OUT=$PROJECT_DIR/HTSeq/$Barcode
	if [ $RUN_HTSEQ -eq 1 ]; then
		mkdir -p $HTSEQ_OUT

		# count based on illumina iGenome GTF
		echo -e "\nsamtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count --stranded=$HTSEQ_STRAND --type=$HTSEQ_TYPE --quiet --idattr=$HTSEQ_IDATTR - $GTF > $HTSEQ_OUT/$SLX.$Barcode.igenome.HTSeq.count.txt"
		time samtools view $TOPHAT_OUT/accepted_hits.bam | htseq-count \
			--stranded=$HTSEQ_STRAND \
			--type=$HTSEQ_TYPE \
			--quiet \
			--idattr=$HTSEQ_IDATTR \
			- $GTF > $HTSEQ_OUT/$SLX.$Barcode.igenome.HTSeq.count.txt
	fi
else
	echo -e "\e[031m Paired-end mode\e[0m\n"
	##################################################
	# 0. Run the inital fastqc 
	##################################################
	if [ $RUN_FASTQC -eq 1 ]; then
		mkdir -p $PROJECT_DIR/FastQC/$Barcode
		FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`
		for FastQ_file in $FASTQ_FILES	# per fastq file (per lane, both forward and reverse)
		do
			echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
			time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
		done
	fi

	Trimmed_Fw_FastQ_array=() #initialise
	Trimmed_Rv_FastQ_array=() #initialise

	FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*r_1.fq.gz | grep -v lost`
	Cell_lane_array=`for i in $FASTQ_FILES; do echo $i | cut -d/ -f7 | cut -d. -f3,4 ; done`
	## For each lane 
	for CELL_LANE in ${Cell_lane_array[*]} 
	do
		Cell=`echo $CELL_LANE | cut -d. -f1`
		Lane=`echo $CELL_LANE | cut -d. -f2`
		if [ $IS_SMARTER1 -eq 1 ];then
			Fw_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_1.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz
			Rv_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_2.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_2.fq.gz
		else
			Fw_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_2.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz
			Rv_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_1.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_2.fq.gz
		fi
		if [ ! -s $Fw_FastQ ];then
			echo -e "\e[031m$Fw_FastQ not found\e[0m\n"
			exit
		fi
		if [ ! -s $Rv_FastQ ];then
			echo -e "\e[031m$Rv_FastQ not found\e[0m\n"
			exit
		fi
		############################################################
		# 1-a. Trimming via trim_galore 	
		############################################################
		# After trim_galore, I am expecting:
		if [ $IS_SMARTER1 -eq 1 ];then
			Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1_val_1.fq.gz 
			Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2_val_2.fq.gz 
		else
			Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2_val_2.fq.gz 
			Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1_val_1.fq.gz 
		fi
		if [ $RUN_TRIM -eq 1 ]; then
			mkdir -p $PROJECT_DIR/Trim/$Barcode
			echo -e "\ntrim_galore $Fw_FastQ $Rv_FastQ"
			time trim_galore \
				--paired \
				--output_dir $PROJECT_DIR/Trim/$Barcode \
				--stringency $TR_STRNCY \
				--adapter $TR_ADAPTOR1 \
				--adapter2 $TR_ADAPTOR2 \
				$Fw_FastQ $Rv_FastQ
		fi # end of run trim
		####################################################################
		# 1-b. Triming Plasma samples prepared by SMARTer kit of clonetech #
		####################################################################
		if [ $RUN_TRIM_SMARTER -eq 1 ]; then
			mkdir -p $PROJECT_DIR/Trim/$Barcode
			if [ ! -s $Trimmed_Fw_FastQ ];then
				echo -e "\ncutadapt $Fw_FastQ $Rv_FastQ"
				time cutadapt \
					-a $TR_ADAPTOR1 \
					-A $TR_ADAPTOR2 \
					-q $TR_QUAL \
					-u $TR_CUT_LEN \
					-O $TR_STRNCY \
					-m $TR_LEN \
					-o $Trimmed_Fw_FastQ \
					-p $Trimmed_Rv_FastQ \
					--trim-n \
					--max-n 5 \
					$Fw_FastQ $Rv_FastQ &> ${Trimmed_Fw_FastQ%.fq.gz}.trimming_report.txt
			fi
		fi # end of run trim

		# scrap trimmed reads
		if [ -s $Trimmed_Fw_FastQ ];then
			Trimmed_Fw_FastQ_array+=($Trimmed_Fw_FastQ) #push it
		else
			echo -e "\e[031m$Trimmed_Fw_FastQ not found\e[0m\n"
			exit
		fi

		if [ -s $Trimmed_Rv_FastQ ];then
			Trimmed_Rv_FastQ_array+=($Trimmed_Rv_FastQ) #push it
		else
			echo -e "\e[031m$Trimmed_Rv_FastQ not found\e[0m\n"
			exit
		fi
	done # endof CELL_LANE
	###################################################################
	# 2. TOPHAT 
	# Map the reads for each sample to the reference genome           #
	# real    61m11.381s (for R3 with 2 threads) 82m22.024s (for R4)
	# output: $TOPHAT_OUT/accepted_hits.bam (sorted by coordinate by default)
	###################################################################
	TOPHAT_OUT=$PROJECT_DIR/TopHat/$Barcode
	if [ $RUN_TOPHAT -eq 1 ]; then
		mkdir -p $TOPHAT_OUT
		Merged_Fw_FastQ=$(printf ",%s" "${Trimmed_Fw_FastQ_array[@]}")
		Merged_Fw_FastQ=${Merged_Fw_FastQ:1} # remove first comma 

		Merged_Rv_FastQ=$(printf ",%s" "${Trimmed_Rv_FastQ_array[@]}")
		Merged_Rv_FastQ=${Merged_Rv_FastQ:1} # remove first comma 

		echo -e "\ntophat2 -p $NT --library-type $LIB_TYPE --output-dir $TOPHAT_OUT --max-multihits $TH_MH $TH_PREFILTER -transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX $BOWTIE2_INDEX_BASE $Merged_Fw_FastQ $Merged_Rv_FastQ"
		time tophat2 \
			--num-threads $NT \
			--library-type $LIB_TYPE \
			--output-dir $TOPHAT_OUT \
			--max-multihits $TH_MH \
			$TH_PREFILTER \
			--transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX \
			$BOWTIE2_INDEX_BASE \
			$Merged_Fw_FastQ $Merged_Rv_FastQ
		# check out result
		if [ ! -s $TOPHAT_OUT/accepted_hits.bam ];then 
			echo -e "\e[031m$TOPHAT_OUT/accepted_hits.bam not found. tophat failed? \e[0m\n"
			exit
		fi
	fi # end of RUN_TOPHAT
	###########################
	## 3. Unmapped BAM2FASTQ ##
	###########################
	if [ $RUN_UNMAP2FA -eq 1 ]; then
		# sort unmapped.bam by readname
		if [ ! -s $TOPHAT_OUT/unmapped.sorted.bam ];then  
			# 1. assuming readnames of unmapped.bam are *NOT* uniq 
			#time samtools view $TOPHAT_OUT/unmapped.bam \
			#	| awk '{print $0}' | sort -S 60G -k 1,1 -k 2,2n -S 60G -T $PROJECT_DIR/scratch/$Barcode -u \
			#	| samtools view -bT $GENOME - > $TOPHAT_OUT/unmapped.sorted.bam 
			# 2. assuming readnames of unmapped.bam are uniq 
			echo -e "samtools sort -n -@ $NT -m 60G -o $TOPHAT_OUT/unmapped.sorted.bam -O bam -T $SLX.$Barcode $TOPHAT_OUT/unmapped.bam"
			time samtools sort -n -@ $NT -m 60G -o $TOPHAT_OUT/unmapped.sorted.bam -O bam -T $SLX.$Barcode $TOPHAT_OUT/unmapped.bam
		fi
		Unmapped_Fw_FastQ=$TOPHAT_OUT/unmapped.sorted_1.fastq.gz # paired read1
		Unmapped_Rv_FastQ=$TOPHAT_OUT/unmapped.sorted_2.fastq.gz # paired read2
		Unmapped_Un_FastQ=$TOPHAT_OUT/unmapped.sorted.fastq.gz   # un-paired 
		if [ ! -s $Unmapped_Fw_FastQ ] && [ ! -s $Unmapped_Rv_FastQ ]; then
			echo -e "bam bam2FastQ --readname --gzip --in $TOPHAT_OUT/unmapped.sorted.bam &> /dev/null"
			time bam bam2FastQ --readname --gzip --in $TOPHAT_OUT/unmapped.sorted.bam &> /dev/null
		fi
	fi
	######################
	# 5. featureCount
	# featureCounts automatically sorts reads by name 
	# if paired reads are not in consecutive positions in the SAM or BAM file
	# -T: No. of core
	# -a: GTF file
	# -Q: min MQ
	# --primary: Count primary alignments only
	# -s: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
	# -p: If specified, fragments (or templates) will be counted
	# -B: Count read pairs that have both ends successfully aligned only.
	# -C: Do not count read pairs that have their two ends mapping to different chromosomes
	######################	
	if [ $RUN_FTCNT -eq 1 ]; then
		FTCNT_OUT=$PROJECT_DIR/featureCount/$Barcode
		mkdir -p $FTCNT_OUT
		FTCNT_FILE=$FTCNT_OUT/$SLX.$Barcode.featureCount.$TR_PREFIX.$ENS_VER.txt
		echo -e "featureCounts -T $NT -a $GTF -Q $FTCNT_MQ -s $FTCNT_STRAND -p -C -o $FTCNT_FILE $TOPHAT_OUT/accepted_hits.bam\n"
		time featureCounts \
			-T $NT \
			-a $GTF \
			-Q $FTCNT_MQ \
			-s $FTCNT_STRAND \
			-p -C \
			-o $FTCNT_FILE \
			$TOPHAT_OUT/accepted_hits.bam
		time awk 'BEGIN{OFS="\t"}NR>2{print $1,$7}' $FTCNT_FILE | sort > ${FTCNT_FILE%.txt}.clean.txt
	fi
fi

echo -e "\e[32m$Barcode done for mapping\e[0m" 
