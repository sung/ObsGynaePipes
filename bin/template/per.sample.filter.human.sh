#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 17/Dec/2015
# Last modified: 17/Dec/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001

mkdir_unless $PROJECT_DIR	

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$PROJECT.$Barcode\e[0m\n"

TOPHAT_OUT=$PROJECT_DIR/TopHat/$Barcode
TOPHAT_OUT2=$PROJECT_DIR/TopHat/$Barcode/unmapped         # unmapped top dir
Unmapped_BAM=$TOPHAT_OUT2/unmapped.bam                    # unmapped BAM (after 2-pass)
# assuming single-end
# [todo] what if paired end?
if [ $IS_SE -eq 1 ]; then
	Unmapped_FastQ=$TOPHAT_OUT2/unmapped.fq.gz
	# covert un-mapped of unmapped bam file to fastq
	if [ ! -s $Unmapped_FastQ ];then 
		if [ -s $Unmapped_BAM ]; then
			echo -e "\nbamutils tofastq $Unmapped_BAM | gzip -9 > $Unmapped_FastQ"
			time bamutils tofastq $Unmapped_BAM | gzip -9 > $Unmapped_FastQ
		else
			echo -e "\e[031m$Unmapped_BAM not found\e[0m\n"
			exit
		fi
		# or using 'awk'
		#time samtools view $Unmapped_BAM | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip -9 > $Unmapped_FastQ 
	fi

	if [ ! -s $Unmapped_FastQ ]; then
		echo -e "\e[031m$Unmapped_Fast not found\e[0m\n"
		exit
	fi

	##################################
	## 1. Preprocess unmapped.fq.gz ##
	## firstly by cutadapt (qc-trim)
	## followed by prinseq (lc-trim)
	## this setting based on SE125
	##################################
	Prep_Unmapped_FastQ=$TOPHAT_OUT2/unmapped.prep.fq.gz
	if [ $RUN_PREPROCESS_UNMAPPED -eq 1 ];then 
		if [ ! -s $Prep_Unmapped_FastQ ];then 
			echo -e "\npreprocessing $Unmapped_FastQ"
			time cutadapt \
				-q 30 \
				-m 50 \
				--max-n 5 \
				--trim-n \
				$Unmapped_FastQ \
				2> ${Unmapped_FastQ%.fq.gz}.cutadapt.log \
				| perl ~/Install/prinseq-lite/prinseq-lite.pl \
				-fastq stdin \
				-min_len 50 \
				-ns_max_n 5 \
				-lc_method dust \
				-lc_threshold 7 \
				-noniupac \
				-trim_tail_left 5 \
				-trim_tail_right 5 \
				-out_good stdout \
				-out_bad null \
				2> ${Unmapped_FastQ%.fq.gz}.cutadapt.prinseq.log \
				| gzip -9 > $Prep_Unmapped_FastQ 
		fi
	fi
	##############################################################
	## 2. Map the preprocessed unmapped reads against human ref ##
	##############################################################
	# via bowtie2 (end-to-end)
	if [ $RUN_BOWTIE_UNMAPPED -eq 1 ]; then
		if [ -s $Prep_Unmapped_FastQ ];then 
			if [ ! -s $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz ];then 
				echo -e "bowtie2 -U $Prep_Unmapped_FastQ | samtools view -f4 | awk | gzip > $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz"
				time bowtie2 \
					--threads $NT \
					-U $Prep_Unmapped_FastQ \
					-x $BOWTIE2_INDEX_BASE2 \
					2> $TOPHAT_OUT2/unmapped.prep.bt2.log \
					| samtools view -f4 - | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip -9 > $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz 
			fi
		else
			echo -e "\e[031m$Prep_Unmapped_FastQ not found.\e[0m\n"
			exit
		fi
	fi
	# via bowtie2 --local
	if [ $RUN_BOWTIE_LOCAL_UNMAPPED -eq 1 ]; then
		if [ -s $Prep_Unmapped_FastQ ];then 
			if [ ! -s $TOPHAT_OUT2/unmapped.prep.bt2.local.fq.gz ];then 
				echo -e "bowtie2 --local -U $Prep_Unmapped_FastQ | samtools view -f4 | awk | gzip > $TOPHAT_OUT2/unmapped.prep.bt2.local.fq.gz"
				time bowtie2 \
					--local \
					--threads $NT \
					-U $Prep_Unmapped_FastQ \
					-x $BOWTIE2_INDEX_BASE2 \
					2> $TOPHAT_OUT2/unmapped.prep.bt2.local.log \
					| samtools view -f4 - | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip -9 > $TOPHAT_OUT2/unmapped.prep.bt2.local.fq.gz 
			fi
		else
			echo -e "\e[031m$Prep_Unmapped_FastQ not found.\e[0m\n"
			exit
		fi
	fi
	# via bbmap
	#http://seqanswers.com/forums/showthread.php?t=42552
	if [ $RUN_BBMAP_UNMAPPED -eq 1 ]; then
		if [ -s $Prep_Unmapped_FastQ ];then 
			time ~/Install/bbmap/bbmap.sh \
				minid=0.95 \
				maxindel=3 \
				bwr=0.16 \
				bw=12 \
				quickmatch fast \
				minhits=2 \
				path=~/Pipelines/ \
				qtrim=rl \
				trimq=10 \
				untrim \
				-Xmx23g \
				in=$Prep_Unmapped_FastQ \
				outu=$TOPHAT_OUT2/unmapped.prep.bbmap.fq.gz \
				outm=$TOPHAT_OUT2/mapped.prep.bbmap.fq.gz
		else
			echo -e "\e[031m$Prep_Unmapped_FastQ not found.\e[0m\n"
			exit
		fi
	fi
else
	echo -e "\e[031m Paired-end mode\e[0m\n"
	Unmapped_Fw_FastQ=$TOPHAT_OUT2/unmapped.sorted_1.fastq.gz # paired unmapped read1 (after 2-pass)
	Unmapped_Rv_FastQ=$TOPHAT_OUT2/unmapped.sorted_2.fastq.gz # paired unmapped read2 (after 2-pass)

	# covert un-mapped of unmapped bam file to fastq
	if [ ! -s $Unmapped_Fw_FastQ ] && [ ! -s $Unmapped_Rv_FastQ ]; then
		# sort unmapped.bam by readname
		if [ ! -s ${Unmapped_BAM%.bam}.sorted.bam ]; then
			echo -e "samtools sort -n -@ $NT -m 60G -o ${Unmapped_BAM%.bam}.sorted.bam -O bam -T $SLX.$Barcode $Unmapped_BAM"
			time samtools sort -n -@ $NT -m 60G -o ${Unmapped_BAM%.bam}.sorted.bam -O bam -T $SLX.$Barcode $Unmapped_BAM
		fi
		echo -e "bam bam2FastQ --readname --gzip --in $TOPHAT_OUT2/unmapped.sorted.bam &> /dev/null"
		time bam bam2FastQ --readname --gzip --in $TOPHAT_OUT2/unmapped.sorted.bam &> /dev/null
	fi
	##################################
	## 1. Preprocess unmapped.fq.gz ##
	## firstly by cutadapt (qc-trim)
	## followed by prinseq (lc-trim)
	## this setting based on PE75
	##################################
	Prep_Unmapped_Fw_FastQ=$TOPHAT_OUT2/unmapped.prep_1.fastq.gz # preprocessed unmapped read1 (after cutadapt & prinseq)
	Prep_Unmapped_Rv_FastQ=$TOPHAT_OUT2/unmapped.prep_2.fastq.gz # preprocessed unmapped read2 (after cutadapt & prinseq)
	Prep_Unmapped_Fw_Sg_FastQ=$TOPHAT_OUT2/unmapped.prep_1_singletons.fastq.gz # preprocessed unmapped singleton read2 (after cutadapt & prinseq)
	Prep_Unmapped_Rv_Sg_FastQ=$TOPHAT_OUT2/unmapped.prep_2_singletons.fastq.gz # preprocessed unmapped read2 (after cutadapt & prinseq)
	Prep_Unmapped_FastQ=$TOPHAT_OUT2/unmapped.prep.fastq.gz # preprocessed unmapped read1 (after cutadapt & prinseq)
	if [ $RUN_PREPROCESS_UNMAPPED -eq 1 ];then 
		##################################
		## Input: Paired unmapped reads ##
		##################################
		if [ ! -s $Prep_Unmapped_Fw_FastQ ] && [ ! -s $Prep_Unmapped_Rv_FastQ ];then 
			echo -e "\ncutadapt $Unmapped_Fw_FastQ $Unmapped_Rv_FastQ"
			time cutadapt \
				-m 20 \
				-q 30 \
				--max-n 0.3 \
				--trim-n \
				-o $TOPHAT_OUT2/unmapped.cutadapt_1.fastq.gz  \
				-p $TOPHAT_OUT2/unmapped.cutadapt_2.fastq.gz \
				$Unmapped_Fw_FastQ $Unmapped_Rv_FastQ &> $TOPHAT_OUT2/unmapped.cutadapt.log

			# decompress cutadapt-processed fq.gz
			time zcat $TOPHAT_OUT2/unmapped.cutadapt_1.fastq.gz > $TOPHAT_OUT2/unmapped.cutadapt_1.fastq
			time zcat $TOPHAT_OUT2/unmapped.cutadapt_2.fastq.gz > $TOPHAT_OUT2/unmapped.cutadapt_2.fastq
			# run prinseq-lite.pl
			echo -e "perl ~/Install/prinseq-lite/prinseq-lite.pl -fastq $TOPHAT_OUT2/unmapped.cutadapt_1.fastq -fastq2 $TOPHAT_OUT2/unmapped.cutadapt_2.fastq"
			time perl ~/Install/prinseq-lite/prinseq-lite.pl \
				-noniupac \
				-min_len 20 \
				-ns_max_p 30 \
				-trim_tail_left 5 \
				-trim_tail_right 5 \
				-lc_threshold 7 \
				-lc_method dust \
				-out_bad null \
				-out_good $TOPHAT_OUT2/unmapped.prep \
				-fastq $TOPHAT_OUT2/unmapped.cutadapt_1.fastq \
				-fastq2 $TOPHAT_OUT2/unmapped.cutadapt_2.fastq \
				2> $TOPHAT_OUT2/unmapped.prinseq.log

			# compress preprocessed paired fastq files
			time gzip -9 -f ${Prep_Unmapped_Fw_FastQ%.gz}
			time gzip -9 -f ${Prep_Unmapped_Rv_FastQ%.gz}
			# compress preprocessed un-paired (i.e. singletons) fastq files
			time gzip -9 -f ${Prep_Unmapped_Fw_Sg_FastQ%.gz}
			time gzip -9 -f ${Prep_Unmapped_Rv_Sg_FastQ%.gz}
			# remove uncompressed cutadapt fastq
			time rm $TOPHAT_OUT2/unmapped.cutadapt_1.fastq*
			time rm $TOPHAT_OUT2/unmapped.cutadapt_2.fastq*
		fi
		######################################################
		## INPUT: Singleton (i.e. un-paired) unmapped reads ##
		######################################################
		if [ ! -s $Prep_Unmapped_FastQ ];then 
			echo -e "\ncutadapt|prinseq-lite.pl $TOPHAT_OUT/unmapped.sorted.fastq.gz $TOPHAT_OUT2/unmapped.sorted.fastq.gz $Prep_Unmapped_Fw_Sg_FastQ $Prep_Unmapped_Rv_Sg_FastQ"
			time zcat $TOPHAT_OUT/unmapped.sorted.fastq.gz $TOPHAT_OUT2/unmapped.sorted.fastq.gz $Prep_Unmapped_Fw_Sg_FastQ $Prep_Unmapped_Rv_Sg_FastQ \
				| cutadapt \
				-m 20 \
				-q 30 \
				--max-n 0.3 \
				--trim-n \
				- \
				2>> $TOPHAT_OUT2/unmapped.cutadapt.log \
				| perl ~/Install/prinseq-lite/prinseq-lite.pl \
				-fastq stdin \
				-noniupac \
				-min_len 20 \
				-ns_max_p 30 \
				-trim_tail_left 5 \
				-trim_tail_right 5 \
				-lc_threshold 7 \
				-lc_method dust \
				-out_good stdout \
				-out_bad null \
				2>> $TOPHAT_OUT2/unmapped.prinseq.log \
				| gzip -9 -f > $Prep_Unmapped_FastQ 
		fi
	fi
	##############################################################
	## 2. Map the preprocessed unmapped reads against human ref ##
	##############################################################
	# via bowtie2 (end-to-end)
	if [ $RUN_BOWTIE_UNMAPPED -eq 1 ]; then
		################################################
		## Input: pre-processed unmapped paired reads ##
		## Output: unmapped.prep.bt2.singletons.fq.gz
		## Output: unmapped.prep.bt2.fq.1.gz
		## Output: unmapped.prep.bt2.fq.2.gz
		## Output: unmapped.prep.bt2.bam
		################################################
		if [ -s $Prep_Unmapped_Fw_FastQ ] && [ -s $Prep_Unmapped_Rv_FastQ ];then 
			echo -e "bowtie2 -1 $Prep_Unmapped_Fw_FastQ -2 $Prep_Unmapped_Rv_FastQ\n"
			time bowtie2 \
				--threads $NT \
				-1 $Prep_Unmapped_Fw_FastQ \
				-2 $Prep_Unmapped_Rv_FastQ \
				-x $BOWTIE2_INDEX_BASE2 \
				--un-gz $TOPHAT_OUT2/unmapped.prep.bt2.singletons.fq.gz \
				--un-conc-gz $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz \
				2> $TOPHAT_OUT2/unmapped.prep.bt2.log | samtools view -bS - > $TOPHAT_OUT2/unmapped.prep.bt2.bam
		else
			echo -e "\e[031m$Prep_Unmapped_Fw_FastQ and $Prep_Unmapped_Rv_FastQ not found.\e[0m\n"
			exit
		fi
		###############################################################
		## Input: pre-processed unmapped un-paired (singleton) reads ##
		## Output: unmapped.prep.bt2.fq.gz
		###############################################################
		echo -e "bowtie2 -U $Prep_Unmapped_FastQ\n"
		if [ -s $Prep_Unmapped_FastQ ];then 
			time bowtie2 \
				--threads $NT \
				-U $Prep_Unmapped_FastQ \
				-x $BOWTIE2_INDEX_BASE2 \
				2>> $TOPHAT_OUT2/unmapped.prep.bt2.log \
				| samtools view -f4 - | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip -9 -f > $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz 
		else
			echo -e "\e[031m$Prep_Unmapped_FastQ not found.\e[0m\n"
			exit
		fi

		echo -e "zcat $TOPHAT_OUT2/unmapped.prep.bt2.singletons.fq.gz $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz | gzip -9 -f > $TOPHAT_OUT2/unmapped.prep.bt2.singletons.merged.fq.gz"
		time zcat $TOPHAT_OUT2/unmapped.prep.bt2.singletons.fq.gz $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz | gzip -9 -f > $TOPHAT_OUT2/unmapped.prep.bt2.singletons.merged.fq.gz
		time rm $TOPHAT_OUT2/unmapped.prep.bt2.singletons.fq.gz $TOPHAT_OUT2/unmapped.prep.bt2.fq.gz
	fi
	########################################
	## FINAL HUMAN-FILTERED UNMAPPED READS #
	########################################
	# $TOPHAT_OUT2/unmapped.prep.bt2.singletons.merged.fq.gz # un-paired (singltons)
	# $TOPHAT_OUT2/unmapped.prep.bt2.fq.1.gz # paired 1
	# $TOPHAT_OUT2/unmapped.prep.bt2.fq.2.gz # paired 2
fi

echo -e "\e[32m$Barcode done for per.sample.filter.human\e[0m" 
