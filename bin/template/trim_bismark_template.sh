#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 11/Apr/2014
# Last modified: 8/May/2014 

source /whale-data/ssg29/lib/sung.sh #defines 'mkdir_unless'

Barcode="MY_BARCODE" # e.g. A001

mkdir_unless $PROJECT_DIR/FastQC/$Barcode
mkdir_unless $PROJECT_DIR/MergedRead/$Barcode
mkdir_unless $PROJECT_DIR/Trim/$Barcode
mkdir_unless $PROJECT_DIR/Bismark/Alignment/$Barcode
mkdir_unless $PROJECT_DIR/Bismark/MethylCall/$Barcode
mkdir_unless $PROJECT_DIR/Coverage/$Barcode
mkdir_unless $PROJECT_DIR/scratch/$Barcode

# per sample 
printf "\n\e[32mBarcode=$Barcode\e[0m\n"

FASTQ_FILES=`ls $FASTQ_DIR/SLX-*$Barcode*.fq.gz | grep -v lost`

# per fastq file
# assuming paired-end
# [todo] what if single end?
for FastQ_file in $FASTQ_FILES
do
	# 20 fq.gz for A001 (10 r_1 and 10 r_2)
	echo -e "\t\e[32mFASTQ=$FastQ_file\e[0m" 
	##################################################
	# 0. run initial fastqc 
	# output: results/SGA_v1/FastQC/A001/
	##################################################
	if [ $RUN_FASTQC -eq 1 ]; then
		if [ -s $FastQ_File ]; then
			echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
			time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
		else
			echo -e "\e[031m$FastQ_file not found\e[0m\n"
			exit
		fi
	fi
		
	if [ $IS_PE -eq 1 ]; then
		# forward-strand
		if [[ $FastQ_file =~ "r_1.fq.gz" ]]; then
			FastQ_fw_array+=($FastQ_file) #push it
			#/whale-data/mdj30/project/SGA/raw_reads/SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz
			#/whale-data/mdj30/project/SGA/raw_reads/SLX-8074.A001.C2N01ACXX.s_2.r_1.fq.gz
		# reverse-strand
		else
			FastQ_rv_array+=($FastQ_file) #push it
		fi
	fi # end of IS_PE
done # end of per FastQ_file 

if [ $IS_PE -eq 1 ]; then
	##################################################
	# 1. merge fq per strand 
	# time: 23m36.347s for A001.r_1.fq.gz (40G)
	# input: FastQ_rv_array and FastQ_fw_array
	# output: results/SGA_v1/MergedRead/A001/A001.r_1.fq.gz
	# output: results/SGA_v1/MergedRead/A001/A001.r_2.fq.gz
	##################################################
	dummy_fw=$(printf " %s" "${FastQ_fw_array[@]}")
	dummy_rv=$(printf " %s" "${FastQ_rv_array[@]}")
	Merged_fw_FastQ=$PROJECT_DIR/MergedRead/$Barcode/$Barcode.r_1.fq.gz
	Merged_rv_FastQ=$PROJECT_DIR/MergedRead/$Barcode/$Barcode.r_2.fq.gz
	if [ $RUN_MERGE_FASTQ -eq 1 ]; then
		## Merge forward read 
		echo -e "cat $dummy_fw\n"
		time cat $dummy_fw > $Merged_fw_FastQ
		## Merge reverse read
		echo -e "cat $dummy_rv\n"
		time cat $dummy_rv > $Merged_rv_FastQ
	fi

	###############################
	# 2. FastQC for merged fastq?
	# output: results/SGA_v1/FastQC/A001/
	###############################
	if [ $RUN_MERGED_FASTQC -eq 1 ]; then
		# forward merged fastq
		if [ -s $Merged_fw_FastQ ]; then
			time fastqc $Merged_fw_FastQ -o $PROJECT_DIR/FastQC/$Barcode 
		else
			echo -e "\e[031m$Merged_fw_FastQ not found\e[0m\n"
			exit
		fi
		# reverse merged fastq
		if [ -s $Merged_rv_FastQ ]; then
			time fastqc $Merged_rv_FastQ -o $PROJECT_DIR/FastQC/$Barcode 
		else
			echo -e "\e[031m$Merged_rv_FastQ not found\e[0m\n"
			exit
		fi
	fi # end of FastQC

	############################################################
	# 3. Trim merged FastQ file 	
	# time: 831m0.690s (13 hrs) for A001.r_1_val_1.fq.gz (24G) 
	# input: .fq.gz
	# output: results/SGA_v1/Trim/A001/A001.r_1_val_1.fq.gz (24G)
	# output: results/SGA_v1/Trim/A001/A001.r_2_val_2.fq.gz (23G)
	# output: results/SGA_v1/Trim/A001/A001.r_1.fq.gz_trimming_report.txt
	# output: results/SGA_v1/Trim/A001/A001.r_2.fq.gz_trimming_report.txt
	############################################################
	if [ $RUN_TRIM -eq 1 ]; then
		if [ -s $Merged_fw_FastQ ] && [ -s $Merged_rv_FastQ ]
		then
			echo -e "\ntrim_galore $Merged_fw_FastQ $Merged_rv_FastQ"
			time trim_galore \
				--paired --trim1 \
				--fastqc \
				--output_dir $PROJECT_DIR/Trim/$Barcode \
				--stringency $TR_STRNCY \
				$Merged_fw_FastQ $Merged_rv_FastQ
		else
			echo -e "\e[031m$Merged_fw_FastQ or $Merged_rv_FastQ not found\e[0m\n"
			exit
		fi
	fi # end of run trim

	Merged_fw_FastQ=${Merged_fw_FastQ/MergedRead/Trim} # 'MergedRead' => 'Trim'
	Merged_rv_FastQ=${Merged_rv_FastQ/MergedRead/Trim}

	Merged_fw_FastQ=${Merged_fw_FastQ%.fq.gz}_val_1.fq.gz # A001.r_1.fq.gz => A001.r_1_val_1.fq.gz 
	Merged_rv_FastQ=${Merged_rv_FastQ%.fq.gz}_val_2.fq.gz # A001.r_2.fq.gz => A001.r_2_val_2.fq.gz 
	###########################
	## filter by read length ##
	# input: results/SGA_v1/Trim/A001/A001.r_1_val_1.fq.gz (24G)
	# input: results/SGA_v1/Trim/A001/A001.r_2_val_2.fq.gz (23G)
	# output: results/SGA_v1/Trim/A001/A001.r_1_val_1.less70.fq.gz (24G)
	###########################
	if [ $RUN_FLT_BY_LEN -eq 1 ]; then
		if [ -s $Merged_fw_FastQ ] && [ -s $Merged_rv_FastQ ];then
			#echo -e "zcat $Merged_fw_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > ${Merged_fw_FastQ%.fq.gz}.less$FLT_FQ_LEN.fq.gz"
			#time zcat $Merged_fw_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > ${Merged_fw_FastQ%.fq.gz}.less$FLT_FQ_LEN.fq.gz
			# butterfly--only
			echo -e "zcat $Merged_fw_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > /disk2/ssg29/$Barcode.r_1_val_1.less$FLT_FQ_LEN.fq.gz"
			time zcat $Merged_fw_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > /disk2/ssg29/$Barcode.r_1_val_1.less$FLT_FQ_LEN.fq.gz 

			#echo -e "zcat $Merged_rv_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > ${Merged_rv_FastQ%.fq.gz}.less$FLT_FQ_LEN.fq.gz"
			#time zcat $Merged_rv_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > ${Merged_rv_FastQ%.fq.gz}.less$FLT_FQ_LEN.fq.gz
			echo -e "zcat $Merged_rv_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > /disk2/ssg29/$Barcode.r_2_val_2.less$FLT_FQ_LEN.fq.gz"
			time zcat $Merged_rv_FastQ | $FLT_FQ_BIN $FLT_FQ_LEN | gzip -vc > /disk2/ssg29/$Barcode.r_2_val_2.less$FLT_FQ_LEN.fq.gz

			#time rm $Merged_fw_FastQ $Merged_rv_FastQ

		else
			echo -e "\e[031m$Merged_fw_FastQ or $Merged_rv_FastQ not found\e[0m\n"
			exit
		fi
	fi # end of run filter by read len 

	Merged_fw_FastQ=${Merged_fw_FastQ%.fq.gz}.less$FLT_FQ_LEN.fq.gz # A001.r_1_val_1.fq.gz => A001.r_1_val_1.less70.fq.gz
	Merged_rv_FastQ=${Merged_rv_FastQ%.fq.gz}.less$FLT_FQ_LEN.fq.gz
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$Barcode.r_1_val_1.less$FLT_FQ_LEN.fq.gz_bismark_bt2_pe.bam # this will be generated after mapping
	###########################
	# 4. Mapping with bismark #
	# time: A001 started Friday and still running on Tuewday... 
	# input: .fq.gz (.fastq)
	# output: results/SGA_v1/Bismark/Alignment/A001/A001.r_1_val_1.fq.gz_bismark_bt2_PE_report.txt
	# output: results/SGA_v1/Bismark/Alignment/A001/A001.r_1_val_1.fq.gz_bismark_bt2_pe.bam
	# [todo] other parameters? (e.g. -L for seed length)
	###########################
	if [ $RUN_BM_ALN -eq 1 ]; then
		# mapping
		if [ -s $Merged_fw_FastQ ] && [ -s $Merged_rv_FastQ ]
		then
			echo -e "\nbismark -1 $Merged_fw_FastQ -2 $Merged_rv_FastQ"
			time bismark \
			--bam --bowtie2 --maxins $BM_MAXINS -p $NT \
			--output_dir $PROJECT_DIR/Bismark/Alignment/$Barcode \
			--temp_dir $PROJECT_DIR/scratch/$Barcode \
			$GENOME_DIR \
			-1 $Merged_fw_FastQ \
			-2 $Merged_rv_FastQ
			#--output_dir $PROJECT_DIR/Bismark/Alignment/$Barcode \
		else
			echo -e "\e[031m$Merged_fw_FastQ or $Merged_rv_FastQ not found\e[0m\n"
			exit
		fi

		# make a flagstat report (take out chimeric reads to get exact number of reads in flagstat -F 0x100)
		if [ -s $My_bam ]; then
			echo -e "samtools view -u -F 0x100 $My_bam | samtools flagstat - > $My_bam.flagstat"
			time samtools view -u -F 0x100 $My_bam | samtools flagstat - > $My_bam.flagstat
		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
	fi #end of bismark alignment

	####################################################
	# 5. Remove reads of poor mapping quality (MAPQ) ##
	# QC: by mapping quality (MAPQ) - remove poorly mapped reads based on MAPQ
	# do not sort bam at this stage as this breaks methyl_extractor later
	# input: BAM from bismark alignment
	# output: A011.r_1_val_1.fq.gz_bismark_bt2_pe.q20.bam
	####################################################
	if [ $RUN_FLT_BY_MAPQ -eq 1 ]; then
		if [ -s $My_bam ]; then
			# filter
			echo -e "samtools view -bh -q $MAPQ $My_bam -o ${My_bam%.bam}.q$MAPQ.bam "
			time samtools view -bh -q $MAPQ $My_bam -o ${My_bam%.bam}.q$MAPQ.bam 

		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi

		if [ -s ${My_bam%.bam}.q$MAPQ.bam ]; then
			# flagstat
			echo -e "samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.bam.flagstat"
			time samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.bam.flagstat

		else
			echo -e "\e[031m ${My_bam%.bam}.q$MAPQ.bam not found\e[0m\n"
			exit
		fi
	fi

	##################################################### 
	# 6. Remove duplicate reads and poorly mapped reads #
	# [todo] any alternative tool? (e.g. picard)
	# time: (NT=20 for A001)
	# input: bam from bismark
	# output: results/SGA_v1/Bismark/Alignment/A011/A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bam
	##################################################### 
	if [ $RUN_BM_DEDUP -eq 1 ]; then
		if [ -s ${My_bam%.bam}.q$MAPQ.bam ]; then
			## input: q$MAPQ.bam
			## output: q$MAPQ.deduplicated.bam
			echo deduplicate_bismark ${My_bam%.bam}.q$MAPQ.bam 
			time deduplicate_bismark \
				--paired --bam \
				${My_bam%.bam}.q$MAPQ.bam
		else
			echo -e "\e[031m ${My_bam%.bam}.q$MAPQ.bam not found\e[0m\n"
			exit
		fi

		# flagstat for q$MAPQ.deduplicated.bam
		# expecting  A011.r_1_val_1.fq.gz_bismark_bt2_pe.q20.deduplicated.bam
		if [ -s ${My_bam%.bam}.q$MAPQ.deduplicated.bam ]; then
			echo -e "samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.deduplicated.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.deduplicated.bam.flagstat"
			time samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.deduplicated.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.deduplicated.bam.flagstat
		else
			echo -e "\e[031m${My_bam%.bam}.q$MAPQ.deduplicated.bam not found\e[0m\n"
			exit
		fi
	fi

	######################
	# 6. MERGE BAM FILE #
	# IF YOU HAVE MORE THAN TWO
	# [todo] how to get $Other_bam??
	This_bam=${My_bam%.bam}.deduplicated.q$MAPQ.bam
	Other_bam=''
	Merged_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$Barcode.merged.bam
	if [ $RUN_MERGE_BAM -eq 1 ]; then
		if [ -s $This_bam ] && [ -s $Other_bam ]; then
			echo -e "samtools merge $Merged_bam $My_bam $Other_bam "
			time samtools merge $Merged_bam $My_bam $Other_bam 
			echo -e "samtools view -uq $MAPQ $Merged_bam "
			time samtools view -uq $MAPQ $Merged_bam 
		else
			echo -e "\e[031m$Other_bam or $My_bam not found\e[0m\n"
			exit
		fi
		My_bam=${Merged_bam%.bam}.q$MAPQ.bam # e.g. A011.merged.q20.bam
	else
		My_bam=${My_bam%.bam}.q$MAPQ.deduplicated.q$MAPQ.bam # e.g. A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.q20.bam
	fi

	#######################
	# 7. Methyl Extractor #
	# time:
	# input: A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.q20.bam 
	# output: results/SGA_v1/Bismark/MethylCall/A001/ 
	# 1) CpG_context_A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.txt.gz
	# 2) Non_CpG_context_A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.txt.gz
	# 3) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.M-bias_R2.png
	# 4) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bam_splitting_report.txt
	# 5) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.M-bias.txt
	# 6) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.M-bias_R1.png
	# 7) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bismark.cov
	# 8) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bedGraph
	# 9) A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.CpG_report.txt
	#
	# --comprehensive: merge different strand (OT, OB, CTOT, CTOB)
	# --merge_non_CpG: CpG vs non-CpG (rather than CpG, CHG, CGG)
	# --bedGraph calls 'bismark2bedGraph'
	# --cytosine_report calls 'bedGraph2cytosine', which is now renamed as 'coverage2cytosine' 
	# The option '--cytosine_report' above run this in one go, but if you prefer run separately:
	# alterantively, time coverage2cytosine --genome_folder $GENOME_DIR --dir $PROJECT_DIR/Bismark/MethylCall/$Barcode --output $Barcode.cov2cytosine.txt $Barcode.cov
	# this produces genome-wide methylation report for all cytosines in the genome 
	# --gzip: This option does not work on bedGraph and genome-wide cytosine reports as they are 'tiny' anyway
	# [todo] other parameters? (e.g. --ignore after M-bias plot)
	#######################
	if [ $RUN_BM_METHYL -eq 1 ]; then
		echo bismark_methylation_extractor
		if [ -s $My_bam ]; then
			time bismark_methylation_extractor \
				--paired-end --no_overlap --comprehensive \
				--merge_non_CpG --report --gzip \
				--ignore_r2 $BM_IGNORE2 \
				--buffer_size $BM_BUFFER_SIZE \
				--bedGraph \
				--cytosine_report \
				--output $PROJECT_DIR/Bismark/MethylCall/$Barcode \
				--genome_folder $GENOME_DIR \
				$My_bam
		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
	fi

	##########################################################################
	# 8. sort the final bam file then index it
	# Input: unsorted bam file 
	# Output: A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.q20.sorted.bam
	##########################################################################
	if [ $RUN_SAM_SORT -eq 1 ]; then
		echo -e "samtools sort $My_bam ${My_bam%.bam}.sorted"
		time samtools sort $My_bam ${My_bam%.bam}.sorted
		echo -e "samtools index ${My_bam%.bam}.sorted.bam "
		time samtools index ${My_bam%.bam}.sorted.bam 
	fi

	###################################################
	# 9. Coverage Report based on the sorted bam file #
	# Input: sorted bam
	# Output: $Barcode.genomeccov.txt (genome coverage report by bedtool)
	# Output: $Barcode.coverage.txt (DoC and BoC)
	###################################################
	if [ -s ${My_bam%.bam}.sorted.bam ]; then
		# get the genome size based on the BAM file
		genomeSize=`samtools view -H results/SGA_v1/Bismark/Alignment/A002/A002.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.q20.sorted.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'`

		time bedtools genomecov \
			-ibam ${My_bam%.bam}.sorted.bam > $PROJECT_DIR/Bismark/Coverage/$Barcode/$Barocode.genomecov.txt

		# Depth of Coveage
		time grep -P '^genome' $PROJECT_DIR/Bismark/Coverage/$Barcode/$Barocode.genomecov.txt | awk '{sum+=$2*$3} END { print "depth of coverage= ",sum/$genomeSize}' > $PROJECT_DIR/Bismark/Coverage/$Barcode/$Barcode.coverage.txt

		# Breadh of Coverage
		time grep -P '^genome\s+[123456789]' $PROJECT_DIR/Bismark/Coverage/$Barcode/$Barocode.genomecov.txt | awk '{sum+=$3} END {print "breadh of coveage= ", sum/$genomeSize}' >> $PROJECT_DIR/Bismark/Coverage/$Barcode/$Barcode.coverage.txt		
	fi

fi # end of IS_PE

# [todo]
# remove merged read later if trimming is enabled?
echo -e "\e[32m$Barcode done\e[0m" 

