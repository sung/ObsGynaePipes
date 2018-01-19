#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 11/Apr/2014
# Last modified: 8/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/lib/modules.sh # defines necessary software to load 

export SLX="MY_SLX" # e.g. SLX-8080 
export Barcode="MY_BARCODE" # e.g. A001
export Lane="MY_LANE" # e.g. s_1 

# '--maxins' from bismark (Default: 500), Michelle identifies this from BioAnaliser
# see data/README_INFO to get the insertion size
BM_MAXINS=`grep $Barcode $LIB_SIZE_FILE | grep ${SLX/SLX-/} | cut -f7` 

if [[ ! $BM_MAXINS =~ ^-?[0-9]+$ ]]; then
	echo "No max lib insertion size found"
	exit
fi

mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"
echo -e  "Max insertion size=$BM_MAXINS"


FASTQ_FILES=`ls $FASTQ_DIR/$SLX.$Barcode*.fq.gz | grep -v lost`

# assuming paired-end
# [todo] what if single end?
if [ $IS_PE -eq 1 ]; then
	NUM_FW=`ls $FASTQ_DIR/$SLX.$Barcode*r_1.fq.gz | grep -v lost | wc -l`
	NUM_RV=`ls $FASTQ_DIR/$SLX.$Barcode*r_2.fq.gz | grep -v lost | wc -l`
	# check NO. of foward fq files are same with reverse
	if [ $NUM_FW -ne $NUM_RV ]; then
		echo -e "\e[031mNO of fw fq.gz ($NUM_FW) different from rv fq.gz ($NUM_RV)\e[0m\n"
		exit
	fi

	##################################################
	# 1. Get Flowcell & Lane info
	##################################################
	for FastQ_file in $FASTQ_FILES	# per fastq file
	do
		echo -e "\t\e[32mFASTQ=$FastQ_file\e[0m" 
		if [ -s $FastQ_file ]; then
			# forward-strand
			if [[ $FastQ_file =~ "r_1.fq.gz" ]]; then
				#[todo] assume the delimiter '.' works OK and lane info comes 4th
				dummy_cell_lane=`echo $FastQ_file | cut -d. -f3,4` # E.g C2N01ACXX.s_1 H8CP2ADXX.s_1 
				Cell_lane_array+=($dummy_cell_lane) # push it
			fi
		else
			echo -e "\e[031m$FastQ_file not found\e[0m\n"
			exit
		fi
	done # end of per FastQ_file 

	##################################################
	# 2. run initial fastqc 
	# output: results/SGA_v1/FastQC/A001/
	##################################################
	if [ $RUN_FASTQC -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/FastQC/$Barcode
		for FastQ_file in $FASTQ_FILES
		do
			if [ -s $FastQ_file ]; then
				echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
				time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
			else
				echo -e "\e[031m$FastQ_file not found\e[0m\n"
				exit
			fi
		done # end of per FastQ_file 
	fi

	############################################################
	# 3. Trim FastQ file per lane 	
	# input: .fq.gz
	# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1_val_1.fq.gz
	# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1_val_2.fq.gz
	# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz_trimming_report.txt
	############################################################
	for My_cell_lane in ${Cell_lane_array[*]} # C2N01ACXX.s_1, C2N01ACXX.s_2, H8CP2ADXX.s_1, H8CP2ADXX.s_2
	do
		Fw_FastQ=$FASTQ_DIR/$SLX.$Barcode.$My_cell_lane.r_1.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz
		Rv_FastQ=$FASTQ_DIR/$SLX.$Barcode.$My_cell_lane.r_2.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_2.fq.gz

		if [ $RUN_TRIM -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/Trim/$Barcode
			if [ -s $Fw_FastQ ] && [ -s $Rv_FastQ ]
			then
				echo -e "\ntrim_galore $Fw_FastQ $Rv_FastQ"
				time trim_galore \
					--paired --trim1 \
					--fastqc \
					--output_dir $PROJECT_DIR/Trim/$Barcode \
					--stringency $TR_STRNCY \
					$Fw_FastQ $Rv_FastQ
			else
				echo -e "\e[031m$Fw_FastQ or $Rv_FastQ not found\e[0m\n"
				exit
			fi
		fi # end of run trim

		# After trim_galore, I am expecting:
		Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$My_cell_lane.r_1_val_1.fq.gz 
		Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$My_cell_lane.r_2_val_2.fq.gz 
		if [ -s $Trimmed_Fw_FastQ ] && [ -s $Trimmed_Rv_FastQ ]
		then
			Trimmed_Fw_array+=($Trimmed_Fw_FastQ) #push it
			Trimmed_Rv_array+=($Trimmed_Rv_FastQ) #push it
		else
			echo -e "\e[031m$Trimmed_Fw_FastQ or $Trimmed_Rv_FastQ not found\e[0m\n"
			exit
		fi
	done # end of My_cell_lane (e.g. C48CWACXX.s_1)

	if [ ${#Trimmed_Fw_array[@]} -ne ${#Trimmed_Rv_array[@]} ];then
		echo -e "\e[031mArray size not same (${#Trimmed_Fw_array[@]} vs. ${#Trimmed_Rv_array[@]}) \e[0m\n"
		exit
	fi

	Merged_fw_FastQ=$(printf ",%s" "${Trimmed_Fw_array[@]}")
	Merged_rv_FastQ=$(printf ",%s" "${Trimmed_Rv_array[@]}")
	Merged_fw_FastQ=${Merged_fw_FastQ:1} # remove the first comma
	Merged_rv_FastQ=${Merged_rv_FastQ:1}

	###################
	## 4. BsExpress
	###################
	if [ $RUN_BS_EXPRESS -eq 1 ]; then
		source $BS_EXP_WRAPPER
	fi # end of BsExpress

	###########################
	# 5. Mapping with bismark #
	# input: .fq.gz (.fastq)
	# output: results/SGA_v1/Bismark/Alignment/A001/A001.r_1_val_1.fq.gz_bismark_bt2_PE_report.txt
	# output: results/SGA_v1/Bismark/Alignment/A001/A001.r_1_val_1.fq.gz_bismark_bt2_pe.bam
	# [todo] other parameters? (e.g. -L for seed length)
	###########################
	#[fixme]
	if [ $RUN_BM_ALN -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Bismark/Alignment/$Barcode
		# mapping
		echo -e "\nbismark -1 $Merged_fw_FastQ -2 $Merged_rv_FastQ"
		time bismark \
		--bam --bowtie2 --maxins $BM_MAXINS -p $NT \
		--temp_dir /disk2/ssg29/scratch \
		--output_dir /disk2/ssg29/Bismark/Alignment/$Barcode \
		$GENOME_DIR \
		-1 $Merged_fw_FastQ \
		-2 $Merged_rv_FastQ
		#--temp_dir $PROJECT_DIR/scratch/$Barcode \
		#--output_dir $PROJECT_DIR/Bismark/Alignment/$Barcode \
		#--temp_dir /disk2/ssg29/scratch \
		#--output_dir /disk2/ssg29/Bismark/Alignment/$Barcode \


	fi #end of bismark alignment


	####################################################
	# 7. Remove reads of poor mapping quality (MAPQ) ##
	# QC: by mapping quality (MAPQ) - remove poorly mapped reads based on MAPQ
	# do not sort bam at this stage as this breaks methyl_extractor later
	# input: BAM from bismark alignment
	# output: A011.r_1_val_1.fq.gz_bismark_bt2_pe.q20.bam
	####################################################

	#fixme
	#BAM_FILES=`ls $PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode*.bam`
	BAM_FILES=`ls /disk2/ssg29/Bismark/Alignment/$Barcode/$SLX.$Barcode*.bam`
	if [ $RUN_FLT_BY_MAPQ -eq 1 ]; then
		for My_bam in $BAM_FILES # per bam file (which is per-lane)
		do
			if [ -s $My_bam ]; then
				# filter
				echo -e "samtools view -bh -q $MAPQ $My_bam -o ${My_bam%.bam}.q$MAPQ.bam "
				time samtools view -bh -q $MAPQ $My_bam -o ${My_bam%.bam}.q$MAPQ.bam 

				# make a flagstat report (take out chimeric reads to get exact number of reads in flagstat -F 0x100)
				echo -e "samtools view -u -F 0x100 $My_bam | samtools flagstat - > $My_bam.flagstat"
				time samtools view -u -F 0x100 $My_bam | samtools flagstat - > $My_bam.flagstat

			else
				echo -e "\e[031m$My_bam not found\e[0m\n"
				exit
			fi

			if [ -s ${My_bam%.bam}.q$MAPQ.bam ]; then
				# flagstat
				echo -e "samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.bam.flagstat"
				time samtools view -u -F 0x100 ${My_bam%.bam}.q$MAPQ.bam | samtools flagstat - > ${My_bam%.bam}.q$MAPQ.bam.flagstat

				My_bam_array+=(${My_bam%.bam}.q$MAPQ.bam) #push it

			else
				echo -e "\e[031m ${My_bam%.bam}.q$MAPQ.bam not found\e[0m\n"
				exit
			fi
		done
	fi

	Merged_perlane_bam=$(printf " %s" "${My_bam_array[@]}")
	###########################
	## Merge per-lane bam file
	###########################
	#fixme
	#My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam 
	My_bam=/disk2/ssg29/Bismark/Alignment/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam 
	echo -e "samtools merge ${My_bam%.bam}.q$MAPQ.bam  $Merged_perlane_bam"
	time samtools merge ${My_bam%.bam}.q$MAPQ.bam  $Merged_perlane_bam

	##################################################
	# 6. run preseq 
	##################################################
	if [ $RUN_PRESEQ -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Preseq/$Barcode
		# sort inital bam file
		if [ -s $My_bam ]; then
			echo -e "samtools sort $My_bam ${My_bam%.bam}.sorted"
			time samtools sort $My_bam ${My_bam%.bam}.sorted
		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
		if [ -s ${My_bam%.bam}.sorted.bam ]; then
			echo -e "preseq lc_extrap -v -B -P -o $PROJECT_DIR/Preseq/$Barcode/$Barcode.preseq.txt ${My_bam%.bam}.sorted.bam"
			time preseq lc_extrap -v -B -P -o $PROJECT_DIR/Preseq/$Barcode/$Barcode.preseq.txt ${My_bam%.bam}.sorted.bam 
		else
			echo -e "\e[031m${My_bam%.bam}.sorted.bam not found\e[0m\n"
			exit
		fi
	fi

	##################################################### 
	# 8. Remove duplicate reads and poorly mapped reads #
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


	##################################
	# 9. MERGE Previous BAM FILE #
	# IF YOU HAVE MORE THAN TWO
	##################################
	This_bam=${My_bam%.bam}.q$MAPQ.deduplicated.bam
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
		My_bam=${My_bam%.bam}.q$MAPQ.deduplicated.bam # e.g. A011.r_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.q20.bam
	fi

	#########################################################################
	# 10. Bismark Methyl Calls 
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
		#fixme
		#mkdir_unless $PROJECT_DIR/Bismark/MethylCall/$Barcode
		mkdir_unless /disk2/ssg29/Bismark/MethylCall/$Barcode
		echo bismark_methylation_extractor
		if [ -s $My_bam ]; then
			time bismark_methylation_extractor \
				--paired-end --no_overlap --comprehensive \
				--merge_non_CpG --report --gzip \
				--ignore_r2 $BM_IGNORE2 \
				--buffer_size $BM_BUFFER_SIZE \
				--bedGraph \
				--cytosine_report \
				--output /disk2/ssg29/Bismark/MethylCall/$Barcode \
				--genome_folder $GENOME_DIR \
				$My_bam
				#[fixme]
				#--output $PROJECT_DIR/Bismark/MethylCall/$Barcode \
			
			CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.CpG_report.txt
			if [ -s $CpG_report ];then
				time cat $CpG_report | awk '{sum+=$4+$5}END{print "CpG Depth of Cov=",sum/NR}' > $PROJECT_DIR/Bismark/MethylCall/$Barcode/$Barcode.CpG.coverage.txt
				depth=(1 3 5 7 10)
				for d in ${depth[*]}
				do
					echo "calculating x$d coverage..."
					time cat $CpG_report | awk "{if(\$4+\$5>=$d){sum++}}END{print \"x$d CpG BoC=\",sum/NR*100}" >> $PROJECT_DIR/Bismark/MethylCall/$Barcode/$Barcode.CpG.coverage.txt

				done

				# make a file for MethylKit format
				time awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($3=="-"){str="R"}}{if($4+$5>0){cov=$4+$5;print $1"."$2,$1,$2,str,cov,$4/cov*100,$5/cov*100}}' $CpG_report > ${CpG_report%.txt}.methylkit.txt
			else
				echo -e "\e[031m$CpG_report not found\e[0m\n"
			fi
		else
			echo -e "\e[031m$My_bam not found\e[0m\n"
			exit
		fi
	fi

	##########################################################################
	# 11. sort the final bam file then index it
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
	# 12. Coverage Report based on the sorted bam file #
	# Input: sorted bam
	# Output: $Barcode.genomeccov.txt (genome coverage report by bedtool)
	# Output: $Barcode.coverage.txt (DoC and BoC)
	###################################################
	if [ $RUN_GENOMECOV -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Coverage/$Barcode
		if [ -s ${My_bam%.bam}.sorted.bam ]; then
			# get the genome size based on the BAM file
			genomeSize=`samtools view -H ${My_bam%.bam}.sorted.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'`
			Genome_Cov=$PROJECT_DIR/Coverage/$Barcode/$Barcode.bedtools.genomecov.txt

			echo -e "bedtools genomecov -ibam ${My_bam%.bam}.sorted.bam > $Genome_Cov"
			time bedtools genomecov -ibam ${My_bam%.bam}.sorted.bam > $Genome_Cov

			if [ -s $Genome_Cov ];then
				# Depth of Coveage
				time grep -P '^genome' $Genome_Cov | awk "{sum+=\$2*\$3} END { print \"depth of coverage= \",sum/$genomeSize}" > $PROJECT_DIR/Coverage/$Barcode/$Barcode.coverage.txt
				depth=(1 3 5 7 10)
				for d in ${depth[*]}
				do
					echo "calculating x$d coverage..."
				# Breadh of Coverage
					time grep -P '^genome' $Genome_Cov | awk "{if(\$2>=$d){sum+=\$3}} END {print \"x$d BoC=\",sum/$genomeSize*100}" >> $PROJECT_DIR/Coverage/$Barcode/$Barcode.coverage.txt

				done
			else
				echo -e "\e[031m$Genome_Cov not found\e[0m\n"
			fi
		else
			echo -e "\e[031m${My_bam%.bam}.sorted.bam not found\e[0m\n"
			exit
		fi
	fi

	##################################################################################
	## 13. MethylExtract
	##################################################################################
	if [ $RUN_METHYL_EXT -eq 1 ]; then
		if [ -s $ME_WRAPPER ]; then
			source $ME_WRAPPER
		else
			echo -e "cannot find $ME_WRAPPER"
			exit
		fi
	fi # end of $RUN_METHYL_EXT -eq 1 

	############
	## 14 BIS-SNP
	############
	if [ $RUN_BIS_SNP -eq 1 ]; then
		if [ -s $BIS_SNP_WRAPPER ]; then
			source $BIS_SNP_WRAPPER
		else
			echo -e "cannot find $BIS_SNP_WRAPPER"
			exit
		fi
	fi

fi # end of IS_PE

# [todo]
# remove merged read later if trimming is enabled?
echo -e "\e[32m$Barcode done\e[0m" 

