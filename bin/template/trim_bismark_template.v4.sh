#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 11/Apr/2014
# Last modified: 11/Jun/2014 

source /whale-data/ssg29/lib/sung.sh #defines 'mkdir_unless'

export SLX="MY_SLX" # e.g. SLX-8080 
export Barcode="MY_BARCODE" # e.g. A001

# '--maxins' from bismark (Default: 500), Michelle identifies this from BioAnaliser
# see data/README_INFO to get the insertion size
BM_MAXINS=`grep $Barcode $LIB_SIZE_FILE | grep ${SLX/SLX-/} | cut -f 7` 

mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
# per sample 
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"
echo -e  "Max insertion size=$BM_MAXINS"

<<fastq
-rw-r--r-- 1 1250 1250 400353641 Jun 10 21:13 /whale-data/mdj30/project/SGA/raw_reads/SLX-8080.A001.000000000-A8R45.s_1.r_1.fq.gz
-rw-r--r-- 1 1250 1250 429073517 Jun 10 21:14 /whale-data/mdj30/project/SGA/raw_reads/SLX-8080.A001.000000000-A8R45.s_1.r_2.fq.gz
fastq


FASTQ_FILES=`ls $FASTQ_DIR/$SLX.$Barcode*.fq.gz | grep -v lost`

# per fastq file
for FastQ_file in $FASTQ_FILES
do
	echo -e "\t\e[32mFASTQ=$FastQ_file\e[0m" 
	##################################################
	# 1. run initial fastqc 
	# output: results/SGA_v1/FastQC/A001/
	##################################################
	if [ $RUN_FASTQC -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/FastQC/$Barcode
		if [ -s $FastQ_File ]; then
			echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
			time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
		else
			echo -e "\e[031m$FastQ_file not found\e[0m\n"
			exit
		fi
	fi
		
	# assuming paired-end
	# [todo] what if single end?
	if [ $IS_PE -eq 1 ]; then
		if [ -s $FastQ_File ]; then
			# forward-strand
			if [[ $FastQ_file =~ "r_1.fq.gz" ]]; then
				FastQ_fw_array+=($FastQ_file) #push it
				#/whale-data/mdj30/project/SGA/raw_reads/SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz
				#/whale-data/mdj30/project/SGA/raw_reads/SLX-8074.A001.C2N01ACXX.s_2.r_1.fq.gz
			# reverse-strand
			else
				FastQ_rv_array+=($FastQ_file) #push it
			fi
		else
			echo -e "\e[031m$FastQ_file not found\e[0m\n"
			exit
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
		mkdir_unless $PROJECT_DIR/MergedRead/$Barcode

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
		mkdir_unless $PROJECT_DIR/Trim/$Barcode
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
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam # this will be generated after mapping
	#My_bam=/disk2/ssg29/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam # this will be generated after mapping
	if [ $RUN_BM_ALN -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Bismark/Alignment/$Barcode
		# mapping
		if [ -s $Merged_fw_FastQ ] && [ -s $Merged_rv_FastQ ]
		then
			echo -e "\nbismark -1 $Merged_fw_FastQ -2 $Merged_rv_FastQ"
			time bismark \
			--bam --bowtie2 --maxins $BM_MAXINS -p $NT \
			--temp_dir $PROJECT_DIR/scratch/$Barcode \
			--output_dir $PROJECT_DIR/Bismark/Alignment/$Barcode \
			$GENOME_DIR \
			-1 $Merged_fw_FastQ \
			-2 $Merged_rv_FastQ
			#--temp_dir $PROJECT_DIR/scratch/$Barcode \
			#--output_dir $PROJECT_DIR/Bismark/Alignment/$Barcode \
			#--temp_dir /disk2/ssg29/scratch \
			#--output_dir /disk2/ssg29 \
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

	####################################################
	# 7. Remove reads of poor mapping quality (MAPQ) ##
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
		mkdir_unless $PROJECT_DIR/Bismark/MethylCall/$Barcode
		echo bismark_methylation_extractor
		if [ -s $My_bam ]; then
			time bismark_methylation_extractor \
				--paired-end --no_overlap --comprehensive \
				--merge_non_CpG --report --gzip \
				--ignore_r2 $BM_IGNORE2 \
				--buffer_size $BM_BUFFER_SIZE \
				--bedGraph \
				--cytosine_report \
				--output /disk2/ssg29 \
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

