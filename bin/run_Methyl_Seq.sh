#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 17/Apr/2014
# Last modified: 18/Jan/2015 
# A wrapper script to run trim_galore and bismark 
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

mkdir_unless $RESULT_DIR # $HOME/scratch/results	
mkdir_unless $RESULT_DIR/$SLX.$VERSION # $HOME/scratch/result/SLX-8080.v1 (aka. $PROJECT_DIR)	
mkdir_unless $TOP/logs
mkdir_unless $BIN_TOP/script
mkdir_unless $BIN_TOP/slurm

function make_run_script(){
    MY_TEMPLATE=$1 # e.g. $MAPPING_T
    MY_APP=$2 # e.g. SLX-8077.trim.bismark.A005.C48CWACXX.s_6.sh
    MY_BARCODE=$3 # e.g. A005 

    MY_CELL=$4 # e.g. C48CWACXX or chr1 (optional)
    MY_LANE=$5 # e.g. s_1 (optional)
    MY_CHUNK=$6 # e.g. 1 (optional)

    cp $MY_TEMPLATE $BIN_TOP/script/$MY_APP

    sed -i "s/MY_SLX/$SLX/" $BIN_TOP/script/$MY_APP
    sed -i "s/MY_VERSION/$VERSION/" $BIN_TOP/script/$MY_APP
    sed -i "s/MY_BARCODE/$MY_BARCODE/" $BIN_TOP/script/$MY_APP
	if [ -n "$MY_CELL" ];then #string is not null.
    	sed -i "s/MY_CELL/$MY_CELL/" $BIN_TOP/script/$MY_APP
	fi
	if [ -n "$MY_LANE" ];then #string is not null.
    	sed -i "s/MY_LANE/$MY_LANE/" $BIN_TOP/script/$MY_APP
	fi
	if [ -n "$MY_CHUNK" ];then #string is not null.
    	sed -i "s/MY_CHUNK/$MY_CHUNK/" $BIN_TOP/script/$MY_APP
	fi
}

<<fastq
SLX-8077.A005.C48CWACXX.s_1.r_1.fq.gz
SLX-8077.A005.C48CWACXX.s_1.r_2.fq.gz
SLX-8077.A005.C48CWACXX.s_2.r_1.fq.gz
SLX-8077.A005.C48CWACXX.s_2.r_2.fq.gz
fastq
# template/trim_bismark_template_hpc.v1.sh (fastqc, trim, mapping <- per lane)
# this function is deprecated
function run_per_lane_mapping(){
	USE_SLURM=$1 # 1 (use slurm) or 0 (do not use slurm)

	if [ ! -d $FASTQ_DIR ];then
		echo "$FASTQ_DIR not found"
		exit
	fi
	for MY_BARCODE in ${ALL_BARCODES[*]}
	do
		printf "\n\e[32mmapping $SLX.$MY_BARCODE\e[0m\n"
		#/home/ssg29/scratch/data/fastq/SLX-8075/SLX-8075.A008.C39MDACXX.s_1.r_1.fq.gz
		FASTQ_FILES=`ls $FASTQ_DIR/${SLX}.${MY_BARCODE}*r_1.fq.gz | grep -v lost`
		Cell_lane_array=`for i in $FASTQ_FILES; do echo $i | cut -d/ -f8 | cut -d. -f3,4 ; done`
		for MY_CELL_LANE in ${Cell_lane_array[*]} 
		do
			# fixme
			# skip for SLX-8077.A005 lane 1,2,3,4
			if [ $MY_CELL_LANE = "C48CWACXX.s_3" ] || [ $MY_CELL_LANE = "C48CWACXX.s_4" ]; then 
				echo skipping lane $MY_CELL_LANE
				continue
			fi

			echo -e "\t\e[32mCELL.LANE=$MY_CELL_LANE\e[0m" 

			MY_CELL=`echo $MY_CELL_LANE | cut -d. -f1`
			MY_LANE=`echo $MY_CELL_LANE | cut -d. -f2`
			MY_APP=$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.mapping.sh
			make_run_script $MAPPING_T $MY_APP $MY_BARCODE $MY_CELL $MY_LANE

			# submit this job to the scheduler
			if [ $USE_SLURM -eq 1 ]; then
				MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.mapping.sh
				make_slurm_script $MY_SLURM_SCRIPT $MY_APP "mapping" $MY_BARCODE $MY_CELL $MY_LANE 1
				cd $TOP/logs
				echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
				sbatch $SLURM_OPT $MY_SLURM_SCRIPT
			# run the script locally
			else
				echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.mapping.log"
				time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.mapping.log &
			fi # end of USE_SLURM
		done # end of MY_CELL_LANE

	done # end of MY_BARCODE
} # end of run_per_lane_mapping

<<chunk
SLX-8075.A011.000000000-A7DHM.s_1.r_1.1.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_2.1.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_1.2.fastq.gz
SLX-8075.A011.000000000-A7DHM.s_1.r_2.2.fastq.gz
chunk
function run_per_lane_per_chunk_mapping(){
	if [ ! -d $FASTQ_DIR ];then
		echo "$FASTQ_DIR not found"
		exit
	fi
	for MY_BARCODE in ${ALL_BARCODES[*]}
	do
		printf "\n\e[32mmapping $SLX.$MY_BARCODE\e[0m\n"
		#/home/ssg29/scratch/data/fastq/SLX-8075/SLX-8075.A008.C39MDACXX.s_1.r_1.fq.gz
		FASTQ_FILES=`ls $FASTQ_DIR/${SLX}.${MY_BARCODE}*r_1.fq.gz | grep -v lost`
		Cell_lane_array=`for i in $FASTQ_FILES; do echo $i | cut -d/ -f8 | cut -d. -f3,4 ; done`
		for MY_CELL_LANE in ${Cell_lane_array[*]}
		do
			echo -e "\t\e[32mCELL.LANE=$MY_CELL_LANE\e[0m" 

			MY_CELL=`echo $MY_CELL_LANE | cut -d. -f1`
			MY_LANE=`echo $MY_CELL_LANE | cut -d. -f2`

			#NUM_CHUNK within config/methyl_seq.config
			for (( MY_CHUNK=1; MY_CHUNK<=$NUM_CHUNK; MY_CHUNK++ ))
			do
				MY_APP=$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CHUNK.mapping.$sh

				make_run_script $MAPPING_CHUNK_T $MY_APP $MY_BARCODE $MY_CELL $MY_LANE $MY_CHUNK

				MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CHUNK.mapping.sh
				make_slurm_script $MY_SLURM_SCRIPT $MY_APP "mapping" $MY_BARCODE $MY_CELL $MY_LANE $MY_CHUNK

				# submit this job
				cd $TOP/logs
				echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
				sbatch $SLURM_OPT $MY_SLURM_SCRIPT
			done # end of MY_CHUNK
		done # end of MY_CELL_LANE

	done # end of MY_BARCODE
} # end of run_per_lane_per_chunk

function run_pipeline(){
	MY_CALLER=$1 # fastq.split, merge.perlane.bam, bismark.caller, methylextract, bissnp, coverage.report
	USE_SLURM=$2 # 1 (use slurm) or 0 (do not use slurm)
	PER_CHR=$3 # 1 (a flag to run by chromosome-wise)

	TEMPLATE=$BIN_TOP/template/$MY_CALLER.sh
	if [ ! -f $TEMPLATE ];then
		echo "$TEMPLATE not found"
		exit
	fi

	printf "\e[32mdoing $MY_CALLER for $SLX\e[0m\n"
	#############################
	# per SLX (count read.base) #
	#############################
	if [[ $MY_CALLER =~ ^per.run ]]; then
		MY_APP=$SLX.$VERSION.$MY_CALLER.sh
		make_run_script $TEMPLATE $MY_APP

		# run the script locally
		echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.log"
		nice -n 10 time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.log &
	############
	# per lane #
	############
	elif [[ $MY_CALLER =~ ^per.lane ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\n\e[32m$MY_CALLER for $SLX.$MY_BARCODE\e[0m\n"
			#/home/ssg29/data/fastq/SLX-8771/SLX-8771.FLD0181.000000000-AD9G6.s_1.r_1.fq.gz
			FASTQ_FILES=`ls $FASTQ_DIR/${SLX}.${MY_BARCODE}.*r_1.fq.gz | grep -v lost`
			Cell_lane_array=`for i in $FASTQ_FILES; do echo $i | cut -d/ -f7 | cut -d. -f3,4 ; done`
			for MY_CELL_LANE in ${Cell_lane_array[*]} 
			do
				MY_CELL=`echo $MY_CELL_LANE | cut -d. -f1`
				MY_LANE=`echo $MY_CELL_LANE | cut -d. -f2`
				MY_APP=$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.sh
				make_run_script $TEMPLATE $MY_APP $MY_BARCODE $MY_CELL $MY_LANE

				# submit this job to the scheduler
				if [ $USE_SLURM -eq 1 ]; then
					MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.sh
					make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $MY_CELL $MY_LANE 1
					cd $TOP/logs
					echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
					sbatch $SLURM_OPT $MY_SLURM_SCRIPT
				# run the script locally
				else
					echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.log"
					nice -n 10 time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.log &
				fi # end of USE_SLURM
			done # end of MY_CELL_LANE
		done
	###########################
	# per sample (or barcode) #
	###########################
	else
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			if [ -n "$PER_CHR" ];then #string is not null.
				for Chr in ${UCSC_CHR[*]}
				do 
					printf "\n\e[32m doing $MY_CALLER for $SLX.$MY_BARCODE.$Chr\e[0m\n"
					MY_APP=$SLX.$MY_BARCODE.$MY_CALLER.$Chr.sh
					MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CALLER.$Chr.sh
					make_run_script $TEMPLATE $MY_APP $MY_BARCODE $Chr
					make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $Chr 
					# submit this job to the scheduler
					if [ $USE_SLURM -eq 1 ]; then
						cd $TOP/logs
						echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
						sbatch $SLURM_OPT $MY_SLURM_SCRIPT
					# run the script locally
					else
						echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CALLER.$Chr.log"
						time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CALLER.$Chr.log &
					fi # end of USE_SLURM
				done
			else
				printf "\n\e[32m doing $MY_CALLER for $SLX.$MY_BARCODE\e[0m\n"
				#FASTQ_FILES=`ls $FASTQ_DIR/${SLX}.${MY_BARCODE}.*r_1.fq.gz | grep -v lost`
				#MY_CELL=`echo $FASTQ_FILES[0] | cut -d/ -f7 | cut -d. -f3`
				MY_APP=$SLX.$MY_BARCODE.$MY_CALLER.sh
				MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CALLER.sh

				make_run_script $TEMPLATE $MY_APP $MY_BARCODE $MY_CELL
				make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $MY_CELL

				# submit this job to the scheduler
				if [ $USE_SLURM -eq 1 ]; then
					cd $TOP/logs
					echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
					sbatch $SLURM_OPT $MY_SLURM_SCRIPT
				# run the script locally
				else
					echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_CALLER.log"
					time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_CALLER.log &
				fi # end of USE_SLURM
			fi # end of -n PER_CHR
		done # end of MY_BARCODE

	fi # end of #MY_CALLER
} # end of run_pipeline

################
## READ COUNT ##
################
#run_pipeline "per.run.count.read.base.fq" 0 # per SLX

###########################################
## Split FASTQ into several chunks       ##
## If SL2 in use (wall-time<36), no need ##
###########################################
##run_pipeline "per.barcode.fastq.split" 1               # per barcode
##run_per_lane_per_chunk_mapping $NUM_CHUNK

#############
## MAPPING ##
#############
#run_pipeline "per.lane.trim.bismark" 1 # per lane
#run_per_lane_mapping 0 # [depreciated; use above] equivalent with "run_per_lane_per_chunk_mapping 1"

#############
## MERGING ##
#############
## A) SINGLE RUN (e.g. SLX.8074)
#run_pipeline "per.barcode.merge.perlane.bam" 1  #run_pipeline "split.bam" 1 # for a patch (now part of merge.perlane.bam)
## B) MULTIPLE RUNS (e.g. SLX.8074.SLX8080)
#run_pipeline "merge.perrun.bam" 0

#############
## SORTING ##
#############
# run this only when bam files are splitted by chr
# this is to sort reads by name to run Bismark later
# and also add RG tag to run BisSNP later
#run_pipeline "sort.per.chr.bam" 0 

#############
## CALLING ##
#############
#run_pipeline "bismark.caller" 0        # per barcode
#run_pipeline "bismark.caller.by.chr" 0 # or run_pipeline "bismark.caller" 1

#run_pipeline "bissnp" 1             # per barcode
									# this is to make a RG.recal.bam to feed next step
#run_pipeline "bissnp.per.chr" 1 1
#run_pipeline "methylextract.per.chr" 1 1 #run_pipeline "methylextract" 0 # [deprecated]

###################
## COVERAGE STAT ##
###################
#run_pipeline "cpg.coverage.report" 0     # per barcode
#run_pipeline "non.cpg.coverage.report.by.chr" 0 1     # per chr CpH only 
