#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 17/Apr/2014
# Last modified: 8/Jun/2014 
# A wrapper script to run trim_galore and bismark 
# Optimised and customised to run at the Darwin HPC

#source $HOME/Methyl-Seq/config/sample.SLX-8077.config # export envrionment variables
source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

SLX="SLX-8075"
ALL_BARCODES=(A008 A009) 
VERSION="v1" 

MAPPING_T=$BIN_TOP/template/trim_bismark_template_hpc.v1.sh
SLURM_T=$BIN_TOP/template/slurm_template.v1.sh

export Other_bam='' # if RUN_MERGE_BAM=1, specify your bam file location

if [ ! -f $MAPPING_T ];then
	echo "$MAPPING_T not found"
	exit
fi
if [ ! -f $MERGE_BAM_T ];then
	echo "$MERGE_BAM_T not found"
	exit
fi
if [ ! -f $SLURM_T ];then
	echo "$SLURM_T not found"
	exit
fi

mkdir_unless $RESULT_DIR	
mkdir_unless $TOP/logs
mkdir_unless $BIN_TOP/script
mkdir_unless $BIN_TOP/slurm

function make_run_script(){
    MY_TEMPLATE=$1 # e.g. $MAPPING_T
    MY_APP=$2 # e.g. SLX-8077.trim.bismark.A005.C48CWACXX.s_6.sh
    MY_BARCODE=$3 # e.g. A005 
    MY_CELL=$4 # e.g. C48CWACXX (optional)
    MY_LANE=$5 # e.g. s_1 (optional)

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
}

function make_slurm_script(){
    MY_SLURM_SCRIPT=$1 # e.g. $BIN_TOP/slurm/slurm.SLX-8077.A005.C48CWACXX.s_5.sh
    MY_APP=$2 # e.g. SLX-8077.trim.bismark.A005.C48CWACXX.s_6.sh
    MY_STEP=$3 # e.g. 'mapping', 'merging', 'bismark.caller', 'methylextract', 'bissnp', 'coverage'
    MY_BARCODE=$4 # e.g. A005 
    MY_CELL=$5 # e.g. C48CWACXX (optional)
    MY_LANE=$6 # e.g. s_1  (optional)

    cp $SLURM_T $MY_SLURM_SCRIPT

    sed -i "s/MY_APP/$MY_APP/" $MY_SLURM_SCRIPT
    sed -i "s/MY_STEP/$MY_STEP/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLX/$SLX/" $MY_SLURM_SCRIPT
    sed -i "s/MY_BARCODE/$MY_BARCODE/" $MY_SLURM_SCRIPT
	if [ -n "$MY_CELL" ];then #string is not null.
    	sed -i "s/MY_CELL/$MY_CELL/" $MY_SLURM_SCRIPT
	fi
	if [ -n "$MY_LANE" ];then #string is not null.
    	sed -i "s/MY_LANE/$MY_LANE/" $MY_SLURM_SCRIPT
	fi

    sed -i "s/MY_SLURM_ACNT/$SLURM_ACNT/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_NODE_CNT/$SLURM_NODE_CNT/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_NT/$SLURM_NODE_NT/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_TIME/$SLURM_TIME/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_MAIL_TYPE/$SLURM_MAIL_TYPE/" $MY_SLURM_SCRIPT
}

# template/trim_bismark_template_hpc.v1.sh (fastqc, trim, mapping <- per lane)
function run_per_lane_mapping(){
	for MY_BARCODE in ${ALL_BARCODES[*]}
	do
		printf "\n\e[32mmapping $SLX.$MY_BARCODE\e[0m\n"

		FASTQ_FILES=`ls $FASTQ_DIR/${SLX}.${MY_BARCODE}*.fq.gz | grep -v lost`
		for FastQ_file in $FASTQ_FILES	# per fastq file
		do
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

		#for i in `ls $FASTQ_DIR/${SLX}.${MY_BARCODE}*r_1.fq.gz`; do echo $i | cut -d/ -f6 | cut -d. -f3,4 ; done
		for MY_CELL_LANE in ${Cell_lane_array[*]} 
		do
			# fixme
			# skip for SLX-8077.A005 lane 1,2,3,4
			if [ $MY_CELL_LANE = "C48CWACXX.s_1" ] || [ $MY_CELL_LANE = "C48CWACXX.s_2" ] || [ $MY_CELL_LANE = "C48CWACXX.s_3" ] || [ $MY_CELL_LANE = "C48CWACXX.s_4" ]; then 
				echo skipping lane $MY_CELL_LANE
				continue
			fi

			echo -e "\t\e[32mCELL.LANE=$MY_CELL_LANE\e[0m" 

			MY_CELL=`echo $MY_CELL_LANE | cut -d. -f1`
			MY_LANE=`echo $MY_CELL_LANE | cut -d. -f2`
			MY_APP=$SLX.trim.bismark.$MY_BARCODE.$MY_CELL.$MY_LANE.sh

			make_run_script $MAPPING_T $MY_APP $MY_BARCODE $MY_CELL $MY_LANE

    		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.mapping.$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.sh
			make_slurm_script $MY_SLURM_SCRIPT $MY_APP "mapping" $MY_BARCODE $MY_CELL $MY_LANE

			# submit this job
			cd $TOP/logs
			echo -e "sbatch $MY_SLURM_SCRIPT"
			sbatch $MY_SLURM_SCRIPT
		done # end of MY_CELL_LANE

	done # end of MY_BARCODE
} # end of run_per_lane_mapping

# template/merge_per_lane_bam.sh 
MERGE_BAM_T=$BIN_TOP/template/merge_per_lane_bam.sh
function run_merge_bam(){
	for MY_BARCODE in ${ALL_BARCODES[*]}
	do
		printf "\n\e[32mmerging $SLX.$MY_BARCODE\e[0m\n"
		MY_APP=$SLX.merge.bam.$MY_BARCODE.sh
		make_run_script $MERGE_BAM_T $MY_APP $MY_BARCODE

		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.merging.$SLX.$MY_BARCODE.sh
		make_slurm_script $MY_SLURM_SCRIPT $MY_APP "merging" $MY_BARCODE 

		#fixme
		# submit this job
		#cd $TOP/logs
		#echo -e "sbatch $MY_SLURM_SCRIPT"
		#sbatch $MY_SLURM_SCRIPT

		echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.merging.sh"
		time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.merging.sh &
	done # end of MY_BARCODE
} # end of run_merge_bam

# 3. template/methyl.call.sh
function run_methyl_call(){
	MY_CALLER=$1
	METHYL_CALL_T=$BIN_TOP/template/$MY_CALLER.sh
	if [ ! -f $METHYL_CALL_T ];then
		echo "$METHYL_CALL_T not found"
		exit
	fi

	for MY_BARCODE in ${ALL_BARCODES[*]}
	do
		printf "\n\e[32mcalling CpG with $MY_CALLER from $SLX.$MY_BARCODE\e[0m\n"

		MY_APP=$SLX.$MY_CALLER.$MY_BARCODE.sh
		make_run_script $METHYL_CALL_T $MY_APP $MY_BARCODE

		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$MY_CALLER..$SLX.$MY_BARCODE.sh
		make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE 

		if [ $MY_CALLER = 'bismark.caller' ]; then
			echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CALLER.log"
			time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CALLER.log &
		elif [ $MY_CALLER = 'methylextract' ] || [ $MY_CALLER = 'bissnp' ]; then
			# submit this job
			cd $TOP/logs
			echo -e "sbatch $MY_SLURM_SCRIPT"
			sbatch $MY_SLURM_SCRIPT
		else
			echo -e "caller should one of 1)bismark.caller, 2)methylextract, or 3)bissnp"
			exit
		fi
	done # end of MY_BARCODE
} # end of run_methyl_call 

# 4. template/coverage.report.sh
COVERAGE_T=$BIN_TOP/template/coverage.report.sh
function run_coverage_report(){
	for MY_BARCODE in ${ALL_BARCODES[*]}
	do
		printf "\n\e[32mDepth of Cov $SLX.$MY_BARCODE\e[0m\n"
		MY_APP=$SLX.coverage.report.$MY_BARCODE.sh
		make_run_script $COVERAGE_T $MY_APP $MY_BARCODE

		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.coverage.$SLX.$MY_BARCODE.sh
		make_slurm_script $MY_SLURM_SCRIPT $MY_APP "coverage" $MY_BARCODE 

		# fixme
		# submit this job
		#cd $TOP/logs
		#echo -e "sbatch $MY_SLURM_SCRIPT"
		#sbatch $MY_SLURM_SCRIPT

		echo -e "bash $BIN_TOP/script/$MY_APP"
		time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.coverage.log &
	done # end of MY_BARCODE
} # end of run_coverage_report 

#run_per_lane_mapping 
run_merge_bam
#run_methyl_call "bismark.caller"
#run_methyl_call "methylextract"
#run_methyl_call "bissnp"
#run_coverage_report

<<fastq
SLX-8077.A005.C48CWACXX.s_1.r_1.fq.gz
SLX-8077.A005.C48CWACXX.s_1.r_2.fq.gz
SLX-8077.A005.C48CWACXX.s_2.r_1.fq.gz
SLX-8077.A005.C48CWACXX.s_2.r_2.fq.gz
fastq

