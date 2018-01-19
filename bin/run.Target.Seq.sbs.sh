#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 30/Jan/2015
# Last modified: 12/Feb/2015
# A wrapper script to run BWA and GATK for targetted resequencing
# Optimised and customised to run at the SBS machines (Jenny's)

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/Pipelines/lib/sung.sh

## initialise 'ALL_BARCODES'
if [ $SLX == 'SLX-9338' ]; then
	if [ $IS_PE -eq 1 ]; then
		ALL_BARCODES=( $(for i in `ls $FASTQ_DIR/$SLX.*.r_1.fq.gz | grep -v lost | grep -v HAL48`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done) )
	else
		ALL_BARCODES=( $(for i in `ls $FASTQ_DIR/$SLX.*.fq.gz | grep -v lost | grep -v HAL48`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done) )
	fi
else
	if [ $IS_PE -eq 1 ]; then
		ALL_BARCODES=( $(for i in `ls $FASTQ_DIR/$SLX.*.r_1.fq.gz | grep -v lost | grep -v HALO16`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done) )
	else
		ALL_BARCODES=( $(for i in `ls $FASTQ_DIR/$SLX.*.fq.gz | grep -v lost | grep -v HALO16`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done) )
	fi
fi

ALL_BARCODES=(HALO1)

#'mkdir_unless' defined within $HOME/lib/sung.sh
mkdir_unless $RESULT_DIR	
mkdir_unless $TOP/logs
mkdir_unless $BIN_TOP/script
mkdir_unless $BIN_TOP/slurm

function make_run_script(){
    MY_TEMPLATE=$1 # e.g. template/tophat.persample or template/edgeR 
    MY_APP=$2 # e.g. SLX-8547.D704_D503.tophat.persample.sh

    MY_BARCODE=$3 # e.g. D704_D503 (optional)
    MY_CELL=$4 # e.g. C48CWACXX (optional)
    MY_LANE=$5 # e.g. s_1 (optional)
    MY_CHUNK=$6 # e.g. 1 (optional)

    cp $MY_TEMPLATE $BIN_TOP/script/$MY_APP

    sed -i "s/MY_SLX/$SLX/" $BIN_TOP/script/$MY_APP
    sed -i "s/MY_VERSION/$VERSION/" $BIN_TOP/script/$MY_APP
	if [ -n "$MY_BARCODE" ];then #string is not null.
    	sed -i "s/MY_BARCODE/$MY_BARCODE/" $BIN_TOP/script/$MY_APP
	fi
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

function make_slurm_script(){
    MY_SLURM_SCRIPT=$1 # e.g. $BIN_TOP/slurm/slurm.SLX-8077.A005.C48CWACXX.s_5.sh
    MY_APP=$2 # e.g. SLX-8077.trim.bismark.A005.C48CWACXX.s_6.sh
    MY_STEP=$3 # e.g. 'mapping', 'merging', 'bismark.caller', 'methylextract', 'bissnp', 'coverage'

    MY_BARCODE=$4 # e.g. A005 
    MY_CELL=$5 # e.g. C48CWACXX or chr1 (optional)
    MY_LANE=$6 # e.g. s_1  (optional)
    MY_CHUNK=$7 # e.g. 1  (optional)

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
	if [ -n "$MY_CHUNK" ];then #string is not null.
    	sed -i "s/MY_CHUNK/$MY_CHUNK/" $MY_SLURM_SCRIPT
	fi

    sed -i "s/MY_SLURM_ACNT/$SLURM_ACNT/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_NODE_CNT/$SLURM_NODE_CNT/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_NT/$SLURM_NODE_NT/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_TIME/$SLURM_TIME/" $MY_SLURM_SCRIPT
    sed -i "s/MY_SLURM_MAIL_TYPE/$SLURM_MAIL_TYPE/" $MY_SLURM_SCRIPT
}

function run_pipeline(){
	MY_CALLER=$1 # 
	USE_SLURM=$2 # 1 (use slurm) or 0 (do not use slurm)
	PER_CHR=$3 # 1 (a flag to run by chromosome-wise)

	TEMPLATE=$BIN_TOP/template/$MY_CALLER.sh # template/tophat.persample.sh # per sample 
	if [ ! -f $TEMPLATE ];then
		echo "$TEMPLATE not found"
		exit
	fi

	printf "\e[32mdoing $MY_CALLER for $SLX\e[0m\n"
	#############################
	# per SLX (count read.base) #
	#############################
	if [[ $MY_CALLER == 'per.run.count.read.base.fq' ]]; then
		MY_APP=$SLX.$VERSION.$MY_CALLER.sh
		make_run_script $TEMPLATE $MY_APP

		# run the script locally
		echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$VERSION.$MY_CALLER.log"
		nice -n 10 time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$VERSION.$MY_CALLER.log &
	############
	# per lane #
	############
	elif [[ $MY_CALLER =~ 'per.lane' ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\n\e[32m$MY_CALLER for $SLX.$MY_BARCODE\e[0m\n"
			#/disk1/ssg29/data/fastq/SLX-9338/SLX-9338.HAL21.000000000-ACU88.s_1.r_1.fq.gz
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
					echo -e "sbatch $MY_SLURM_SCRIPT"
					sbatch $MY_SLURM_SCRIPT
				# run the script locally
				else
					echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.log"
					nice -n 10 time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.log &
				fi # end of USE_SLURM
			done # end of MY_CELL_LANE
		done
	###########################
	# per sample (or barcode) #
	###########################
	elif [[ $MY_CALLER =~ 'per.barcode' ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\e[32mdoing $SLX.$MY_BARCODE\e[0m\n"
			FASTQ_FILES=`ls $FASTQ_DIR/${SLX}.${MY_BARCODE}.*r_1.fq.gz | grep -v lost`
			MY_CELL=`echo $FASTQ_FILES[0] | cut -d/ -f7 | cut -d. -f3`
			MY_APP=$SLX.$MY_BARCODE.$MY_CALLER.$MY_CELL.sh
			make_run_script $TEMPLATE $MY_APP $MY_BARCODE $MY_CELL

			# submit this job to the scheduler
			if [ $USE_SLURM -eq 1 ]; then
				MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CALLER.$MY_CELL.sh
				make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $MY_CELL
				cd $TOP/logs
				echo -e "sbatch $MY_SLURM_SCRIPT"
				sbatch $MY_SLURM_SCRIPT
			# run the script locally
			else
				echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CALLER.log"
				nice -n 10 time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.log &
			fi
		done # end of MY_BARCODE
	else
		echo "Not supported -- either 'count.read.base' or per.lane or per.barcode\n"
		exit
	fi # end of if MY_CALLER
} # end of run_methyl_call 

#run_pipeline "per.run.count.read.base.fq" 0 # per SLX
#run_pipeline "per.lane.trim.bwa" 0 # per lane
#run_pipeline "per.barcode.gatk" 0 # per barcode (sample)
#run_pipeline "per.barcode.coverage" 0 # per barcode (sample)
#run_pipeline "per.barcode.gatk.somatic" 0 # per barcode (sample)
#run_pipeline "per.barcode.mutect" 0 # per barcode (sample)
run_pipeline "per.barcode.somatic.snp" 0 # per barcode (sample)
