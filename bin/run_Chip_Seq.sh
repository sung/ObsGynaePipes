#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 26/Feb/2015
# Last modified: 26/Feb/2015 
# Optimised and customised to run at the Darwin HPC

export SL='SL2' # SL2 (paid ad-hoc), SL3 (free)
source $HOME/Pipelines/config/chip_seq.config 
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless', 'make_run_script'

############
## Claire ##
############
SLX='SLX-8766'

VERSION=$SPECIES."v1" 

FASTQ_DIR=$HOME/data/fastq/$SLX
if [ ! -d $FASTQ_DIR ];then
	echo "$FASTQ_DIR not found"
	exit
fi

ALL_BARCODES=( $(for i in `ls $FASTQ_DIR/$SLX.*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done | uniq) )
#ALL_BARCODES=(A002 A005 A007 A012 A013) # Tov21G
#ALL_BARCODES=(A014 A015 A018 A004 A006) # Jhoc5

mkdir_unless $RESULT_DIR	
mkdir_unless $TOP/logs
mkdir_unless $BIN_TOP/script
mkdir_unless $BIN_TOP/slurm

function run_pipeline(){
	MY_CALLER=$1 # 
	USE_SLURM=$2 # 1 (use slurm) or 0 (do not use slurm)
	PER_CHR=$3 # 1 (a flag to run by chromosome-wise)

	TEMPLATE=$BIN_TOP/template/$MY_CALLER.sh # template/tophat.persample.sh # per sample 
	if [ ! -f $TEMPLATE ];then
		echo "$TEMPLATE not found"
		exit
	fi

	printf "\e[32mDoing $MY_CALLER for $SLX\e[0m\n"
	#############################
	# per SLX (count read.base) #
	#############################
	if [[ $MY_CALLER == 'per.run.count.read.base.fq' ]]; then
		MY_APP=$SLX.$VERSION.$MY_CALLER.sh
		make_run_script $TEMPLATE $MY_APP

		# run the script locally
		echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$VERSION.$MY_CALLER.log"
		time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$VERSION.$MY_CALLER.log &
	############
	# per lane #
	############
	elif [[ $MY_CALLER =~ 'per.lane' ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\n\e[32m$MY_CALLER for $SLX.$MY_BARCODE\e[0m\n"
			#/home/ssg29/data/fastq/SLX-9338/SLX-9338.HAL21.000000000-ACU88.s_1.r_1.fq.gz
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
				MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CALLER.sh
				make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $MY_CELL
				cd $TOP/logs
				echo -e "sbatch $MY_SLURM_SCRIPT"
				sbatch $MY_SLURM_SCRIPT
			# run the script locally
			else
				echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CALLER.log"
				time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CALLER.log &
			fi
		done # end of MY_BARCODE
	else
		echo "Not supported -- either 'count.read.base' or per.lane or per.barcode\n"
		exit
	fi # end of if MY_CALLER
} # end of run_methyl_call 

#run_pipeline "per.run.count.read.base.fq" 0 # per SLX

run_pipeline "per.lane.trim.bowtie" 1 # per lane 

#run_pipeline "per.barcode.peak.call" 0 # per barcode (sample) 
