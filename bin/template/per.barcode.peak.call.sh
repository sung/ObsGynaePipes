#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 26/Feb/2015
# Last modified: 27/Feb/2015 
# Optimised and customised to run at the SBS and HPC machine

source $HOME/Pipelines/config/chip_seq.config 
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless', 'make_run_script'

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir_unless $PROJECT_DIR	

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Cell.$Barcode\e[0m\n"

BAM_FILES=`ls $PROJECT_DIR/Bowtie/$Barcode/$SLX.$Barcode.$Cell.*sorted.bam`
query_bam=$(printf " %s" "${BAM_FILES[@]}")

for My_control in ${CONTROL[*]}
do
	mkdir_unless $PROJECT_DIR/MACS2
	mkdir_unless $PROJECT_DIR/MACS2/$Barcode

	if [[ $Barcode == $My_control ]];then
		echo -e "\e[031mThis Barcode($Barcode) is a control\e[0m\n"
		exit
	else
		CTRL_BAM_FILES=`ls $PROJECT_DIR/Bowtie/$My_control/$SLX.$My_control.$Cell.*sorted.bam`
		ctrl_bam=$(printf " %s" "${CTRL_BAM_FILES[@]}")
		MACS_PREFIX=$SLX.$Barcode.$My_control

		if [ $RUN_MACS_CALLPEAK -eq 1 ];then
			echo -e "macs2 callpeak --treatment $query_bam -c $ctrl_bam -g $MACS_GENOME -n $MACS_PREFIX"
			time macs2 callpeak \
				--treatment $query_bam \
				--control $ctrl_bam \
				--gsize $MACS_GENOME \
				--nomodel \
				--extsize $MACS_SHIFT_SIZE \
				--call-summits \
				--bdg \
				--SPMR \
				--outdir $PROJECT_DIR/MACS2/$Barcode \
				--name $SLX.$Barcode.$My_control 2> $PROJECT_DIR/MACS2/$Barcode/$SLX.$Barcode.$My_control.macs2.log 
				#--broad (not for H3K4me3)
		fi

		if [ $RUN_MACS_BDGCMP -eq 1 ];then
			MACS_OUT=$PROJECT_DIR/MACS2/$Barcode/$MACS_PREFIX"_treat_pileup.bdg"
			MACS_CTRL_OUT=$PROJECT_DIR/MACS2/$Barcode/$MACS_PREFIX"_control_lambda.bdg"

			if [ -s $MACS_OUT ] && [ -s $MACS_CTRL_OUT ]; then
				echo -e "macs2 bdgcmp -m FE"
				time macs2 bdgcmp \
					--treatment $MACS_OUT \
					--control $MACS_CTRL_OUT \
					-o $PROJECT_DIR/MACS2/$Barcode/$SLX.$Barcode.$My_control.FE.bdg \
					-m FE

				echo -e "macs2 bdgcmp -m logLR"
				time macs2 bdgcmp \
					--treatment $MACS_OUT \
					--control $MACS_CTRL_OUT \
					-o $PROJECT_DIR/MACS2/$Barcode/$SLX.$Barcode.$My_control.logLR.bdg \
					-m logLR \
					--pvalue $MACS_PVALUE
			else
				echo -e "\e[031m$MACS_OUT or $MACS_CTRL_OUT not found\e[0m\n"
				exit
			fi
		fi
	fi
done

echo -e "\e[32m$Barcode done for per.barcode.peak.call\e[0m" 
