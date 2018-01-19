#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Jul/2014
# Last modified: 15/Jul/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
export Barcode="MY_BARCODE" # e.g. A001
export Cell="MY_CELL" # e.g. C48CWACXX  
export Lane="MY_LANE" # e.g. s_1 


################
## 14 BIS-SNP ##
################
if [ $RUN_BIS_SNP -eq 1 ]; then
	if [ -s $BIS_SNP_WRAPPER ]; then
		mkdir_unless $PROJECT_DIR/BisSNP
		mkdir_unless $PROJECT_DIR/BisSNP/$Barcode
		source $BIS_SNP_WRAPPER
	else
		echo -e "cannot find $BIS_SNP_WRAPPER"
		exit
	fi
	echo -e "\e[32m$SLX.$Barcode done for CpG calling by Bis-SNP\e[0m" 
fi

if [ $RUN_BIS_SNP_NON_CPG -eq 1 ]; then
	if [ -s $BIS_SNP_NON_CPG_WRAPPER ]; then
		mkdir_unless $PROJECT_DIR/BisSNP
		mkdir_unless $PROJECT_DIR/BisSNP/$Barcode
		source $BIS_SNP_NON_CPG_WRAPPER
	else
		echo -e "cannot find $BIS_SNP_NON_CPG_WRAPPER"
		exit
	fi
	echo -e "\e[32m$SLX.$Barcode done for Non-CpG calling by Bis-SNP\e[0m" 
fi
