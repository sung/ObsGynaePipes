#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Jul/2014
# Last modified: 15/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
export Barcode="MY_BARCODE" # e.g. A001

	#######################
	## 13. MethylExtract ##
	#######################
	if [ $RUN_METHYL_EXT -eq 1 ]; then
		if [ -s $ME_WRAPPER ]; then
			mkdir_unless $PROJECT_DIR/MethylExtract
			mkdir_unless $PROJECT_DIR/MethylExtract/$Barcode
			source $ME_WRAPPER
		else
			echo -e "cannot find $ME_WRAPPER"
			exit
		fi
		echo -e "\e[32m$SLX.$Barcode done for MethylExtract \e[0m" 
	fi # end of $RUN_METHYL_EXT -eq 1 
