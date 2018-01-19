#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 2/Jun/2014
# Last modified: 2/Jun/2014 
# To run bedtools genomecov

PROJECT_DIR=$RESULT_DIR/$PROJECT 

	# CoverageBed
	mkdir_unless $PROJECT_DIR/Coverage/$Sample
	time bedtools genomecov \
			-ibam $TOPHAT_OUT/accepted_hits.bam > $PROJECT_DIR/Coverage/$Sample/$Sample.genomecov.txt
