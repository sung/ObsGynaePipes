#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Sep/2015
# Last modified: 26/Nov/2015
# Optimised and customised to run at the Darwin HPC

###############################
# See bin/merge.junctions.sh ##
###############################

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT

mkdir_unless $PROJECT_DIR	

RAW_BEDS=$PROJECT_DIR/$PROJECT.merged.junctions.bed
RAW_JUNCS=$PROJECT_DIR/$PROJECT.merged.junctions.junc # also defined within the config file for project-wise

echo -e "merging juction file $RAW_JUNCS\n"

echo -e "cat $PROJECT_DIR/TopHat/*/junctions.bed > $RAW_BEDS\n"
time cat $PROJECT_DIR/TopHat/*/junctions.bed > $RAW_BEDS

echo -e "bed_to_juncs < $RAW_BEDS | sort -k 1,4 -u | sort -k 1,1 > $RAW_JUNCS\n"
time bed_to_juncs < $RAW_BEDS | sort -k 1,4 -u | sort -k 1,1 > $RAW_JUNCS

echo -e "gzip -9 -f $RAW_BEDS\n"
gzip -9 -f $RAW_BEDS

echo -e "ALL done\n"
