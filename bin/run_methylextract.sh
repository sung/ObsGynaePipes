#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 22/may/2014
# Last modified: 22/May/2014 
# A wrapper script to run MethylExtract 
# http://bioinfo2.ugr.es/MethylExtract/Manual.html

source /whale-data/ssg29/Methyl-Seq/config/methyl_seq.config # export envrionment variables

ME_T=$BIN_TOP/template/methylextract.template.sh

if [ !  -f $ME_T ];then
	echo "$ME_T not found"
	exit
fi

function make_methylextract_template(){
    MY_BARCODE=$1 # e.g. A001 

    MY_RUN_SCRIPT=$BIN_TOP/script/$PROJECT.methylextract.${MY_BARCODE}.sh
    cp $ME_T $MY_RUN_SCRIPT

    sed -i "s/MY_BARCODE/$MY_BARCODE/" $MY_RUN_SCRIPT
}

ALL_BARCODES=(A011) 

for MY_BARCODE in ${ALL_BARCODES[*]}
do
	printf "\n\e[32mrunning $MY_BARCODE\e[0m\n"
	make_methylextract_template $MY_BARCODE
	time bash $BIN_TOP/script/$PROJECT.methylextract.${MY_BARCODE}.sh >& $TOP/logs/$PROJECT.methylextract.$MY_BARCODE.log & 
	#time bash $BIN_TOP/script/$PROJECT.methylextract.${MY_BARCODE}.sh 
done

