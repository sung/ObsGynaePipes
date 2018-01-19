#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 20/Apr/2014
# Last modified: 20/Jun/2014 
# A wrapper script to run bis-snp 
# http://epigenome.usc.edu/publicationdata/bissnp2011/

source /whale-data/ssg29/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source /whale-data/ssg29/lib/sung.sh #defines 'mkdir_unless'

PROJECT='SGA_v3'
export PROJECT_DIR=$RESULT_DIR/$PROJECT 
BISSNP_T=$BIN_TOP/template/bissnp.template.sh

if [ !  -f $BISSNP_T ];then
	echo "$BISSNP_T not found"
	exit
fi

function make_bissnp_template(){
    MY_SLX=$1 # e.g. SLX-8080
    MY_BARCODE=$2 # e.g. A001 

    MY_RUN_SCRIPT=$BIN_TOP/script/$PROJECT.bissnp.${MY_SLX}.${MY_BARCODE}.sh
    cp $BISSNP_T $MY_RUN_SCRIPT

    sed -i "s/MY_SLX/$MY_SLX/" $MY_RUN_SCRIPT
    sed -i "s/MY_BARCODE/$MY_BARCODE/" $MY_RUN_SCRIPT
}

#ALL_BARCODES=(SLX-8074.A001) 

for MY_CODE in ${ALL_BARCODES[*]}
do
	MY_BARCODE=(${MY_CODE//./ }) # split string and make an array

	printf "\n\e[32mrunning $MY_CODE\e[0m\n"
	NOW=$(date +"%F-%H-%M")
	make_bissnp_template ${MY_BARCODE[0]} ${MY_BARCODE[1]}
	nice -n 10 time bash $BIN_TOP/script/$PROJECT.bissnp.${MY_CODE}.sh >& $TOP/logs/$PROJECT.bissnp.$MY_CODE.$NOW.log &
	printf "\n\e[32mLog=$TOP/logs/$PROJECT.bissnp.$MY_CODE.$NOW.log\e[0m\n"
done

