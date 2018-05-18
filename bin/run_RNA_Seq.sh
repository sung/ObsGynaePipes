#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Jul/2014
# Last modified: 19/Jan/2016
# a Wrapper script to run tophat and cufflink
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # calls global.config, slurm.config, samples.config 
source $HOME/Pipelines/lib/sung.sh # defines user defined functions (e.g. make_run_script)
# 'make_slurm_script' defined within config/slurm.config # defines 'make_slurm_script'

# $CASECTR_PAIRS defined within config/samples.config
<<DEGONLY 
#edgeR.complex.sh & edgeR.simple.sh
for MY_CASECTR in ${CASECTR_PAIRS[*]}
do
	dummy_array=(${MY_CASECTR//,/ }) # (D709_D501,D709_D502,D101_D102) => (D709_D501 D709_D502 D101_D102)
	for dummy_element in ${dummy_array[*]}
	do
		ALL_BARCODES+=($dummy_element) #push it
	done
done
DEGONLY

mkdir -p $RESULT_DIR
mkdir -p $TOP/logs
mkdir -p $BIN_TOP/script
mkdir -p $BIN_TOP/slurm

# this function is also defined within lib/sung.sh
function run_pipeline(){
	MY_CALLER=$1 # 
	USE_SLURM=$2 # 1 (use slurm) or 0 (do not use slurm)

	TEMPLATE=$BIN_TOP/template/$MY_CALLER.sh # template/tophat.persample.sh # per sample 
	if [ ! -f $TEMPLATE ];then
		echo "$TEMPLATE not found"
		exit
	fi

	printf "\e[32mdoing $MY_CALLER for $SLX\e[0m\n"
	#################################
	# per SLX run (count read.base) #
	#################################
	if [[ $MY_CALLER =~ ^per.run ]]; then
		MY_APP=$SLX.$VERSION.$COHORT.$MY_CALLER.sh
		make_run_script $TEMPLATE $MY_APP

		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$VERSION.$COHORT.$MY_CALLER.sh
		make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER

		# submit this job to the scheduler
		if [ $USE_SLURM -eq 1 ]; then
			cd $TOP/logs
			echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
			sbatch $SLURM_OPT $MY_SLURM_SCRIPT
		else
			# run the script locally
			echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.$COHORT.log"
			time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.$COHORT.log &
		fi
	###########################
	# per sample (or barcode) #
	###########################
	elif [[ $MY_CALLER =~ ^per.sample ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\e[32mdoing $SLX.$VERSION.$MY_BARCODE\e[0m\n"
			#FASTQ_FILES=`ls $HOME/data/fastq/$SLX/${SLX}.${MY_BARCODE}.*r_1.fq.gz | grep -v lost`
			#MY_CELL=`echo $FASTQ_FILES[0] | cut -d/ -f7 | cut -d. -f3`
            MY_CELL="XXX"
			MY_APP=$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.$MY_CELL.sh
			make_run_script $TEMPLATE $MY_APP $MY_BARCODE $MY_CELL

			# submit this job to the scheduler
			if [ $USE_SLURM -eq 1 ]; then
				MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.$MY_CELL.sh
				make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $MY_CELL
				cd $TOP/logs
				echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
				sbatch $SLURM_OPT $MY_SLURM_SCRIPT
			# run the script locally
			else
				echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.log"
				time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.log &
			fi
		done # end of MY_BARCODE
	############
	# per lane #
	############
	elif [[ $MY_CALLER =~ '^per.lane' ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\n\e[32m$MY_CALLER for $SLX.$MY_BARCODE\e[0m\n"
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
	#####################
	# per batch (edgeR) #
	#####################
	elif [[ $MY_CALLER =~ ^edgeR ]]; then
		MY_APP=$SLX.$VERSION.$MY_CALLER.sh
		make_run_script $TEMPLATE $MY_APP

		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$VERSION.$MY_CALLER.sh
		make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER

		# submit this job to the scheduler
		if [ $USE_SLURM -eq 1 ]; then
			cd $TOP/logs
			echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
			sbatch $SLURM_OPT $MY_SLURM_SCRIPT
		# run the script locally
		else
			echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.log"
			time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.log &
		fi
	else
		echo -e "\e[031m$MY_CALLER not found\e[0m\n"
		exit
	fi # end of if MY_CALLER
} # end of run_pipeline

######################################################################################
# 1. build the transcriptome index data files for the given annotation and then exit #
# once and for all
######################################################################################
if [ $RUN_TH_TRINDEX -eq 1 ]; then
	mkdir -p $TR_INDEX_DIR
	echo -e "\ntophat2 -p 8 -G $GTF --transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX $BOWTIE2_INDEX_BASE"
	time tophat2 -p 8 -G $GTF --transcriptome-index=$TR_INDEX_DIR/$TR_PREFIX $BOWTIE2_INDEX_BASE
fi
if [ $RUN_STAR_INDEX -eq 1 ]; then
	mkdir -p $STAR_INDEX_DIR
	echo -e "STAR --runMode genomeGenerate\n"
	time STAR --runMode genomeGenerate \
		--genomeDir $STAR_INDEX_DIR \
		--genomeFastaFiles $STAR_GENOME.gz \
		--sjdbGTFfile $STAR_GTF \
		--runThreadN 4
fi

#run_pipeline "per.run.count.read.base.fq" 0 # per SLX

#####################
## for m/total-RNA ##
#####################
#run_pipeline "per.sample.samtools.cram" 1 # per barcode (sample)
#run_pipeline "per.sample.tophat" 1 # per barcode (sample)
#run_pipeline "per.run.merge.junctions" 0 # per run 
#run_pipeline "per.sample.tophat.unmapped" 1 # per barcode (sample)
#run_pipeline "per.sample.bwa.rna.seq" 1 # per barcode (sample) # this is a pre-processing of CIRI2
run_pipeline "per.sample.CIRI" 1 # per barcode (sample) # this is a pre-processing of CIRI2
#run_pipeline "per.sample.salmon" 1 # per barcode (sample)
#run_pipeline "per.sample.filter.human" 1 # per barcode (sample) for Marcus
#run_pipeline "per.sample.dupRadar" 1 # per barcode (sample)
#run_pipeline "per.sample.dupRadar.9" 1 # per barcode (sample)

###################################
## Below two for variant calling ##
###################################
#run_pipeline "per.sample.star2pass" 0 # per barcode (sample)
#run_pipeline "per.sample.gatk.rnaseq" 1 # per barcode (sample)


###################
## for small-RNA ##
###################
#run_pipeline "per.sample.mirdeep2" 1 # per barcode (sample)
#run_pipeline "per.sample.seqcluster" 1 # deprecated (use mirdeep2)

###############
## STRINGTIE ##
###############
#run_pipeline "per.sample.stringtie" 1 # per barcode (sample) (< 1hr per-sample)

##############
## CUFFLINK ##
##############
# alternative transcriptome assembly method of stringtie
#run_pipeline "per.sample.cufflink" 1 # per barcode (sample)

#################
## CUFFCOMPARE ##
#################
#merge individual assembly into one
#run_pipeline "per.run.cuffcompare" 0 # per run (either 'case', 'control' or 'all')

#####################
## STRINGTIE-MERGE ##
#####################
# alternative of cuffcompare above
# not implemented yet
#run_pipeline "per.run.stringtie.merge" 0 # per run (either 'case' or 'control')

##############
## BEDGRAPH ##
##############
#run_pipeline "per.run.bedgraph" 1 # per run (e.g. 'SGA', 'PET', 'CTL', 'CTLM', 'CTLF')

## DEPRECATED BELOW
## NO NOT USE
#########
## DEG ##
#########
#run_pipeline "edgeR" 0
##run_pipeline "edgeR.simple" 0
##run_pipeline "edgeR.complex" 0
