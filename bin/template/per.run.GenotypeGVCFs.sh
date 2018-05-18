#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 9/Apr/2018
# Last modified: 9/Apr/2018
# Optimised and customised to run at the UIS CSD3 peta4 
# Based on GATK3

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/Pipelines/lib/sung.sh # defines user defined functions (e.g. make_run_script)
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT

COHORT="MY_COHORT" # CTLF, CTLM, POPS, CTL, SGA, PET
echo -e "COHORT=$COHORT"

mkdir -p $PROJECT_DIR	
mkdir -p $PROJECT_DIR/GATK

################
## GVCF Files ##
################
ls $PROJECT_DIR/GATK/*/*.HC.g.vcf.gz > $PROJECT_DIR/GATK/HC.g.vcf.list 

if [ -s $PROJECT_DIR/GATK/HC.g.vcf.list ];then
    My_VCF=$PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.vcf.gz 
	echo -e "java $GATK_HOME/GenomeAnalysisTK.jar -T GenotypeGVCFs --variant $PROJECT_DIR/GATK/HC.g.vcf.list -o $My_VCF 2> ${My_VCF%.gz}.log\n" 
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$HOME/scratch/temp $GATK_HOME/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -nt $NT \
		-R $GENOME \
        --dbsnp $DBSNP \
        --variant $PROJECT_DIR/GATK/HC.g.vcf.list \
        -o $My_VCF 2> ${My_VCF%.gz}.log
	echo -e "\e[032mDone: GATK GenotypeGVCFs\e[0m\n"
else
	echo -e "\e[031$PROJECT_DIR/GATK/HC.g.vcf.list not found\e[0m\n"
	exit

fi

echo -e "GenotypeGVCFs done for $SLX"
