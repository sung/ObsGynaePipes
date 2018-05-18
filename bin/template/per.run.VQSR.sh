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

##########
## VQSR ##
##########
#https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
#https://software.broadinstitute.org/gatk/documentation/article.php?id=39
#https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php

My_VCF=$PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.vcf.gz 
if [ ! -s $My_VCF ];then
    echo -e "\e[031m$My_VCF not found\e[0m\n"
    exit
fi

###################
# 1. VQSR for SNP #
###################
	echo -e "java GATK -T VariantRecalibrator -input $My_VCF 2> $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_SNP.log\n"
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$HOME/scratch/temp $GATK_HOME/GenomeAnalysisTK.jar \
        -T VariantRecalibrator \
        -nt $NT \
		-R $GENOME \
        -input $My_VCF \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_HAPMAP \
        -resource:omni,known=false,training=true,truth=true,prior=12.0 $GATK_OMNI \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_1K \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
        -an DP \
        -an QD \
        -an FS \
        -an SOR \
        -an MQ \
        -an MQRankSum \
        -an ReadPosRankSum \
        -an InbreedingCoeff \
        -mode SNP \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        -recalFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_SNP.recal \
        -tranchesFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_SNP.tranches \
        -rscriptFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_SNP_plots.R \
        2> $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_SNP.log

########################
# 2. ApplyVQSR for SNP #
########################
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$HOME/scratch/temp $GATK_HOME/GenomeAnalysisTK.jar \
        -T ApplyRecalibration \
        -nt $NT \
		-R $GENOME \
        -input $My_VCF \
        -mode SNP \
        --ts_filter_level 99.0 \
        -recalFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_SNP.recal \
        -tranchesFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_SNP.tranches \
        -o ${My_VCF%.vcf.gz}.recalSNP.rawINDEL.vcf.gz
#####################
# 3. VQSR for Indel #
#####################
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$HOME/scratch/temp $GATK_HOME/GenomeAnalysisTK.jar \
        -T VariantRecalibrator \
        -nt $NT \
		-R $GENOME \
        -input ${My_VCF%.vcf.gz}.recalSNP.rawINDEL.vcf.gz \
        -resource:mills,known=false,training=true,truth=true,prior=12.0 $GATK_MILLS \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
        -an QD \
        -an DP \
        -an FS \
        -an SOR \
        -an MQRankSum \
        -an ReadPosRankSum \
        -an InbreedingCoeff \
        -mode INDEL \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        --maxGaussians 4 \
        -recalFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_INDEL.recal \
        -tranchesFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_INDEL.tranches \
        -rscriptFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_INDEL_plots.R \
        2> $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_INDEL.log
##########################
# 4. ApplyVQSR for Indel #
##########################
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$HOME/scratch/temp $GATK_HOME/GenomeAnalysisTK.jar \
        -T ApplyRecalibration \
        -nt $NT \
		-R $GENOME \
        -input ${My_VCF%.vcf.gz}.recalSNP.rawINDEL.vcf.gz \
        -mode INDEL \
        --ts_filter_level 99.0 \
        -recalFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_INDEL.recal \
        -tranchesFile $PROJECT_DIR/GATK/$SLX.bwa.samblaster.sorted.recal.HC.GenotypeGVCFs.recalibrate_INDEL.tranches \
        -o ${My_VCF%.vcf.gz}.recalSNP.recalINDEL.vcf.gz
echo -e "VQSR done for $SLX"
