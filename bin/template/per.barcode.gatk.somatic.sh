#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 25/Mar/2015
# Last modified: 25/Mar/2015 

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/lib/sung.sh

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

#Control Barcode
My_control=`grep -P "$SLX\t$Barcode\t" $NPOOL | cut -f3`
if [ ! $My_control ]; then
	echo -e "\e[31mNo control Barcode found: $Barcode must be control\e[0m" 
	exit
else
	#################
	## Somatic SNP ##
	#################
	if [ $RUN_GATK_SOMATIC_SNP -eq 1 ]; then
		My_UG_SNP_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.snp.filtered.vcf.gz 
		My_UG_CTR_SNP_VCF=$PROJECT_DIR/GATK/$My_control/$SLX.$My_control.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.snp.filtered.vcf.gz 

		My_HC_SNP_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.snp.filtered.vcf.gz
		My_HC_CTR_SNP_VCF=$PROJECT_DIR/GATK/$My_control/$SLX.$My_control.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.snp.filtered.vcf.gz 

		if [ -s $My_UG_SNP_VCF ] && [ -s $My_UG_CTR_SNP_VCF ] && [ -s $My_HC_SNP_VCF ] && [ -s $My_HC_CTR_SNP_VCF ]; then
			echo -e "\e[32m$Barcode done for somatic snp\e[0m" 
			My_UG_SOMATIC_SNP_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.somatic.snp.vcf.gz 
			My_HC_SOMATIC_SNP_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.somatic.snp.vcf.gz 
			My_SOMATIC_SNP_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.somatic.snp.vcf.gz 
			#somatic mutations by substracting varaints of normal sample from that of tumor sample
			time bcftools isec -n-1 -w1 -o $My_UG_SOMATIC_SNP_VCF -O z $My_UG_SNP_VCF $My_UG_CTR_SNP_VCF 
			time bcftools isec -n-1 -w1 -o $My_HC_SOMATIC_SNP_VCF -O z $My_HC_SNP_VCF $My_HC_CTR_SNP_VCF 

			#PASS only
			time bcftools view -fPASS -o $My_UG_SOMATIC_SNP_VCF -O z $My_UG_SOMATIC_SNP_VCF 
			time bcftools view -fPASS -o $My_HC_SOMATIC_SNP_VCF -O z $My_HC_SOMATIC_SNP_VCF 

			#Index it
			time tabix -f -p vcf $My_UG_SOMATIC_SNP_VCF
			time tabix -f -p vcf $My_HC_SOMATIC_SNP_VCF

			#Concatenate (unify) HC and UG
			time bcftools concat -a -D -o $My_SOMATIC_SNP_VCF -O z $My_HC_SOMATIC_SNP_VCF $My_UG_SOMATIC_SNP_VCF
			time tabix -f -p vcf $My_SOMATIC_SNP_VCF
		else
			echo -e "\e[031m$My_UG_SNP_VCF or $My_UG_CTR_SNP_VCF or $My_HC_SNP_VCF or $My_HC_CTR_SNP_VCF not found\e[0m\n"
			exit
		fi
	fi
	###################
	## Somatic Indel ##
	###################
	if [ $RUN_GATK_SOMATIC_Indel -eq 1 ]; then
		My_UG_Indel_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.indel.filtered.vcf.gz 
		My_UG_CTR_Indel_VCF=$PROJECT_DIR/GATK/$My_control/$SLX.$My_control.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.indel.filtered.vcf.gz 

		My_HC_Indel_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.indel.filtered.vcf.gz
		My_HC_CTR_Indel_VCF=$PROJECT_DIR/GATK/$My_control/$SLX.$My_control.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.indel.filtered.vcf.gz 

		if [ -s $My_UG_Indel_VCF ] && [ -s $My_UG_CTR_Indel_VCF ] && [ -s $My_HC_Indel_VCF ] && [ -s $My_HC_CTR_Indel_VCF ]; then
			echo -e "\e[32m$Barcode done for somatic indel\e[0m" 
			My_UG_SOMATIC_Indel_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.somatic.indel.vcf.gz 
			My_HC_SOMATIC_Indel_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.somatic.indel.vcf.gz 
			My_SOMATIC_Indel_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.somatic.indel.vcf.gz 
			#somatic mutations by substracting varaints of normal sample from that of tumor sample
			time bcftools isec -n-1 -w1 -o $My_UG_SOMATIC_Indel_VCF -O z $My_UG_Indel_VCF $My_UG_CTR_Indel_VCF 
			time bcftools isec -n-1 -w1 -o $My_HC_SOMATIC_Indel_VCF -O z $My_HC_Indel_VCF $My_HC_CTR_Indel_VCF 

			#PASS only
			time bcftools view -fPASS -o $My_UG_SOMATIC_Indel_VCF -O z $My_UG_SOMATIC_Indel_VCF 
			time bcftools view -fPASS -o $My_HC_SOMATIC_Indel_VCF -O z $My_HC_SOMATIC_Indel_VCF 

			#Index it
			time tabix -f -p vcf $My_UG_SOMATIC_Indel_VCF
			time tabix -f -p vcf $My_HC_SOMATIC_Indel_VCF

			#intersection of HC and UG
			time bcftools isec -n=2 -w1 -o $My_SOMATIC_Indel_VCF -O z $My_HC_SOMATIC_Indel_VCF $My_UG_SOMATIC_Indel_VCF
			time tabix -f -p vcf $My_SOMATIC_Indel_VCF
		else
			echo -e "\e[031m$My_UG_Indel_VCF or $My_UG_CTR_Indel_VCF or $My_HC_Indel_VCF or $My_HC_CTR_Indel_VCF not found\e[0m\n"
			exit
		fi
	fi
fi

echo -e "\e[32m$Barcode done for per.barcode.gatk.somatic\e[0m" 
