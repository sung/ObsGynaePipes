#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 26/Mar/2015
# Last modified: 26/Mar/2015 

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
	#####################################
	## Merge GATK Somatic SNP & MUTECT ##
	#####################################
	mkdir_unless $PROJECT_DIR/SomaticVariant
	mkdir_unless $PROJECT_DIR/SomaticVariant/$Barcode

	RENAME=$PROJECT_DIR/SomaticVariant/$Barcode/$Barcode.rename.txt
	time echo "$Barcode $Barcode.gatk" > $RENAME 

	# Input VCF
	MUTECT_VCF=$PROJECT_DIR/MUTECT/$Barcode/$SLX.$Barcode.mutect.merged.vcf.gz
	GATK_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.somatic.snp.vcf.gz 
	New_GATK_VCF=$PROJECT_DIR/SomaticVariant/$Barcode/$SLX.$Barcode.gatk.somatic.snp.vcf.gz 

	# Change the header of GATK VCF
	time bcftools reheader -s $RENAME -o $New_GATK_VCF $GATK_VCF

	# Index the new GATK
	time tabix -f -p vcf $New_GATK_VCF

	# Output
	Common_VCF=$PROJECT_DIR/SomaticVariant/$Barcode/$SLX.$Barcode.mutect.gatk.common.vcf.gz
	Unified_VCF=$PROJECT_DIR/SomaticVariant/$Barcode/$SLX.$Barcode.mutect.gatk.merged.vcf.gz

	if [ -s $MUTECT_VCF ] && [ -s $New_GATK_VCF ]; then
		# common varaints between MUTECT & Somatic GATK
		time bcftools isec -n=2 -w1 -o $Common_VCF -O z $MUTECT_VCF $New_GATK_VCF
		time tabix -f -p vcf $Common_VCF
		# unified variants between MUTECT & Somatic GATK
		time bcftools merge --info-rules - --output-type z --output $Unified_VCF $MUTECT_VCF $New_GATK_VCF 
		time tabix -f -p vcf $Unified_VCF
		# remove renamed GATK VCF
		time rm $New_GATK_VCF*
	else
		echo -e "\e[031m$MUTECT_VCF or $New_GATK_VCF expected but not found\e[0m\n"
		exit
	fi
fi

echo -e "\e[32m$Barcode done for per.barcode.somatic.snp\e[0m" 
