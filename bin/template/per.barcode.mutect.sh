#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 19/Mar/2015
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

#######
# PON #
#######
#for i in ${CONTROL[*]}; do echo $i; bgzip -f $PROJECT_DIR/GATK/$SLX.$i.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UnifiedGenotyper.snp.filtered.vcf; tabix -f -p vcf $PROJECT_DIR/GATK/$SLX.$i.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UnifiedGenotyper.snp.filtered.vcf.gz; done
#Merged_vcf=$(printf " %s" `for i in ${CONTROL[*]}; do echo $PROJECT_DIR/GATK/$SLX.$i.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UnifiedGenotyper.snp.filtered.vcf.gz; done`)
#Merged_control=$(printf ".%s" "${CONTROL[@]}") # .HAL40.HAL41.HAL42.HAL43
#Merged_control=${Merged_control#.} # HAL40.HAL41.HAL42.HAL43
#time bcftools merge  $Merged_control > $PROJECT_DIR/GATK/$SLX.$Merged_control.UG.vcf

if [ ! $My_control ]; then
	echo -e "\e[31mNo control Barcode found: $Barcode must be control\e[0m" 
	exit
else
	################
	## 1. Mutect  ##
	################
	mkdir_unless $PROJECT_DIR/MUTECT
	mkdir_unless $PROJECT_DIR/MUTECT/$Barcode

	Merged_control=$(printf ".%s" "${CONTROL[@]}") # .HAL40.HAL41.HAL42.HAL43
	Merged_control=${Merged_control#.} # HAL40.HAL41.HAL42.HAL43

	# SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.realigned.recalibrated.bam
	My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.bam
	My_control_bam=$PROJECT_DIR/BWA/$My_control/$SLX.$My_control.$Cell.bwa.sorted.realigned.recalibrated.bam
	My_PON_vcf=$PROJECT_DIR/GATK/$SLX.$Merged_control.UG.vcf
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi
	if [ ! -s $My_control_bam ];then
		echo -e "\e[031m$My_control_bam is expected but not found\e[0m\n"
		exit
	fi
	if [ ! -s $My_PON_vcf ];then
		echo -e "\e[031m$My_PON_vcf is expected but not found\e[0m\n"
		exit
	fi
	# SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.realigned.recalibrated.bam.bai
	if [ ! -s $My_bam".bai" ];then
		echo -e "samtools index $My_bam\n"
		time samtools index $My_bam
	fi
	if [ ! -s $My_control_bam".bai" ];then
		echo -e "samtools index $My_control_bam\n"
		time samtools index $My_control_bam
	fi

	Mutect_pair_vcf=$PROJECT_DIR/MUTECT/$Barcode/$SLX.$Barcode.mutect.pair.vcf.gz
	Mutect_pair_out=$PROJECT_DIR/MUTECT/$Barcode/$SLX.$Barcode.mutect.pair.out.txt
	if [ $RUN_MUTECT_PAIR -eq 1 ]; then
		echo -e "java -Xmx$GATK_MEM -jar $MUTECT $My_bam...\n"
		time java -Xmx$GATK_MEM -jar $MUTECT \
			--analysis_type MuTect \
			--reference_sequence $GENOME \
			--intervals $TARGET \
			--dbsnp $DBSNP \
			--cosmic $COSMIC \
			--input_file:normal $My_control_bam \
			--input_file:tumor $My_bam \
			--out $Mutect_pair_out \
			--vcf $Mutect_pair_vcf
		if [ -s $Mutect_pair_vcf ]; then
			time tabix -f -p vcf $Mutect_pair_vcf
		else
			echo -e "\e[031m$Mutect_pair_vcf expected but not found\e[0m\n"
			exit
		fi
		# filter only PASS
		#time bcftools view -fPASS $Mutect_pair_vcf > ${Mutect_pair_vcf%.vcf.gz}.passed.vcf
	fi

	Mutect_pon_vcf=$PROJECT_DIR/MUTECT/$Barcode/$SLX.$Barcode.mutect.pon.vcf.gz
	Mutect_pon_out=$PROJECT_DIR/MUTECT/$Barcode/$SLX.$Barcode.mutect.pon.out.txt
	if [ $RUN_MUTECT_PON -eq 1 ]; then
		echo -e "java -Xmx$GATK_MEM -jar $MUTECT $My_bam...\n"
		time java -Xmx$GATK_MEM -jar $MUTECT \
			--analysis_type MuTect \
			--reference_sequence $GENOME \
			--intervals $TARGET \
			--dbsnp $DBSNP \
			--cosmic $COSMIC \
			--input_file:tumor $My_bam \
			--normal_panel $My_PON_vcf \
			--out $Mutect_pon_out \
			--vcf $Mutect_pon_vcf
		if [ -s $Mutect_pon_vcf ]; then
			time tabix -f -p vcf $Mutect_pon_vcf
		else
			echo -e "\e[031m$Mutect_pon_vcf expected but not found\e[0m\n"
			exit
		fi
		# filter only PASS
		#time bcftools view -fPASS $Mutect_pon_vcf > ${Mutect_pon_vcf%.vcf.gz}.passed.vcf
	fi

	Mutect_merged_vcf=$PROJECT_DIR/MUTECT/$Barcode/$SLX.$Barcode.mutect.merged.vcf.gz
	if [ $RUN_MUTECT_MERGE -eq 1 ]; then
		if [ -s $Mutect_pair_vcf ] && [ -s $Mutect_pon_vcf ]; then
			echo -e "bcftools merge --force-samples -fPASS --info-rules - --output-type z --output $Mutect_merged_vcf $Mutect_pair_vcf $Mutect_pon_vcf\n"
			time bcftools merge --force-samples -fPASS --info-rules - --output-type z --output $Mutect_merged_vcf $Mutect_pair_vcf $Mutect_pon_vcf
			if [ -s $Mutect_merged_vcf ]; then
				time tabix -f -p vcf $Mutect_merged_vcf
			else
				echo -e "\e[031m$Mutect_merged_vcf expected but not found\e[0m\n"
				exit
			fi
		else
			echo -e "\e[031m$Mutect_pair_vcf or $Mutect_pon_vcf not found\e[0m\n"
			exit
		fi
	fi
fi

echo -e "\e[32m$Barcode done for per.barcode.mutect\e[0m" 
