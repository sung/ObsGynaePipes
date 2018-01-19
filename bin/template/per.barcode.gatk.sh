#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 30/Jan/2015
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

# SLX-7634.HALO1.000000000-A5NK2.bwa.bam
My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.bam 
if [ -s $My_bam ];then
	if [ ! -s ${My_bam%.bam}.sorted.bam ];then
		echo -e "samtools sort $My_bam ${My_bam%.bam}.sorted\n"
		time samtools sort $My_bam ${My_bam%.bam}.sorted
	fi

	if [ -s ${My_bam%.bam}.sorted.bam ];then
		# remove unsorted bam file
		echo -e "rm -f $My_bam\n"
		time rm -f $My_bam
	else
		echo -e "\e[031m${My_bam%.bam}.sorted.bam is expected but not found\e[0m\n"
		exit
	fi
else
	#SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam
	if [ ! -s ${My_bam%.bam}.sorted.bam ];then
		#SLX-7634.HALO1.000000000-A5NK2.s_1.bwa.bam
		My_bam_array=`ls $PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.s_*.bwa.bam`
		Dummy_bam_files=$(printf " %s" "${My_bam_array[@]}")
		#My_bam_array is not null
		if [ -n $My_bam_array ]; then
			# two or more bam files
			if [ ${#My_bam_array[@]} -gt 1 ]; then
				# merge bam files
				echo -e "samtools merge $My_bam $Dummy_bam_files"
				time samtools merge $My_bam $Dummy_bam_files 
			else
				echo -e "cp $Dummy_bam_files $My_bam\n"
				time cp $Dummy_bam_files $My_bam
			fi

			if [ -s $My_bam ];then
				echo -e "samtools sort $My_bam ${My_bam%.bam}.sorted\n"
				time samtools sort $My_bam ${My_bam%.bam}.sorted

				if [ -s ${My_bam%.bam}.sorted.bam ];then
					# remove unsorted bam file
					echo -e "rm -f $My_bam\n"
					time rm -f $My_bam
				fi
			else
				echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
				exit
			fi
		fi # end of My_bam_array is not null
	else
		echo -e "\e[032m${My_bam%.bam}.sorted.bam found\e[0m\n"
	fi
fi # end of if [ -s $My_bam ];then

#SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam
My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.bam 
if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
	exit
fi
#SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam.bai
if [ ! -s $My_bam".bai" ];then
	echo -e "samtools index $My_bam\n"
	time samtools index $My_bam
fi

#########################
## 1. Indel Realignment
#########################
My_interval=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.gatk.bam.intervals
if [ $RUN_GATK_INTERVAL -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/GATK
	mkdir_unless $PROJECT_DIR/GATK/$Barcode
	#Creating Intervals
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator\n"
	time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-I $My_bam \
		-R $GENOME \
		-nt $NT \
		-L $TARGET \
		--filter_mismatching_base_and_quals \
		-known $GATK_MILLS \
		-known $GATK_10KIndel \
		--out $My_interval
	echo -e "\e[032mDone Creating Interval\e[0m\n"
fi

if [ $RUN_GATK_REALIGN -eq 1 ]; then
	if [ -s $My_interval ]; then
		echo -e "time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T IndelRealigner\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-I $My_bam \
			-R $GENOME \
			-targetIntervals $My_interval \
			-known $GATK_MILLS \
			-known $GATK_10KIndel \
			-model USE_READS \
			--filter_mismatching_base_and_quals \
			--filter_bases_not_stored \
			--out ${My_bam%.bam}.realigned.bam
		echo -e "\e[032mDone Indel Realigning\e[0m\n"
	
		if [ -s ${My_bam%.bam}.realigned.bam ]; then
			if [ -s ${My_bam%.bam}.realigned.bai ]; then
				echo -e "mv ${My_bam%.bam}.realigned.bai ${My_bam%.bam}.realigned.bam.bai\n"
				mv ${My_bam%.bam}.realigned.bai ${My_bam%.bam}.realigned.bam.bai
			else
				echo -e "samtools index ${My_bam%.bam}.realigned.bam\n"
				time samtools index ${My_bam%.bam}.realigned.bam 
			fi
		else
			echo -e "\e[031m${My_bam%.bam}.realigned.bam not found\e[0m\n"
			exit
		fi
	else
		echo -e "\e[031m$My_interval not found\e[0m\n"
		exit
	fi
fi

My_recalTable=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.gatk.realigned.bam.recalibTable.grp
if [ $RUN_GATK_BASE_RECAL -eq 1 ]; then
	if [ -s ${My_bam%.bam}.realigned.bam ]; then
		###BaseRecalibrator
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar --T BaseRecalibrator\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			-I $My_bam \
			-L $TARGET \
			-R $GENOME \
			-nct $NT \
			-I ${My_bam%.bam}.realigned.bam \
			--knownSites $DBSNP \
			--knownSites $GATK_MILLS \
			--knownSites $GATK_10KIndel \
			-cov ReadGroupCovariate \
			-cov QualityScoreCovariate \
			-cov CycleCovariate \
			-cov ContextCovariate \
			--out $My_recalTable 
		echo -e "\e[032mDone: BaseRecalibrator\e[0m\n"
	else
		echo -e "\e[031m${My_bam%.bam}.realigned.bam not found\e[0m\n"
		exit
	fi
fi

if [ $RUN_GATK_PRNT_READ -eq 1 ];then
	if [ -s ${My_bam%.bam}.realigned.bam ] & [ -s $My_recalTable ]; then
		###PrintReads
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T PrintReads\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T PrintReads \
			-I ${My_bam%.bam}.realigned.bam \
			-R $GENOME \
			-nct $NT \
			-BQSR $My_recalTable \
			-S LENIENT \
			--out  ${My_bam%.bam}.realigned.recalibrated.bam
	else
		echo -e "\e[031m${My_bam%.bam}.realigned.bam not found\e[0m\n"
		exit

	fi

	if [ -s ${My_bam%.bam}.realigned.recalibrated.bam ]; then
		if [ -s ${My_bam%.bam}.realigned.recalibrated.bai ]; then
			mv ${My_bam%.bam}.realigned.recalibrated.bai ${My_bam%.bam}.realigned.recalibrated.bam.bai
		else
			echo -e "samtools index ${My_bam%.bam}.realigned.recalibrated.bam\n"
			time samtools index ${My_bam%.bam}.realigned.recalibrated.bam 
		fi

		#SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.bam
		#echo -e "rm $My_bam\n"
		#time rm $My_bam

		echo -e "rm ${My_bam%.bam}.realigned.bam\n"
		time rm ${My_bam%.bam}.realigned.bam
		time rm ${My_bam%.bam}.realigned.bam.bai
	else
		echo -e "\e[031m${My_bam%.bam}.realigned.recalibrated.bam not found\e[0m\n"
		exit
	fi
fi

#SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.realigned.recalibrated.bam
My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.bam 
if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
	exit
fi
################################################
# Remove reads of poor mapping quality (MAPQ) ##
# QC: by mapping quality (MAPQ) - remove poorly mapped reads based on MAPQ
####################################################
if [ ! -s ${My_bam%.bam}.OnTarget.q$MAPQ.bam ];then
	### Generate read on Target MAP qual > 8
	echo -e "samtools view -uq $MAPQ $My_bam | intersectBed -abam stdin -b $TARGET -u > ${My_bam%.bam}.OnTarget.q$MAPQ.bam\n"
	time samtools view -uq $MAPQ $My_bam | intersectBed -abam stdin -b $TARGET -u > ${My_bam%.bam}.OnTarget.q$MAPQ.bam
	echo -e "samtools index ${My_bam%.bam}.OnTarget.q$MAPQ.bam\n"
	time samtools index ${My_bam%.bam}.OnTarget.q$MAPQ.bam
fi

if [ ! -s ${My_bam%.bam}.OnTarget.q$MAPQ.bam.flagstat ]; then
	### Generating flagstat ## omit chimeric alignment reads -F 0x100
	echo -e "samtools view -u -F 0x100 ${My_bam%.bam}.OnTarget.q$MAPQ.bam | samtools flagstat - > ${My_bam%.bam}.OnTarget.q$MAPQ.bam.flagstat\n"
	time samtools view -u -F 0x100 ${My_bam%.bam}.OnTarget.q$MAPQ.bam | samtools flagstat - > ${My_bam%.bam}.OnTarget.q$MAPQ.bam.flagstat
fi

#SLX-7634.HALO1.000000000-A5NK2.bwa.sorted.realigned.recalibrated.OnTarget.q1.bam
My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.bam 
if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
	exit
fi

My_UG_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.snp.indel.vcf.gz
My_UG_SNP_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.snp.vcf.gz 
My_UG_Indel_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.indel.vcf.gz 
My_UG_metrics=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.UG.snp.indel.metrics
if [ $RUN_GATK_UG -eq 1 ];then
	###Calling UnifiedGenotyper  ###
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper\n"
	time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T UnifiedGenotyper \
		-R $GENOME \
		-nt $NT \
		-I $My_bam \
		-L $TARGET \
		-l INFO \
		-baq CALCULATE_AS_NECESSARY \
		-minIndelCnt $GATK_MIN_INDEL_CNT \
		-glm BOTH \
		--dbsnp $DBSNP \
		-A Coverage \
		-A AlleleBalance \
		-G Standard \
		--min_base_quality_score $GATK_MIN_BASE_QUAL \
		-stand_call_conf $GATK_CALL_CONF \
		-stand_emit_conf $GATK_EMIT_CONF \
		--downsample_to_coverage $GATK_DCOVG \
		--downsampling_type BY_SAMPLE \
		--computeSLOD  \
		--metrics_file $My_UG_metrics \
		-o $My_UG_VCF

	echo -e "\e[032mDone: GATK UG (SNP & Indels initial calling)\e[0m\n"

	if [ -s $My_UG_VCF ];then
		#### print Indel variant
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T SelectVariants --selectTypeToInclude INDEL\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T SelectVariants \
			-R $GENOME \
			-nt $NT \
			--variant $My_UG_VCF \
			--selectTypeToInclude INDEL \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-o $My_UG_Indel_VCF 
		echo -e "\e[032mDone : GATK UG Indels calling\e[0m\n"

		#### print SNP variant
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T SelectVariants --selectTypeToInclude SNP\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T SelectVariants \
			-R $GENOME \
			-nt $NT \
			--variant $My_UG_VCF \
			--selectTypeToInclude SNP \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-o $My_UG_SNP_VCF
		echo -e "\e[032mDone : GATK UG SNP calling\e[0m\n"
	else
		echo -e "\e[031m$My_UG_VCF is expected but not found\e[0m\n"
		exit
	fi

	if [ -s $My_UG_Indel_VCF ];then
		###Indel filter
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T VariantFiltration\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $GENOME \
			--variant $My_UG_Indel_VCF \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-baq CALCULATE_AS_NECESSARY \
			--filterExpression "QD < $GATK_QD" \
			--filterName QDFilter \
			-o ${My_UG_Indel_VCF%.vcf.gz}.filtered.vcf.gz
		echo -e "\e[032mDone : GATK UG Indels filtering\e[0m\n"
		time tabix -f -p vcf ${My_UG_Indel_VCF%.vcf.gz}.filtered.vcf.gz
	else
		echo -e "\e[031m$My_UG_Indel_VCF is expected but not found\e[0m\n"
		exit
	fi

	if [ -s $My_UG_SNP_VCF ];then
		### SNPs Filter
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T VariantFiltration\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $GENOME \
			--variant $My_UG_SNP_VCF \
			--mask $My_UG_Indel_VCF \
			--maskName InDel \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-baq CALCULATE_AS_NECESSARY \
			--filterExpression "QD < $GATK_QD" \
			--filterExpression "MQ < $GATK_MQ" \
			--filterName QDFilter \
			--filterName MQFilter \
			-o ${My_UG_SNP_VCF%.vcf.gz}.filtered.vcf.gz
		echo -e "\e[032mDone : GATK UG SNP filtering\e[0m\n"
		time tabix -f -p vcf ${My_UG_SNP_VCF%.vcf.gz}.filtered.vcf.gz
	else
		echo -e "\e[031m$My_UG_SNP_VCF is expected but not found\e[0m\n"
		exit
	fi

	echo -e "\e[032mGATK UG done\e[0m\n"
fi

#My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.bam 
My_HC_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.snp.indel.vcf.gz 
My_HC_Indel_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.indel.vcf.gz
My_HC_SNP_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.HC.snp.vcf.gz
if [ $RUN_GATK_HC -eq 1 ];then
	###Calling Initial HaplotypeCaller###
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T HaplotypeCaller\n"
	time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R $GENOME \
		-nct $NT \
		-I $My_bam \
		-L $TARGET \
		-l INFO \
		--dbsnp $DBSNP \
		-A Coverage \
		-A AlleleBalance \
		-G Standard \
		--min_base_quality_score $GATK_MIN_BASE_QUAL \
		-stand_call_conf $GATK_CALL_CONF \
		-stand_emit_conf $GATK_EMIT_CONF \
		--downsample_to_coverage $GATK_DCOVG \
		--downsampling_type BY_SAMPLE \
		-o $My_HC_VCF

	echo -e "\e[032mDone: GATK HC (SNP & Indels initial calling)\e[0m\n"
	
	if [ -s $My_HC_VCF ];then
		### Variant Annotate
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T VariantAnnotator\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T VariantAnnotator \
			-R $GENOME \
			-nt $NT \
			-I $My_bam \
			--variant $My_HC_VCF \
			-L $TARGET \
			-A Coverage \
			-A AlleleBalance \
			-A HaplotypeScore \
			-A InbreedingCoeff \
			-A HomopolymerRun \
			-A HardyWeinberg \
			-A GCContent \
			--dbsnp $DBSNP \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-o ${My_HC_VCF%.vcf.gz}.varanno.vcf.gz
	else
		echo -e "\e[031m$My_HC_VCF is expected but not found\e[0m\n"
		exit
	fi

	if [ -s ${My_HC_VCF%.vcf.gz}.varanno.vcf.gz ];then
		#### print Indel variant
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T SelectVariants --selectTypeToInclude INDEL\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T SelectVariants \
			-R $GENOME \
			-nt $NT \
			--variant ${My_HC_VCF%.vcf.gz}.varanno.vcf.gz \
			--selectTypeToInclude INDEL \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-o $My_HC_Indel_VCF 
		echo -e "\e[032mDone : GATK HC Indels calling\e[0m\n"

		#### print SNP variant
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T SelectVariants --selectTypeToInclude SNP\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T SelectVariants \
			-R $GENOME \
			-nt $NT \
			--variant ${My_HC_VCF%.vcf.gz}.varanno.vcf.gz \
			--selectTypeToInclude SNP \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-o $My_HC_SNP_VCF
		echo -e "\e[032mDone : GATK HC SNP calling\e[0m\n"

	else
		echo -e "\e[031m${My_HC_VCF%.vcf.gz}.varanno.vcf.gz is expected but not found\e[0m\n"
		exit
	fi

	if [ -s $My_HC_Indel_VCF ];then
		###Indel filter
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T VariantFiltration\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $GENOME \
			--variant $My_HC_Indel_VCF \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-baq CALCULATE_AS_NECESSARY \
			--filterExpression "QD < $GATK_QD" \
			--filterName QDFilter \
			-o ${My_HC_Indel_VCF%.vcf.gz}.filtered.vcf.gz
		echo -e "\e[032mDone : GATK HC Indels filtering\e[0m\n"
		time tabix -f -p vcf ${My_HC_Indel_VCF%.vcf.gz}.filtered.vcf.gz
	else
		echo -e "\e[031m$My_HC_Indel_VCF is expected but not found\e[0m\n"
		exit
	fi

	if [ -s $My_HC_SNP_VCF ];then
		### SNPs Filter
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T VariantFiltration\n"
		time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $GENOME \
			--variant $My_HC_SNP_VCF \
			--mask $My_HC_Indel_VCF \
			--maskName InDel \
			--downsample_to_coverage $GATK_DCOVG \
			--downsampling_type BY_SAMPLE \
			-baq CALCULATE_AS_NECESSARY \
			--filterExpression "QD < $GATK_QD" \
			--filterExpression "MQ < $GATK_MQ" \
			--filterName QDFilter \
			--filterName MQFilter \
			-o ${My_HC_SNP_VCF%.vcf.gz}.filtered.vcf.gz
		echo -e "\e[032mDone : GATK HC SNP filtering\e[0m\n"
		time tabix -f -p vcf ${My_HC_SNP_VCF%.vcf.gz}.filtered.vcf.gz
	else
		echo -e "\e[031m$My_HC_SNP_VCF is expected but not found\e[0m\n"
		exit
	fi

	echo -e "\e[032mGATK HC done\e[0m\n"
fi

echo -e "\e[32m$Barcode done for per.barcode.gatk\e[0m" 
