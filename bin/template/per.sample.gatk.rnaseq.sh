#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 17/Aug/2015
# Last modified: 2/Jun/2016
# Copied mostly from ~/Pipelines/bin/template/per.barcode.gatk.sh (for amplicon-based targetted resequencing)
# https://www.broadinstitute.org/gatk/guide/article?id=3891
# https://www.broadinstitute.org/gatk/guide/article?id=1247

source $HOME/RNA-Seq/config/rna_seq.config # export envrionment variables
source $HOME/lib/sung.sh
source $HOME/config/sung.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir_unless $PROJECT_DIR
mkdir_unless $PROJECT_DIR/GATK
mkdir_unless $PROJECT_DIR/scratch
mkdir_unless $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Cell.$Barcode\e[0m\n"

#My_bam=$PROJECT_DIR/STAR/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.bam # STAR 2-pass
#My_bam=$PROJECT_DIR/TopHat/$Barcode/accepted_hits.bam                # TopHat 1-pass
My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.bam          # TopHat 2-pass
########################################
## 0. Add RG tag for TopHat bam files ##
########################################
if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
	exit
fi

if [ $RUN_PICARD_RG -eq 1 ]; then
	echo -e "java -Xmx$RG_MEM -jar $PICARD/AddOrReplaceReadGroups.jar $My_bam\n"
	time java -Xmx$RG_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $PICARD/AddOrReplaceReadGroups.jar \
		INPUT=$My_bam \
		OUTPUT=${My_bam%.bam}.RG.bam \
		RGID=$SLX \
		RGLB=$SLX \
		RGPL=illumina \
		RGPU=run \
		RGSM=$Barcode \
		RGCN=$RG_CN \
		VALIDATION_STRINGENCY=$RG_STRINGENCY \
		SORT_ORDER=$RG_ORDER
fi

#input: merged_accepted_hits.RG.bam
#output: merged_accepted_hits.RG.markDup.bam
######################
## 1. MarkDuplicate ##
######################
if [ $RUN_PICARD_MD -eq 1 ]; then
	My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.bam 
	if [ ! -s $My_bam ]; then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi

	echo -e "java -Xmx$RG_MEM -jar $PICARD/MarkDuplicates.jar $My_bam\n"
	time java -Xmx$RG_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $PICARD/MarkDuplicates.jar \
		INPUT=$My_bam \
		OUTPUT=${My_bam%.bam}.markDup.bam \
		CREATE_INDEX=true \
		METRICS_FILE=${My_bam%.bam}.markDup.metrics.txt \
		VALIDATION_STRINGENCY=$RG_STRINGENCY
	echo -e "\e[32m$SLX.$Barcode done for Picard MarkDuplicate\e[0m" 
fi

if [ -s ${My_bam%.bam}.markDup.bai ]; then
	mv ${My_bam%.bam}.markDup.bai ${My_bam%.bam}.markDup.bam.bai
fi

######################################################
# 2. Split'N'Trim and reassign mapping qualities     #
# Splits reads that contain Ns in their CIGAR string #
######################################################
#input: merged_accepted_hits.RG.markDup.bam
#output: merged_accepted_hits.RG.markDup.split.bam
#https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_SplitNCigarReads.php
if [ $RUN_GATK_SPLIT -eq 1 ]; then
	My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.bam 
	if [ ! -s $My_bam ]; then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T SplitNCigarReads -I $My_bam\n"
	time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T SplitNCigarReads \
		-R $STAR_GENOME \
		-I $My_bam \
		-o ${My_bam%.bam}.split.bam \
		-rf ReassignOneMappingQuality \
		-RMQF 255 \
		-RMQT 60 \
		-U ALLOW_N_CIGAR_READS
	echo -e "\e[32m$SLX.$Barcode done for GATK SplitNCigarReads\e[0m" 
fi

if [ -s ${My_bam%.bam}.split.bai ]; then
	mv ${My_bam%.bam}.split.bai ${My_bam%.bam}.split.bam.bai
fi

# input: merged_accepted_hits.RG.markDup.split.bam
# output: merged_accepted_hits.RG.markDup.split.realigned.bam
##########################
## 3. Indel Realignment ##
##########################
if [ $RUN_GATK_REALIGN -eq 1 ]; then
	My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.split.bam 
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi

	mkdir_unless $PROJECT_DIR/GATK/$Barcode
	My_interval=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.gatk.bam.intervals # will be made
	# Creating Intervals
	if [ ! -s $My_interval ];then
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $My_bam\n"
		time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T RealignerTargetCreator \
			-nt $NT \
			-I $My_bam \
			-R $STAR_GENOME \
			--filter_mismatching_base_and_quals \
			-known $GATK_MILLS \
			--out $My_interval
			#-known $GATK_10KIndel \
			#-L $TARGET \
		echo -e "\e[032mDone Creating Interval\e[0m\n"
	fi
	if [ -s $My_interval ]; then
		echo -e "time java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I $My_bam\n"
		time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-I $My_bam \
			-R $STAR_GENOME \
			-targetIntervals $My_interval \
			-known $GATK_MILLS \
			-model USE_READS \
			--filter_mismatching_base_and_quals \
			--filter_bases_not_stored \
			--out ${My_bam%.bam}.realigned.bam
			#-known $GATK_10KIndel \
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

# input: merged_accepted_hits.RG.markDup.split.realigned.bam
# output: merged_accepted_hits.RG.markDup.split.realigned.recalibrated.bam
###################################
## 3. Base-quality Recalibration ##
###################################
if [ $RUN_GATK_BASE_RECAL -eq 1 ]; then
	My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.split.realigned.bam 
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi

	mkdir_unless $PROJECT_DIR/GATK/$Barcode
	My_recalTable=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.$Cell.gatk.realigned.bam.recalibTable.grp # will be made
	###BaseRecalibrator
	if [ -s $My_bam ]; then
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar --T BaseRecalibrator -I $My_bam \n"
		time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			-nct $NT \
			-I $My_bam \
			-R $STAR_GENOME \
			--knownSites $DBSNP \
			--knownSites $GATK_MILLS \
			-cov ReadGroupCovariate \
			-cov QualityScoreCovariate \
			-cov CycleCovariate \
			-cov ContextCovariate \
			--out $My_recalTable 
			#--knownSites $GATK_10KIndel \
			#-L $TARGET \
		echo -e "\e[032mDone: BaseRecalibrator\e[0m\n"
	else
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi
	###PrintReads
	if [ -s $My_bam ] & [ -s $My_recalTable ]; then
		echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T PrintReads -I $My_bam\n"
		time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-T PrintReads \
			-nct $NT \
			-I $My_bam \
			-R $STAR_GENOME \
			-BQSR $My_recalTable \
			-S LENIENT \
			--out  ${My_bam%.bam}.recalibrated.bam
	else
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit

	fi
	## Index bam 
	if [ -s ${My_bam%.bam}.recalibrated.bam ]; then
		if [ -s ${My_bam%.bam}.recalibrated.bai ]; then
			mv ${My_bam%.bam}.recalibrated.bai ${My_bam%.bam}.recalibrated.bam.bai
		else
			echo -e "samtools index ${My_bam%.bam}.recalibrated.bam\n"
			time samtools index ${My_bam%.bam}.recalibrated.bam 
		fi
	else
		echo -e "\e[031m${My_bam%.bam}.recalibrated.bam not found\e[0m\n"
		exit
	fi
fi

# input: merged_accepted_hits.RG.markDup.split.realigned.recalibrated.bam
# output: $SLX.$Barcode.tophat.HC.snp.indel.vcf.gz
My_HC_VCF=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.tophat.HC.snp.indel.vcf.gz
#####################
## Variant Calling ##
#####################
#https://www.broadinstitute.org/gatk/guide/article?id=3891
if [ $RUN_GATK_HC -eq 1 ];then
	My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.split.realigned.recalibrated.bam
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi

	mkdir_unless $PROJECT_DIR/GATK/$Barcode
	###Calling Initial HaplotypeCaller###
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T HaplotypeCaller -I $My_bam\n"
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode $GATK_HOME/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-nct $NT \
		-R $STAR_GENOME \
		-I $My_bam \
		--dbsnp $DBSNP \
		-dontUseSoftClippedBases \
		-stand_call_conf $GATK_CALL_CONF \
		-stand_emit_conf $GATK_EMIT_CONF \
		-o $My_HC_VCF
	echo -e "\e[032mDone: GATK HC (SNP & Indels initial calling)\e[0m\n"
	
	if [ ! -s $My_HC_VCF ];then
		echo -e "\e[031m$My_HC_VCF is expected but not found\e[0m\n"
		exit
	fi

	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T VariantFiltration --variant $My_HC_VCF\n"
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode $GATK_HOME/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R $STAR_GENOME \
		--variant $My_HC_VCF \
		-window 35 -cluster 3 \
		-filterName FS -filter "FS > $GATK_FS" \
		-filterName QD -filter "QD < $GATK_QD" \
		-o ${My_HC_VCF%.vcf.gz}.filtered.vcf.gz
	echo -e "\e[032mGATK HC done\e[0m\n"

	if [ ! -s ${My_HC_VCF%.vcf.gz}.filtered.vcf.gz ];then
		echo -e "\e[031m${My_HC_VCF%.vcf.gz}.filtered.vcf.gz is expected but not found\e[0m\n"
		exit
	fi
	#Index
	echo -e "tabix -f -p vcf ${My_HC_VCF%.vcf.gz}.filtered.vcf.gz\n"
	time tabix -f -p vcf ${My_HC_VCF%.vcf.gz}.filtered.vcf.gz

	#https://www.broadinstitute.org/gatk/guide/article?id=1255
	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T SelectVariants --variant ${My_HC_VCF%.vcf.gz}.filtered.vcf.gz\n"
	time java -Xmx$GATK_MEM -jar -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode $GATK_HOME/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-nt $NT \
		-R $STAR_GENOME \
		--variant ${My_HC_VCF%.vcf.gz}.filtered.vcf.gz \
		--selectTypeToExclude INDEL \
		--restrictAllelesTo BIALLELIC \
		-select "vc.getGenotype(\"$Barcode\").isHet()" \
		-o ${My_HC_VCF%.snp.indel.vcf.gz}.filtered.het.snp.vcf.gz
		#-select "GT == 0/1" (not valid)
	echo -e "\e[032mGATK HC done\e[0m\n"

	if [ ! -s ${My_HC_VCF%.snp.indel.vcf.gz}.filtered.het.snp.vcf.gz ];then
		echo -e "\e[031m${My_HC_VCF%.snp.indel.vcf.gz}.filtered.het.snp.vcf.gz is expected but not found\e[0m\n"
		exit
	fi
	#Index
	echo -e "tabix -f -p vcf ${My_HC_VCF%.snp.indel.vcf.gz}.filtered.het.snp.vcf.gz\n"
	time tabix -f -p vcf ${My_HC_VCF%.snp.indel.vcf.gz}.filtered.het.snp.vcf.gz

	echo -e "\e[032mDone : GATK HC Indels filtering\e[0m\n"
fi

##################
## ASE Counting ##
##################
#https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php
My_ASE=$PROJECT_DIR/GATK/$Barcode/$SLX.$Barcode.tophat.HC.ASE.txt 
if [ $RUN_GATK_ASE -eq 1 ];then
	My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.split.realigned.recalibrated.bam
	if [ -s $My_bam ];then
		# Remove unnecessary bam files
		echo -e "rm $PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.bam\n"
		time rm $PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.bam*
		time rm $PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.bam*
		time rm $PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.split.bam*
		time rm $PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.RG.markDup.split.realigned.bam*
	else
		echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
		exit
	fi

	echo -e "java -Xmx$GATK_MEM -jar $GATK_HOME/GenomeAnalysisTK.jar -T ASEReadCounter -I $My_bam\n"
	time java -Xmx$GATK_MEM -Djava.io.tmpdir=$PROJECT_DIR/scratch/$Barcode -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-T ASEReadCounter \
		-R $STAR_GENOME \
		-I $My_bam \
		-sites ${My_HC_VCF%.snp.indel.vcf.gz}.filtered.het.snp.vcf.gz \
		-U ALLOW_N_CIGAR_READS \
		-minDepth 10 \
		--minMappingQuality 10 \
		--minBaseQuality 2 \
		-o $My_ASE
		#-drf DuplicateRead \
fi

echo -e "\e[32m$Barcode done for per.barcode.gatk\e[0m" 
