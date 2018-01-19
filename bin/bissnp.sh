#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Reference: http://epigenome.usc.edu/publicationdata/bissnp2011/BisSNP-UserGuide-latest.pdf
# Called by bin/template/bissnp.sh
# First created: 20/May/2014
# Last modified: 15/Apr/2015

#source ~/Pipelines/bin/bissnp.SLX-10409.v2.patch.sh
if [ $IS_PE -eq 1 ]; then
	if [ $IS_FLD -eq 1 ]; then
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.sorted.bam # Fludigm
	else
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.sorted.bam # WGBS/SureSelect
	fi
else
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.sorted.bam 
fi

############################################################
# 1. Add RG (Read Group) to bam header
# 'CREATE_INDEX=true' will make unknown error - disabled
# input: my.q1.deduplicated.sorted.bam
# output: my.q1.deduplicated.sorted.RG.bam
# output: my.q1.deduplicated.sorted.RG.bam.bai
############################################################
if [ $RUN_PICARD_RG -eq 1 ]; then
	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	else
		echo java -Xmx$RG_MEM -jar $PICARD/AddOrReplaceReadGroups.jar $My_bam
		time java -Xmx$RG_MEM -jar $PICARD/AddOrReplaceReadGroups.jar \
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

	if [ -s ${My_bam%.bam}.RG.bam ]; then
		echo -e "samtools index ${My_bam%.bam}.RG.bam "
		time samtools index ${My_bam%.bam}.RG.bam 
	else
		echo -e "\e[031m${My_bam%.bam}.RG.bam not found\e[0m\n"
		exit
	fi
	echo -e "\e[32m$SLX.$Barcode done for Picard RG\e[0m" 
fi

#########################
## 2. Indel Realignment
#  If your mapping tools allows gapped alignment, like Bismark with Bowtie2 or Novoaligner, you could omit this step.
#########################

#######################
## 3. Mark duplicates
# so skip this step if done by 'deduplicate_bismark' previously
#######################

###################################
## 4. Base Quality Recalibration
###################################
#4-1. Count Covariant
# Currently, Bis-SNP only allows recalibration on 3 covariants: 
# 1) ReadGroupCovariate, 
# 2) QualityScoreCovariate and 
# 3) CycleCovariate
if [ $RUN_BS_BASE_RECAL -eq 1 ]; then
	if [ -s ${My_bam%.bam}.RG.bam ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T BisulfiteCountCovariates ${My_bam%.bam}.RG.bam"
		time java -Xmx$BIS_MEM -jar $BIS_SNP \
			-T BisulfiteCountCovariates \
			-nt $NT \
			-R $REF_SEQ \
			-I ${My_bam%.bam}.RG.bam \
			-knownSites $DBSNP \
			-cov ReadGroupCovariate \
			-cov QualityScoreCovariate \
			-cov CycleCovariate \
			-recalFile $PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.recalFile_before.csv
			#-L $TARGET \
	else
		echo -e "\e[031m${My_bam%.bam}.RG.bam not found\e[0m\n"
		exit
	fi
fi

if [ $RUN_BS_BASE_RECAL_BAM -eq 1 ]; then
	#4-2. Write recalibrated base quality score into BAM file
	# -maxQ: the maximum base quality score in the original BAM file 
	if [ -s $PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.recalFile_before.csv ]; then
		# BisulfiteTableRecalibration does not support "-nt"
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T BisulfiteTableRecalibration -I ${My_bam%.bam}.RG.bam"
		time java -Xmx$BIS_MEM -jar $BIS_SNP \
			-T BisulfiteTableRecalibration \
			-R $REF_SEQ \
			-I ${My_bam%.bam}.RG.bam \
			-recalFile $PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.recalFile_before.csv \
			-maxQ $BIS_MAXQ \
			-o ${My_bam%.bam}.RG.recal.bam 
			#-L $TARGET \
	else
		echo -e "\e[031m $PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.recalFile_before.csv not found\e[0m\n"
		exit
	fi
	echo -e "\e[32m$SLX.$Barcode done for BisSNP base quality score recal\e[0m" 
fi #end of $RUN_BS_BASE_RECAL

#4-3. Re-Count Covariant
# [optional] This step is to validate that if the recalibration step is correct
if [ $RUN_BS_OPT -eq 1 ]; then
	if [ -s ${My_bam%.bam}.RG.recal.bam ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T BisulfiteCountCovariates -I ${My_bam%.bam}.RG.recal.bam" 
		time java -Xmx$BIS_MEM -jar $BIS_SNP \
			-T BisulfiteCountCovariates \
			-nt $NT \
			-R $REF_SEQ \
			-I ${My_bam%.bam}.RG.recal.bam \
			-knownSites $DBSNP \
			-cov ReadGroupCovariate \
			-cov QualityScoreCovariate \
			-cov CycleCovariate \
			-recalFile $PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.recalFile_after.csv
	fi
	#4-4. Generate recalibration plot
	# [optional] 
	echo -e " java -Xmx$BIS_MEM -jar $BIS_SNP_COVARI -recalFile $Barcode/$SLX.$Barcode.recalFile_before.csv"
	time java -Xmx$BIS_MEM -jar $BIS_SNP_COVARI \
		-recalFile $PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.recalFile_before.csv \
		-outputDir $PROJECT_DIR/BisSNP/$Barcode \
		-ignoreQ 5 \
		--max_quality_score $BIS_MAXQ &

	echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP_COVARI -recalFile $Barcode/$SLX.$Barcode.recalFile_after.csv"
	time java -Xmx$BIS_MEM -jar $BIS_SNP_COVARI \
		-recalFile $PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.recalFile_after.csv \
		-outputDir $PROJECT_DIR/BisSNP/$Barcode \
		-ignoreQ 5 \
		--max_quality_score $BIS_MAXQ &
	wait
	echo -e "\e[32m$SLX.$Barcode done for BisSNP re-count covariant\e[0m" 
fi

# moved to bin/bissnp.per.chr.sh
##############################################
## 5. Bis-SNP genotying (BisulfiteGenotyper)
## output: snp.raw.vcf and cpg.raw.vcf
## The chromosome order of dbSNP VCF file should be the same as your reference genome file and input BAM file's header
## -C (or --cytosine contexts acquired)
## Default: ‘-C CG,1 -C CH,1’
## -out modes (or --output modes)
## DEFAULT FOR TCGA: emit all of CpG sites above emit threshold into vcf1 file, and all of SNP sites above
## emit threshold into vcf2 file
## EMIT VARIANT AND CYTOSINES: emit all of Cytosine sites above emit threshold into vcf1 file, and all of SNP
## sites above emit threshold into vcf2 file.
##############################################
VFN1=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.cpg.raw.vcf
VFN2=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.snp.raw.vcf 
if [ $RUN_BS_GENO -eq 1 ]; then
	# e.g. SLX-8074.SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.sorted.RG.recal.bam
	if [ -s ${My_bam%.bam}.RG.recal.bam ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T BisulfiteGenotyper -I ${My_bam%.bam}.RG.recal.bam"
		if [ $IS_FLD -eq 1 ]; then
			time java -Xmx$BIS_MEM -jar $BIS_SNP \
				-T BisulfiteGenotyper \
				-R $REF_SEQ \
				-I ${My_bam%.bam}.RG.recal.bam \
				-D $DBSNP \
				-nt $NT \
				-L $TARGET \
				-dcov $BIS_DCOV \
				-toCoverage $BIS_MAXC \
				-stand_call_conf $BIS_CALL_CONF \
				-stand_emit_conf $BIS_EMIT_CONF \
				-mmq $BIS_MMQ \
				-mbq $BIS_MBQ \
				-vfn1 $VFN1 \
				-vfn2 $VFN2
		else
			time java -Xmx$BIS_MEM -jar $BIS_SNP \
				-T BisulfiteGenotyper \
				-R $REF_SEQ \
				-I ${My_bam%.bam}.RG.recal.bam \
				-D $DBSNP \
				-nt $NT \
				-toCoverage $BIS_MAXC \
				-stand_call_conf $BIS_CALL_CONF \
				-stand_emit_conf $BIS_EMIT_CONF \
				-mmq $BIS_MMQ \
				-mbq $BIS_MBQ \
				-vfn1 $VFN1 \
				-vfn2 $VFN2
		fi
	else
		echo -e "\e[031m${My_bam%.bam}.RG.bam not found\e[0m\n"
		exit
	fi
	echo -e "\e[32m$SLX.$Barcode done for BisSNP genotyping\e[0m" 
fi

## Sort VCF before 'VCFpostprocess'
if [ $RUN_BS_SORT -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/scratch
	mkdir_unless $PROJECT_DIR/scratch/$Barcode
	if [ -s $VFN1 ]; then
		echo -e "perl $BIS_SNP_SORTER $VFN1\n"
		perl $BIS_SNP_SORTER \
			--k 1 --c 2 \
			--tmp $PROJECT_DIR/scratch/$Barcode \
			$VFN1 \
			$GENOME_DIR/hg19b.fa.fai > ${VFN1%.vcf}.sorted.vcf &
	else
		echo -e "$VFN1 not found"
	fi

	if [ -s $VFN2 ]; then
		echo -e "perl $BIS_SNP_SORTER $VFN2\n"
		perl $BIS_SNP_SORTER \
			--k 1 --c 2 \
			--tmp $PROJECT_DIR/scratch/$Barcode \
			$VFN2 \
			$GENOME_DIR/hg19b.fa.fai > ${VFN2%.vcf}.sorted.vcf &
	else
		echo -e "$VFN2 not found"
	fi
	wait
	echo -e "\e[32m$SLX.$Barcode done for BisSNP vcf sorting\e[0m" 
fi

##################################################
## 6. Filter Fake SNPs
## output: snp.filtered.vcf and cpg.filtered.vcf
## VCFpostprocess does not support "-nt"
## By default, it filter out SNPs with 
## 1) Genotype quality score filter for heterozygous SNP (less than 20), (-qual or --genotype qual) 
## 2) reads coverage more than 120, (-maxCov or --max coverage)
## 3) strand bias more than -0.02, (-sb or --strand bias)
## 4) quality score by depth less than 1.0, (-qd or --quality by depth)
## 5) mapping quality zero reads fraction more than 0.1 and (-mq0 or --mapping quality zero)
## 6) 2 SNPs within the same 20 bp window (-minSNPinWind)
##################################################
if [ $RUN_BS_FILTER -eq 1 ]; then
	# SNP.vcf
	if [ -s ${VFN2%.vcf}.sorted.vcf ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T VCFpostprocess -oldVcf ${VFN2%.vcf}.sorted.vcf"
		if [ $IS_FLD -eq 1 ]; then
			time java -Xmx$BIS_MEM -jar $BIS_SNP \
				-T VCFpostprocess \
				-R $REF_SEQ \
				-L $TARGET \
				-maxCov $BIS_MAXC \
				-oldVcf ${VFN2%.vcf}.sorted.vcf \
				-newVcf ${VFN2%.raw.vcf}.filtered.vcf \
				-snpVcf ${VFN2%.vcf}.sorted.vcf \
				-o ${VFN2%.raw.vcf}.filtered.summary.txt
		else
			time java -Xmx$BIS_MEM -jar $BIS_SNP \
				-T VCFpostprocess \
				-R $REF_SEQ \
				-maxCov $BIS_MAXC \
				-oldVcf ${VFN2%.vcf}.sorted.vcf \
				-newVcf ${VFN2%.raw.vcf}.filtered.vcf \
				-snpVcf ${VFN2%.vcf}.sorted.vcf \
				-o ${VFN2%.raw.vcf}.filtered.summary.txt
		fi
	else
		echo -e "${VFN2%.vcf}.sorted.vcf not found"
	fi
	# CpG.vcf
	if [ -s ${VFN1%.vcf}.sorted.vcf ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T VCFpostprocess -oldVcf ${VFN1%.vcf}.sorted.vcf"
		if [ $IS_FLD -eq 1 ]; then
			time java -Xmx$BIS_MEM -jar $BIS_SNP \
				-T VCFpostprocess \
				-R $REF_SEQ \
				-L $TARGET \
				-maxCov $BIS_MAXC \
				-oldVcf ${VFN1%.vcf}.sorted.vcf \
				-newVcf ${VFN1%.raw.vcf}.filtered.vcf \
				-snpVcf ${VFN2%.vcf}.sorted.vcf \
				-o ${VFN1%.raw.vcf}.filtered.summary.txt
		else
			time java -Xmx$BIS_MEM -jar $BIS_SNP \
				-T VCFpostprocess \
				-R $REF_SEQ \
				-maxCov $BIS_MAXC \
				-oldVcf ${VFN1%.vcf}.sorted.vcf \
				-newVcf ${VFN1%.raw.vcf}.filtered.vcf \
				-snpVcf ${VFN2%.vcf}.sorted.vcf \
				-o ${VFN1%.raw.vcf}.filtered.summary.txt
		fi
	else
		echo -e "${VFN1%.vcf}.sorted.vcf not found"
	fi 
	wait
	echo -e "\e[32m$SLX.$Barcode done for BisSNP vcf filtering\e[0m" 
fi

echo -e "\e[32m$SLX.$Barcode done for BisSNP\e[0m" 

##
## 7. Generate bed file or wig file for SNP/DNA methylation visualization 
