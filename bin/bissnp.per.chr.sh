#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Reference: http://epigenome.usc.edu/publicationdata/bissnp2011/BisSNP-UserGuide-latest.pdf
# Called by bin/template/trim_bismark_template.v{n}.sh
# First created: 20/May/2014
# Last modified: 1/Aug/2014 

if [ $IS_PE -eq 1 ]; then
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.sorted.bam 
else
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.sorted.bam 
fi

if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam not found\e[0m\n"
	exit
fi

############################################################
# 1. Add RG (Read Group) to bam header
# 'CREATE_INDEX=true' will make unknown error - disabled
# output: my.RG.bam
############################################################
# should have been done by bin/bissnp.sh

#########################
## 2. Indel Realignment
#  If your mapping tools allows gapped alignment, like Bismark with Bowtie2 or Novoaligner, you could omit this step.
#########################
# should have been done by bin/bissnp.sh

#######################
## 3. Mark duplicates
# so skip this step if done by 'deduplicate_bismark' previously
#######################

###################################
## 4. Base Quality Recalibration
###################################
#4-1. Count Covariant
#4-2. Write recalibrated base quality score into BAM file
#4-3. Re-Count Covariant
#4-4. Generate recalibration plot
# should have been done by bin/bissnp.sh

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
## EMIT VARIANT AND CYTOSINES: emit all of Cytosine sites above emit threshold into vcf1 file, and all of SNP
## sites above emit threshold into vcf2 file.
##############################################
if [ $RUN_BS_GENO_BY_CHR -eq 1 ]; then
	# e.g. SLX-8074.SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.sorted.RG.recal.bam
	# e.g. SLX-8075.SLX-8077.SLX-8081.A002.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.sorted.RG.recal.bam
	if [ -s ${My_bam%.bam}.RG.recal.bam ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T BisulfiteGenotyper -I ${My_bam%.bam}.RG.recal.bam"
		time java -Xmx$BIS_MEM -jar $BIS_SNP \
			-T BisulfiteGenotyper \
			-nt $NT \
			--logging_level $BIS_LOG_LEVEL \
			-L $Chr \
			-R $REF_SEQ \
			-I ${My_bam%.bam}.RG.recal.bam \
			-D $DBSNP \
			-vfn1 $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.vcf \
			-vfn2 $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.vcf \
			-stand_call_conf $BIS_CALL_CONF \
			-stand_emit_conf $BIS_EMIT_CONF \
			-mmq $BIS_MMQ \
			-mbq $BIS_MBQ \
			-cgi $BIS_CGI_FILE \
			-coverage $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.coverage.txt 
			echo -e "\e[32m$SLX.$Barcode.$Chr done for BisSNP genotyping\e[0m" 
	else
		echo -e "\e[031m${My_bam%.bam}.RG.bam not found\e[0m\n"
		exit
	fi
fi

## Sort VCF before 'VCFpostprocess'
if [ $RUN_BS_SORT_BY_CHR -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/scratch
	mkdir_unless $PROJECT_DIR/scratch/$Barcode
	if [ -s $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.vcf ]; then
		perl $BIS_SNP_SORTER \
			--k 1 --c 2 \
			--tmp $PROJECT_DIR/scratch/$Barcode \
			$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.vcf \
			$GENOME_DIR/hg19b.fa.fai > $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.sorted.vcf &
	else
		echo -e "$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.vcf not found"
		exit
	fi

	if [ -s $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.vcf ]; then
		perl $BIS_SNP_SORTER \
			--k 1 --c 2 \
			--tmp $PROJECT_DIR/scratch/$Barcode \
			$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.vcf \
			$GENOME_DIR/hg19b.fa.fai > $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.sorted.vcf &
	else
		echo -e "$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.vcf not found"
		exit
	fi
	wait
	echo -e "\e[32m$SLX.$Barcode.$Chr done for BisSNP vcf sorting\e[0m" 
fi

##################################################
## 6. Filter Fake SNPs
## output: snp.filtered.vcf and cpg.filtered.vcf
## VCFpostprocess does not support "-nt"
##################################################
if [ $RUN_BS_FILTER_BY_CHR -eq 1 ]; then
	if [ -s $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.sorted.vcf ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T VCFpostprocess -oldVcf $Barcode.$Chr.snp.raw.vcf"
		time java -Xmx$BIS_MEM -jar $BIS_SNP \
			-T VCFpostprocess \
			-R $REF_SEQ \
			--logging_level $BIS_LOG_LEVEL \
			-oldVcf $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.sorted.vcf \
			-newVcf $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.filtered.vcf \
			-snpVcf $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.sorted.vcf \
			-o $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.filter.summary.txt &
	else
		echo -e "$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.sorted.vcf not found"
		exit
	fi

	if [ -s $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.vcf ]; then
		echo -e "java -Xmx$BIS_MEM -jar $BIS_SNP -T VCFpostprocess -oldVcf $Barcode.$Chr.cpg.raw.vcf"
		time java -Xmx$BIS_MEM -jar $BIS_SNP \
			-T VCFpostprocess \
			-R $REF_SEQ \
			--logging_level $BIS_LOG_LEVEL \
			-oldVcf $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.sorted.vcf \
			-newVcf $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.filtered.vcf \
			-snpVcf $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.raw.sorted.vcf \
			-o $PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.filter.summary.txt &
	else
		echo -e "$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.raw.sorted.vcf not found"
		exit
	fi 
	wait
	echo -e "\e[32m$SLX.$Barcode.$Chr done for BisSNP vcf filtering\e[0m" 
fi

##
## 7. Generate bed file or wig file for SNP/DNA methylation visualization 
