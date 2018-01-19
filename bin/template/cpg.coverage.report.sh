#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Jul/2014
# Last modified: 16/Apr/2015
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

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

###########
## Input ##
###########
if [ $IS_PE -eq 1 ]; then
	if [ $IS_FLD -eq 1 ]; then
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.sorted.bam 
		CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.CpG_report.txt
	else
		My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.sorted.bam 
		CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.CpG_report.txt
	fi
else
	# PNAS-2013.SRR530647.r_1_trimmed.fq.gz_bismark_bt2.q1.deduplicated.bam
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.sorted.bam 
	CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.CpG_report.txt
fi
BS_CpG_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.cpg.filtered.vcf 
BS_CpH_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.non.cpg.filtered.vcf 
BS_SNP_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.snp.filtered.vcf 

############
## OUTPUT ##
############
# Depth of Coverage & Breadth of Coverage
CpG_COV=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.CpG.coverage.txt
BS_CpG_COV=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.CpG.coverage.txt
BS_HOMO_CpG_COV=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.homo.CpG.coverage.txt
BS_CpH_COV=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.CpH.coverage.txt
# Percentage of methylated CpG
CpG_PCT_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.CpG.methylation.percentage.txt
BS_CpG_PCT_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.CpG.methylation.percentage.txt
BS_CpH_PCT_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.CpH.methylation.percentage.txt

###################################################
# 1. Coverage Report based on the sorted bam file #
# Input: sorted bam
# Output: $Barcode.genomeccov.txt (genome coverage report by bedtool)
# Output: $Barcode.coverage.txt (DoC and BoC)
###################################################
if [ $RUN_GENOME_COV -eq 1 ]; then
	echo -e "\e[32m$SLX.$Barcode Genome (or Target) Coverage\e[0m" 
	mkdir_unless $PROJECT_DIR/Coverage
	mkdir_unless $PROJECT_DIR/Coverage/$Barcode

	if [ ! -s $My_bam ];then
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi

	# targeted sequencing (i.e. Fluidigm)
	if [ -n "$TARGET" ];then
		# /Bismakr/Alignment => /Coverage	
		My_coverageBed=${My_bam/Bismark\/Alignment/Coverage}
		### CoverageStas [ PCR/optical duplicateReads were removed]
		echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET > $My_coverageBed.CoverageStats \n"
		time samtools view -u $My_bam | coverageBed -abam stdin -b $TARGET > $My_coverageBed.CoverageStats

		### CoveragePerBase [ PCR/optical duplicateReads were removed] -d: perbase
		echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET -d > $My_coverageBed.CoveragePerBase \n"
		time samtools view -u $My_bam | coverageBed -abam stdin -b $TARGET -d > $My_coverageBed.CoveragePerBase

		### CoverageHistogram [ PCR/optical duplicateReads were removed]
		echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET -hist > $My_coverageBed.CoverageHist \n"
		time samtools view -u $My_bam | coverageBed -abam stdin -b $TARGET -hist > $My_coverageBed.CoverageHist

		# per-amplicon
		if [ -n "$TARGET2" ];then
			echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET2 -s > $My_coverageBed.CoverageStats.per.amplicon \n"
			time samtools view -u $My_bam | coverageBed -abam stdin -b $TARGET2 -s > $My_coverageBed.CoverageStats.per.amplicon

			echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET2 -d -s > $My_coverageBed.CoveragePerBase.per.amplicon \n"
			time samtools view -u $My_bam | coverageBed -abam stdin -b $TARGET2 -d -s > $My_coverageBed.CoveragePerBase.per.amplicon

			echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET2 -hist -s > $My_coverageBed.CoverageHist.per.amplicon \n"
			time samtools view -u $My_bam | coverageBed -abam stdin -b $TARGET2 -hist -s > $My_coverageBed.CoverageHist.per.amplicon
		fi
	# WGBS
	else
		# e.g. SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.sorted.bam
		# get the genome size based on the BAM file
		genomeSize=`samtools view -H $My_bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'`
		Genome_Cov=$PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.bedtools.genomecov.txt

		echo -e "bedtools genomecov -ibam $My_bam > $Genome_Cov"
		time bedtools genomecov -ibam $My_bam > $Genome_Cov

		if [ -s $Genome_Cov ];then
			# Genome Depth of Coveage =
			# sum(depth*num-of-base) / genome-size
			time grep -P '^genome' $Genome_Cov | awk "{sum+=\$2*\$3} END { print \"depth of coverage= \",sum/$genomeSize}" > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.coverage.txt
			for d in ${DEPTH[*]}
			do
				echo "calculating x$d coverage..."
			# Genome Breadh of Coverage at depth-x (%-genome-covered-at-depth-x) =
			# sum(num-of-base-more-than-depth-x) / genome-size
				time grep -P '^genome' $Genome_Cov | awk "{if(\$2>=$d){sum+=\$3}} END {print \"x$d BoC=\",sum/$genomeSize*100}" >> $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.coverage.txt
			done
		else
			echo -e "\e[031m$Genome_Cov not found\e[0m\n"
		fi
	fi
	echo -e "\e[32m$SLX.$Barcode done for Genome (or Target) Coverage\e[0m" 
fi # end of RUN_GENOME_COV

###########################################
## Re-format Bismark to BED or methylkit ## 
###########################################
if [ $RUN_FORMAT_BED_BISMARK -eq 1 ];then
	echo -e "\e[32m$SLX.$Barcode Bismark 2 BED \e[0m" 
	if [ -s $CpG_report ];then
		# Bismakr CpG-report to bed format
		time awk 'BEGIN{OFS="\t"}{depth=$4+$5;if(depth==0)met=0;else met=$4/depth*100;print $1,$2-1,$2,met,depth,$3}' $CpG_report > ${CpG_report%.CpG_report.txt}.bed
		if [ -n "$TARGET" ];then
			# sort by position
			# time sort -k1,1 -k2,2n ${CpG_report%.CpG_report.txt}.bed > ${CpG_report%.CpG_report.txt}.sorted.bed 
			# intersect bed file by target regions
			time bedtools intersect -a ${CpG_report%.CpG_report.txt}.bed -b $TARGET > ${CpG_report%.CpG_report.txt}.ontarget.bed 
			# sort by position
			time sort -k1,1 -k2,2n ${CpG_report%.CpG_report.txt}.ontarget.bed > ${CpG_report%.CpG_report.txt}.ontarget.sorted.bed 
		fi
	else
		echo -e "\e[031m$CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
		exit
	fi # end of $CpG_report
	echo -e "\e[32m$SLX.$Barcode done for converting Bismark to BED\e[0m" 
fi # end of RUN_FORMAT_BED_BISMARK

# [todo] Bismark CpG calls in methylkit format?
if [ $RUN_FORMAT_METHYLKIT_BISMARK -eq 1 ];then
	echo -e "\e[32m$SLX.$Barcode Bismark CpG % Methylation\e[0m" 
	if [ -s $CpG_report ];then
		# make a file for MethylKit format
		echo -e "making ${CpG_report%.txt}.methylkit.txt from Bismark"
		time awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($3=="-"){str="R"}}{if($4+$5>0){cov=$4+$5;print $1"."$2,$1,$2,str,cov,$4/cov*100,$5/cov*100}}' $CpG_report > ${CpG_report%.txt}.methylkit.txt
	else
		echo -e "\e[031m$CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
		exit
	fi # end of $CpG_report
	echo -e "\e[32m$SLX.$Barcode done for converting Bismark to MethylKit\e[0m" 
fi # end of RUN_FORMAT_METHYLKIT_BISMARK


#############################################
## 2. CpG Coverage based on Bismark Caller ##
#############################################
if [ $RUN_CPG_COV_BISMARK -eq 1 ]; then
	echo -e "\e[32m$SLX.$Barcode Bismark CpG Depth of Coverage\e[0m" 
	###################################################
	## if bismark methyl extract was done per sample ##
	###################################################
	if [ -s $CpG_report ];then
		if [ -n "$TARGET" ];then
			# Depth of CpG Coveage = (NO of methylated C + NO of unmetnylated C) / Total NO of CpG
			time awk "{sum+=\$5}END{print \"CpG Depth of Cov=\",sum/NR,sum,NR}" ${CpG_report%.CpG_report.txt}.ontarget.sorted.bed > $CpG_COV
			for d in ${DEPTH[*]}
			do
				echo "calculating x$d coverage..."
				# Depth of CpG Coveage = (NO of methylated C + NO of unmetnylated C) / Total NO of CpG
				time awk "{if(\$5>=$d){sum++}}END{print \"x$d CpG BoC=\",sum/NR*100,sum,NR}" ${CpG_report%.CpG_report.txt}.ontarget.sorted.bed >> $CpG_COV
			done
		else
			# Depth of CpG Coveage = (NO of methylated C + NO of unmetnylated C) / Total NO of CpG
			time awk "{sum+=\$4+\$5}END{print \"CpG Depth of Cov=\",sum/NR,sum,NR}" $CpG_report > $CpG_COV
			for d in ${DEPTH[*]}
			do
				echo "calculating x$d coverage..."
				# Breadth of CpG Coverge at depth x = 
				# NO of CpG called at depth x / Total NO of CpG  
				time awk "{if(\$4+\$5>=$d){sum++}}END{print \"x$d CpG BoC=\",sum/NR*100,sum,NR}" $CpG_report >> $CpG_COV
			done
		fi
	#######################################################
	## or bismark methyl extract was done per chromosome ##
	## Then, there would not be $CpG_report
	## So, make this by merging per chromosome CpG_report
	#######################################################
	else
		for Chr in ${UCSC_CHR[*]}
		do
			if [ $IS_PE -eq 1 ]; then
				#SLX-8080.A001.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.chr1.sorted.CpG_report.txt
				Chr_CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$Chr/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.sorted.CpG_report.txt
			else
				Chr_CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$Chr/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.$Chr.sorted.CpG_report.txt
			fi
			if [ -s $Chr_CpG_report ];then
				if [ $Chr = 'chr1' ]; then
					echo -e "scrapping $Chr"...
					awk "{if(\$1==\"$Chr\"){print \$0}}" $Chr_CpG_report > $CpG_report
				else
					echo -e "scrapping $Chr"...
					awk "{if(\$1==\"$Chr\"){print \$0}}" $Chr_CpG_report >> $CpG_report
				fi
			else
				echo -e "\e[031m$Chr_CpG_report not found\e[0m\n"
			fi
		done # end of Chr 
		#################################
		## now I am expecting $CpG_report 
		#################################
		if [ -s $CpG_report ];then
			if [ -n "$TARGET" ];then
				# Bismark CpG-report to bed format
				time awk 'BEGIN{OFS="\t"}{depth=$4+$5;if(depth==0)met=0;else met=$4/depth*100;print $1,$2-1,$2,met,depth,$3}' $CpG_report > ${CpG_report%.CpG_report.txt}.bed
				# sort by position
				#time sort -k1,1 -k2,2n ${CpG_report%.CpG_report.txt}.bed > ${CpG_report%.CpG_report.txt}.sorted.bed 
				# intersect bed file by target regions
				time bedtools intersect -a ${CpG_report%.CpG_report.txt}.bed -b $TARGET > ${CpG_report%.CpG_report.txt}.ontarget.bed 
				# sort by position
				time sort -k1,1 -k2,2n ${CpG_report%.CpG_report.txt}.ontarget.bed > ${CpG_report%.CpG_report.txt}.ontarget.sorted.bed 
				# Depth of CpG Coveage = (NO of methylated C + NO of unmetnylated C) / Total NO of CpG
				time awk "{sum+=\$5}END{print \"CpG Depth of Cov=\",sum/NR,sum,NR}" ${CpG_report%.CpG_report.txt}.ontarget.sorted.bed > $CpG_COV
				for d in ${DEPTH[*]}
				do
					echo "calculating x$d coverage..."
					# Depth of CpG Coveage = (NO of methylated C + NO of unmetnylated C) / Total NO of CpG
					time awk "{if(\$5>=$d){sum++}}END{print \"x$d CpG BoC=\",sum/NR*100,sum,NR}" ${CpG_report%.CpG_report.txt}.ontarget.sorted.bed >> $CpG_COV
				done
			else
				# Depth of CpG Coveage = (NO of methylated C + NO of unmetnylated C) / Total NO of CpG
				time awk "{sum+=\$4+\$5}END{print \"CpG Depth of Cov=\",sum/NR,sum,NR}" $CpG_report > $CpG_COV
				for d in ${DEPTH[*]}
				do
					echo "calculating x$d coverage..."
					time awk "{if(\$4+\$5>=$d){sum++}}END{print \"x$d CpG BoC=\",sum/NR*100,sum,NR}" $CpG_report >> $CpG_COV
				done
			fi
		else
			echo -e "\e[031m$CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
			exit
		fi # end of $CpG_report
	fi # end of $CpG_report
	echo -e "\e[32m$SLX.$Barcode done for Bismark CpG Depth of Coverage\e[0m" 
fi # end of RUN_CPG_COV_BISMARK 

##############################################################
## 3. Genome-wide % CpG methylation based on Bismark Caller ##
##############################################################
if [ $RUN_CPG_MET_PCT_BISMARK -eq 1 ]; then
	if [ -s $CpG_report ];then
		for d in ${DEPTH[*]}
		do
			echo "computing % CpG methylation at x$d coverage..."
			# % methylation = 
			# NO of methylated-C / NO of total CpG (from hg19)
			if [ $d -eq 1 ]; then
				time awk  "{if(\$4>0 && \$4+\$5>=$d){sum++}}END{print \"x$d % met-CpG=\",sum/NR*100,sum,NR}" $CpG_report > $CpG_PCT_report
			else
				time awk  "{if(\$4>0 && \$4+\$5>=$d){sum++}}END{print \"x$d % met-CpG=\",sum/NR*100,sum,NR}" $CpG_report >> $CpG_PCT_report
			fi
			# % methylation = 
			# NO of methylated-C / NO of CpG called 
			time awk  "{if(\$4+\$5>=$d){total++}}{if(\$4>0 && \$4+\$5>=$d){sum++}}END{print \"x$d % met-called-CpG=\",sum/total*100,sum,total}" $CpG_report >> $CpG_PCT_report
		done
	else
		echo -e "\e[031m$CpG_report not found\e[0m\n"
		exit
	fi
	echo -e "\e[32m$SLX.$Barcode done for CpG % Methylation \e[0m" 
fi # end of RUN_CPG_MET_PCT

################################
## 4. MERGE BIS-SNP VCF FILES ##
## MADE BY BIS-SNP            ##
## SEE README AS WELL         ##
################################
if [ $RUN_MERGE_BS_VCF -eq 1 ]; then
	for Chr in ${UCSC_CHR[*]}
	do
		#SLX-8074.SLX-8080.A012.chrM.cpg.filtered.vcf
		#SLX-8074.SLX-8080.A012.chrM.snp.filtered.vcf
		Chr_CpG_report=$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.cpg.filtered.vcf
		Chr_SNP_report=$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.snp.filtered.vcf
		if [ -s $Chr_CpG_report ];then
			echo -e "scrapping $Chr"...
			if [ $Chr = 'chr1' ]; then
				cat $Chr_CpG_report > $BS_CpG_report 
				cat $Chr_SNP_report > $BS_SNP_report
			else
				grep -v ^# $Chr_CpG_report >> $BS_CpG_report 
				grep -v ^# $Chr_SNP_report >> $BS_SNP_report
			fi
		else
			echo -e "\e[031m$Chr_CpG_report not found\e[0m\n"
			exit
		fi
	done # end of Chr 
	echo -e "\e[32m$SLX.$Barcode done for Merging VCF (based on Bis-SNP)\e[0m" 
fi

#######################################
## 5.1 CpG Coverage based on Bis-SNP ##
#######################################
if [ $RUN_CPG_COV_BISSNP -eq 1 ]; then
	if [ ! -s $BS_CpG_report ];then
		echo -e "\e[031m$BS_CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
		exit
	fi
	# Depth of CpG Coveage = (NO of methylated C + NO of unmetnylated C) / NO of CpG called
	# Homo-CpG (CG)
	echo "calculating BisSNP Homo-CpG depth of coverage..."
	time grep -v ^# $BS_CpG_report | awk '{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>0){if(fm[5]=="CG"){num++;sum+=fm[4]+fm[6]}}}END{print "CpG Depth of Cov=",sum/num,sum,num}' > $BS_HOMO_CpG_COV
	NUM_CG=`awk '{print $7}' $BS_HOMO_CpG_COV`
	for d in ${DEPTH[*]}
	do
		# Breadth of CpG Coverge at depth x = NO of CpG called at depth x / NO of CpG called 
		echo "calculating x$d coverage..."
		time grep -v ^# $BS_CpG_report | awk "{split(\$10,fm,\":\")}{if(fm[5]==\"CG\"){if(\$7==\"PASS\" && fm[4]+fm[6]>=$d){sum++;}}}END{print \"x$d CpG BoC=\",sum/$NUM_CG*100,sum,$NUM_CG}" >> $BS_HOMO_CpG_COV 
	done
	echo -e "\e[32m$SLX.$Barcode done for Homo-CpG Depth of Cov (based on Bis-SNP)\e[0m" 

	# Homo/Hetero-CpG (CG, CR, CS, CK)
	echo "calculating BisSNP Homo/Hetero-CpG depth of coverage..."
	time grep -v ^# $BS_CpG_report | awk '{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>0){if(fm[5]=="CG" || fm[5]=="CR" || fm[5]=="CS" || fm[5]=="CK"){num++;sum+=fm[4]+fm[6]}}}END{print "CpG Depth of Cov=",sum/num,sum,num}' > $BS_CpG_COV
	NUM_CG=`awk '{print $7}' $BS_CpG_COV`
	for d in ${DEPTH[*]}
	do
		# Breadth of CpG Coverge at depth x = NO of CpG called at depth x / NO of CpG called 
		echo "calculating x$d coverage..."
		#time grep -v ^# $BS_CpG_report | awk "{if(\$10~\":CG:\"){split(\$10,fm,\":\");if(fm[4]+fm[6]>=$d){sum++}}}END{print \"x$d CpG BoC=\",sum/$NUM_CG*100,sum,$NUM_CG}" >> $BS_CpG_COV 
		time grep -v ^# $BS_CpG_report | awk "{split(\$10,fm,\":\")}{if(fm[5]==\"CG\" || fm[5]==\"CR\" || fm[5]==\"CS\" || fm[5]==\"CK\"){if(\$7==\"PASS\" && fm[4]+fm[6]>=$d){sum++;}}}END{print \"x$d CpG BoC=\",sum/$NUM_CG*100,sum,$NUM_CG}" >> $BS_CpG_COV 
	done
	echo -e "\e[32m$SLX.$Barcode done for Homo/Hetero-CpG Depth of Cov (based on Bis-SNP)\e[0m" 
	# Homo/Hetero-CpG
fi

#######################################
## 5.2 CpH Coverage based on Bis-SNP ##
## CpH should be called using non-default BisSNP options
## -C CHH,1 -C CHG,1 -out_modes EMIT_ALL_CYTOSINES
#######################################
if [ $RUN_CPH_COV_BISSNP -eq 1 ]; then
	if [ ! -s $BS_CpH_report ];then
		echo -e "\e[031m$BS_CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
		exit
	fi
	# Avg. Depth of CpH Coveage = 
	# (NO of methylated C + NO of unmetnylated C) / NO of CpH called
	echo "calculating BisSNP CpH average depth of coverage..."

	time grep -v ^# $BS_CpH_report | awk '{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>0){if(fm[5]=="CHH" || fm[5]=="CHG"){num++;sum+=fm[4]+fm[6]}}}END{print "CpH Depth of Cov=",sum/num,sum,num}' > $BS_CpH_COV

	NUM_CH=`awk '{print $7}' $BS_CpH_COV`
	for d in ${DEPTH[*]}
	do
		# Breadth of CpH Coverge at depth x = 
		# NO of CpH called at depth x / NO of CpH called 
		echo "calculating x$d coverage..."
		time grep -v ^# $BS_CpH_report | awk "{split(\$10,fm,\":\")}{if(fm[5]==\"CHH\" || fm[5]==\"CHG\"){if(\$7 ==\"PASS\" && fm[4]+fm[6]>=$d){sum++;}}}END{print \"x$d CpH BoC=\",sum/$NUM_CH*100,sum,$NUM_CH}" >> $BS_CpH_COV 
	done
	echo -e "\e[32m$SLX.$Barcode done for CpH Depth of Cov (based on Bis-SNP)\e[0m" 
fi


########################################################
## 6.1 Genome-wide % CpG methylation based on BIS-SNP ##
########################################################
if [ $RUN_CPG_MET_PCT_BISSNP -eq 1 ]; then
	if [ ! -s $BS_CpG_report ];then
		echo -e "\e[031m$BS_CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
		exit
	fi
	# % methylation at-depth-x = sum-of-met-c-at-depth-x / sum-of-cpg-at-depth-x
	for d in ${DEPTH[*]}
	do
		echo "computing BisSNP % CpG methylation at x$d coverage..."
		if [ $d -eq 1 ]; then
			time grep -v ^# $BS_CpG_report | awk "{split(\$10,fm,\":\")}{if(\$7==\"PASS\" && fm[4]+fm[6]>=$d){if(fm[5]==\"CG\"){sum+=fm[4]+fm[6];met+=fm[4]}}}END{print \"x$d % met-called-CpG=\",met/sum*100,met,sum}" > $BS_CpG_PCT_report
		else
			time grep -v ^# $BS_CpG_report | awk "{split(\$10,fm,\":\")}{if(\$7==\"PASS\" && fm[4]+fm[6]>=$d){if(fm[5]==\"CG\"){sum+=fm[4]+fm[6];met+=fm[4]}}}END{print \"x$d % met-called-CpG=\",met/sum*100,met,sum}" >> $BS_CpG_PCT_report
		fi
	done
	echo -e "\e[32m$SLX.$Barcode done for CpG % Methylation (based on Bis-SNP)\e[0m" 
fi # end of RUN_CPG_MET_PCT

########################################################
## 6.2 Genome-wide % CpH methylation based on BIS-SNP ##
## CpH should be called using non-default BisSNP options
## -C CHH,1 -C CHG,1 -out_modes EMIT_ALL_CYTOSINES
########################################################
if [ $RUN_CPH_MET_PCT_BISSNP -eq 1 ]; then
	if [ ! -s $BS_CpH_report ];then
		echo -e "\e[031m$BS_CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
		exit
	fi
	# % methylation = NO of methylated-C / NO of CpH called 
	for d in ${DEPTH[*]}
	do
		echo "computing BisSNP % CpH methylation at x$d coverage..."
		if [ $d -eq 1 ]; then
			time grep -v ^# $BS_CpH_report | awk "{split(\$10,fm,\":\")}{if(\$7==\"PASS\" && fm[4]+fm[6]>=$d){if(fm[5]==\"CHH\" || fm[5]==\"CHG\"){sum+=fm[4]+fm[6];met+=fm[4]}}}END{print \"x$d % met-called-CpH=\",met/sum*100,met,sum}" > $BS_CpH_PCT_report
		else
			time grep -v ^# $BS_CpH_report | awk "{split(\$10,fm,\":\")}{if(\$7==\"PASS\" && fm[4]+fm[6]>=$d){if(fm[5]==\"CHH\" || fm[5]==\"CHG\"){sum+=fm[4]+fm[6];met+=fm[4]}}}END{print \"x$d % met-called-CpH=\",met/sum*100,met,sum}" >> $BS_CpH_PCT_report
		fi
	done
	echo -e "\e[32m$SLX.$Barcode done for CpH % Methylation (based on Bis-SNP)\e[0m" 
fi # end of RUN_CPG_MET_PCT

########################################
## 7. Convert VCF to MethylKit format ## 
########################################
if [ $RUN_FORMAT_METHYLKIT_BISSNP -eq 1 ]; then
	# CS=-;Context=CG,C;DP=12;MQ0=0;NS=1;REF=0 ($8 INFO)
	# GT:BQ:BRC6:CM:CP:CU:DP:DP4:GP:GQ:SS ($9 FORMAT)
	# 0/0:22.6,22.0:0,1,0,10,0,0:0:CG:1:12:10,0,0,1:0,58,248:58.47:5 ($10)
	# CM = num_c
	# CU = num_t
	# %methyl = num_c/(num_c+num_t)

	# CpG
	if [ -s $BS_CpG_report ]; then
		# 1. Homo CG
		#echo -e "homo-CG"
		time grep -v ^# $BS_CpG_report | awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($8~"CS=-;"){str="R"}}{split($10,fm,":")}{if($7 =="PASS" && fm[5] =="CG" && fm[4]+fm[6]>0){cov=fm[4]+fm[6];print $1"."$2,$1,$2,str,cov,fm[4]/cov*100,fm[6]/cov*100}}' > ${BS_CpG_report%.vcf}.homo.CG.methylkit.txt

		# 2. Hetero CG (CR, CS, CK, MG?, SG?, YG?)
		#echo -e "het-CG"
		time grep -v ^# $BS_CpG_report | awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($8~"CS=-;"){str="R"}}{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>0){if(fm[5]=="CR" || fm[5]=="CS" || fm[5]=="CK"){cov=fm[4]+fm[6];print $1"."$2,$1,$2,str,cov,fm[4]/cov*100,fm[6]/cov*100}}}' > ${BS_CpG_report%.vcf}.het.CG.methylkit.txt

		# 3. Homo/Hetero merged CG
		#echo -e "homo/het-CG"
		time grep -v ^# $BS_CpG_report | awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($8~"CS=-;"){str="R"}}{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>0){if(fm[5]=="CG" || fm[5]=="CR" || fm[5]=="CS" || fm[5]=="CK"){cov=fm[4]+fm[6];print $1"."$2,$1,$2,str,cov,fm[4]/cov*100,fm[6]/cov*100}}}' > ${BS_CpG_report%.vcf}.CG.methylkit.txt
	fi

	# CpH
	if [ -s $BS_CpH_report ]; then
		# 4. Non-CG context (CY, CW, CM, CH, MH?, MR?, MS?, MY?, SH?, SR?, YH?, YK?, YM?, YR?, YS?, YY?)
		## CpH should be called using non-default BisSNP options
		## -C CHH,1 -C CHG,1 -out_modes EMIT_ALL_CYTOSINES
		echo -e "Non-CG"
		time grep -v ^# $BS_CpH_report | awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($8~"CS=-;"){str="R"}}{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>0){if(fm[5]=="CHH" || fm[5]=="CHG"){cov=fm[4]+fm[6];print $1"."$2,$1,$2,str,cov,fm[4]/cov*100,fm[6]/cov*100}}}' > ${BS_CpH_report%.vcf}.methylkit.txt
	fi

	echo -e "\e[32m$SLX.$Barcode done for vcf2methylkit (based on Bis-SNP)\e[0m" 
fi

##################################
## 8. Convert VCF to BED format ## 
##################################
if [ $RUN_FORMAT_VCF2BED_BISSNP -eq 1 ]; then  # convert VCF to MethylKit input
	if [ -s $BS_CpG_report ]; then
		# make a file for MethylKit format
		echo -e "making ${BS_CpG_report%.vcf}.bed from BisSNP"

		# CS=-;Context=CG,C;DP=12;MQ0=0;NS=1;REF=0 ($8 INFO)
		# GT:BQ:BRC6:CM:CP:CU:DP:DP4:GP:GQ:SS ($9 FORMAT)
		# 0/0:22.6,22.0:0,1,0,10,0,0:0:CG:1:12:10,0,0,1:0,58,248:58.47:5 ($10)
		# CM = num_c
		# CU = num_t
		# %methyl = num_c/(num_c+num_t)

		# 1. Homo CG
		time grep -v ^# $BS_CpG_report | awk -F'\t' "BEGIN{OFS=FS; print \"track name=$SLX.$Barcode.cpg.filtered.homo.CG.bed type=bedDetail description='Homo-CpG methylation level' visibility=3\"}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>0){if(fm[5]==\"CG\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpG_report%.vcf}.homo.CG.bed

		# 2. Hetero CG (CR, CS, CK, MG?, SG?, YG?)
		time grep -v ^# $BS_CpG_report | awk -F'\t' "BEGIN{OFS=FS; print \"track name=$SLX.$Barcode.cpg.filtered.het.CG.bed type=bedDetail description='Hetero-CpG methylation level' visibility=3\"}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>0){if(fm[5]==\"CR\" || fm[5]==\"CS\" || fm[5]==\"CK\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpG_report%.vcf}.het.CG.bed

		# 3. Homo/Hetero merged CG
		time grep -v ^# $BS_CpG_report | awk -F'\t' "BEGIN{OFS=FS; print \"track name=$SLX.$Barcode.cpg.filtered.CG.bed type=bedDetail description='CpG methylation level' visibility=3\"}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>0){if(fm[5]==\"CG\" || fm[5]==\"CR\" || fm[5]==\"CS\" || fm[5]==\"CK\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpG_report%.vcf}.CG.bed
	fi

	if [ -s $BS_CpH_report ]; then
		# 4. Non-CG context (CHH, CHG)
		## CpH should be called using non-default BisSNP options
		## -C CHH,1 -C CHG,1 -out_modes EMIT_ALL_CYTOSINES
		echo -e "Non-CG"
		time grep -v ^# $BS_CpH_report | awk -F'\t' "BEGIN{OFS=FS; print \"track name=$SLX.$Barcode.non.cpg.filtered.bed type=bedDetail description='non-CpG methylation level' visibility=3\"}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>0){if(fm[5]==\"CHH\" || fm[5]==\"CHG\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpH_report%.vcf}.bed
	fi

	echo -e "\e[32m$SLX.$Barcode done for vcf2bed (based on Bis-SNP)\e[0m" 
fi

echo -e "\e[32m$SLX.$Barcode all done for cpg.coverage.report\e[0m" 
