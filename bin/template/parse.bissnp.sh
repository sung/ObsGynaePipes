#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Jul/2014
# Last modified: 15/Jul/2014 
# Optimised and customised to run at the Darwin HPC

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

export SLX="MY_SLX" # e.g. SLX-8080 
export PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
export Barcode="MY_BARCODE" # e.g. A001

if [ $IS_PE -eq 1 ]; then
	# e.g. SLX-8075.A006.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.CpG_report.txt
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.bam 
	CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.CpG_report.txt
else
	# PNAS-2013.SRR530647.r_1_trimmed.fq.gz_bismark_bt2.q1.deduplicated.bam
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.bam 
	CpG_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.CpG_report.txt
fi
BS_CpG_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.cpg.filtered.vcf 
BS_SNP_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.snp.filtered.vcf 

#########
# OUTPUT
# Depth of Coverage * Breadth of Coverage
CpG_COV=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.CpG.coverage.txt
BS_CpG_COV=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.CpG.coverage.txt

# Percentage of methylated CpG
CpG_PCT_report=$PROJECT_DIR/Bismark/MethylCall/$Barcode/$SLX.$Barcode.CpG.methylation.percentage.txt
BS_CpG_PCT_report=$PROJECT_DIR/BisSNP/$Barcode/$SLX.$Barcode.CpG.methylation.percentage.txt

#############
## BIS-SNP ##
#############
if [ -s $BS_CpG_report ];then
	############################################
	## 1. Convert VCF to TAB format for mysql ##
	############################################
	if [ $RUN_FORMAT_BISSNP_CPG2SQL -eq 1 ]; then
		echo -e "making ${BS_CpG_report%.vcf}.tab.txt from BisSNP"
		# 3. Homo/Hetero merged CG
		time grep ^chr $BS_CpG_report | awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($8~"CS=-;"){str="R"}}{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>0){if(fm[5]=="CG" || fm[5]=="CR" || fm[5]=="CS" || fm[5]=="CK" || fm[5]=="MG" || fm[5]=="SG" || fm[5]=="YG"){cov=fm[4]+fm[6];print $1"."$2,$1,$2,str,cov,fm[4]/cov*100,fm[6]/cov*100}}}' > ${BS_CpG_report%.vcf}.CG.methylkit.txt

		# 2. Hetero CG (CR, CS, CK, MG?, SG?, YG?)
		time grep ^chr $BS_CpG_report | awk -F'\t' "BEGIN{OFS=FS}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>0){if(fm[5]==\"CR\" || fm[5]==\"CS\" || fm[5]==\"CK\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpG_report%.vcf}.het.CG.bed

		echo -e "\e[32m$SLX.$Barcode done for vcf2bed (based on Bis-SNP)\e[0m" 
	fi
else
	echo -e "\e[031m$BS_CpG_report not found even after merging per chromosome. Merging failed?\e[0m\n"
	exit
fi # end of BS_CpG_report
