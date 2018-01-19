#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 15/Jul/2014
# Last modified: 19/Feb/2016
# Optimised and customised to run at the Darwin HPC
# CpH only (CHH & CHG called by BisSNP)

source $HOME/Methyl-Seq/config/methyl_seq.config # export envrionment variables
source $HOME/lib/sung.sh #defines 'mkdir_unless'

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
Barcode="MY_BARCODE" # e.g. A001
Chr="MY_CELL" # e.g. chr1 

mkdir_unless $PROJECT_DIR	

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Cell.$Barcode.$Chr\e[0m\n"

###########
## Input ##
###########
BS_CpH_report=$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.non.cpg.filtered.vcf 

############
## OUTPUT ##
############
BS_CpH_COV=$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.CpH.coverage.txt
BS_CpH_PCT_report=$PROJECT_DIR/BisSNP/$Barcode/$Chr/$SLX.$Barcode.$Chr.CpH.methylation.percentage.txt

#######################################
## 5.1 CpG Coverage based on Bis-SNP ##
#######################################
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

	# CpH
	if [ -s $BS_CpH_report ]; then
		# 4. Non-CG context (CY, CW, CM, CH, MH?, MR?, MS?, MY?, SH?, SR?, YH?, YK?, YM?, YR?, YS?, YY?)
		## CpH should be called using non-default BisSNP options
		## -C CHH,1 -C CHG,1 -out_modes EMIT_ALL_CYTOSINES
		echo -e "Non-CG"
		time grep -v ^# $BS_CpH_report | awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($8~"CS=-;"){str="R"}}{split($10,fm,":")}{if($7 =="PASS" && fm[4]+fm[6]>$MIN_DEPTH){if(fm[5]=="CHH" || fm[5]=="CHG"){cov=fm[4]+fm[6];print $1"."$2,$1,$2,str,cov,fm[4]/cov*100,fm[6]/cov*100}}}' > ${BS_CpH_report%.vcf}.methylkit.txt
	fi

	echo -e "\e[32m$SLX.$Barcode done for vcf2methylkit (based on Bis-SNP)\e[0m" 
fi

##################################
## 8. Convert VCF to BED format ## 
##################################
if [ $RUN_FORMAT_VCF2BED_BISSNP -eq 1 ]; then  # convert VCF to MethylKit input
	if [ -s $BS_CpH_report ]; then
		# 4. Non-CG context (CHH, CHG)
		## CpH should be called using non-default BisSNP options
		## -C CHH,1 -C CHG,1 -out_modes EMIT_ALL_CYTOSINES
		#echo -e "Non-CG (CHH and CHG)"
		#time grep -v ^# $BS_CpH_report | awk -F'\t' "BEGIN{OFS=FS; print \"track name=$SLX.$Barcode.non.cpg.filtered.bed type=bedDetail description='non-CpG methylation level' visibility=3\"}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>$MIN_DEPTH){if(fm[5]==\"CHH\" || fm[5]==\"CHG\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpH_report%.vcf}.bed
		echo -e "Non-CG (CHH)"
		time grep -v ^# $BS_CpH_report | awk -F'\t' "BEGIN{OFS=FS; print \"track name=$SLX.$Barcode.non.cpg.filtered.bed type=bedDetail description='non-CpG (CHH) methylation level' visibility=3\"}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>$MIN_DEPTH){if(fm[5]==\"CHH\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpH_report%.vcf}.CHH.bed
		echo -e "Non-CG (CHG)"
		time grep -v ^# $BS_CpH_report | awk -F'\t' "BEGIN{OFS=FS; print \"track name=$SLX.$Barcode.non.cpg.filtered.bed type=bedDetail description='non-CpG (CHG) methylation level' visibility=3\"}{split(\$10,fm,\":\"); split(fm[3],brc,\",\")}{if(\$7 ==\"PASS\" && fm[4]+fm[6]>$MIN_DEPTH){if(fm[5]==\"CHG\"){if(\$8~\"CS=-;\") str=\"-\"; else str=\"+\"; chr=\$1; start=\$2-1; end=\$2; cov=fm[4]+fm[6]; met=fm[4]/cov; numA=brc[5]; numG=brc[4]; if(met<=0.1) rgb=\"0,240,0\"; else if(met<=0.2 && met>0.1) rgb=\"30,210,0\"; else if(met<=0.3 &&met>0.2) rgb=\"60,180,0\"; else if(met<=0.5 && met>0.3) rgb=\"90,150,0\"; else if(met<=0.6 && met>0.5) rgb=\"120,120,0\"; else if(met<=0.7 && met>0.6) rgb=\"150,90,0\"; else if(met<=0.8 && met>0.7) rgb=\"180,60,0\"; else rgb=\"210,0,0\"; print chr,start,end,met*100,cov,str,start,end,rgb,numA,numG}}}" > ${BS_CpH_report%.vcf}.CHG.bed
	fi

	echo -e "\e[32m$SLX.$Barcode done for vcf2bed (based on Bis-SNP)\e[0m" 
fi

echo -e "\e[32m$SLX.$Barcode all done for cpg.coverage.report\e[0m" 
