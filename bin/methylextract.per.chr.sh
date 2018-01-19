#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Reference: http://bioinfo2.ugr.es/MethylExtract/Manual.html
# Called by bin/template/trim_bismark_template.v{n}.sh
# First created: 1/Jul/2014
# Last modified: 8/Aug/2014 

# input
if [ $IS_PE -eq 1 ]; then
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q$MAPQ.deduplicated.$Chr.bam 
else
	My_bam=$PROJECT_DIR/Bismark/Alignment/$Barcode/$SLX.$Barcode.r_1_trimmed.fq.gz_bismark_bt2.q$MAPQ.deduplicated.$Chr.bam 
fi
if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam not found\e[0m\n"
	exit
fi

mkdir_unless $PROJECT_DIR/MethylExtract/$Barcode/$Chr
##################################################################################
## 13-1. MethylExtract
## Input: dedup bam from bismark
## Output: 
##################################################################################
if [ $RUN_METHYL_EXT_MAIN -eq 1 ]; then
	# if there is no sym link
	#if [ ! -h ${My_bam/Bismark\/Alignment/MethylExtract} ];then
	if [ ! -h ${My_bam/Bismark\/Alignment\/$Barcode/MethylExtract\/$Barcode\/$Chr} ];then
		ln -s $My_bam ${My_bam/Bismark\/Alignment\/$Barcode/MethylExtract\/$Barcode\/$Chr}
	fi

	if [ -h ${My_bam/Bismark\/Alignment\/$Barcode/MethylExtract\/$Barcode\/$Chr} ];then
		echo -e "MethylExtract.pl $SLX.$Barcode.$Chr"
		time echo | perl $ME_DIR/MethylExtract.pl \
			seq=$REF_SEQ \
			inDir=$PROJECT_DIR/MethylExtract/$Barcode/$Chr \
			outDir=$PROJECT_DIR/MethylExtract/$Barcode/$Chr \
			flagW=$ME_FLAGW flagC=$ME_FLAGC \
			peOverlap=$ME_PE_OVERLAP \
			minDepthSNV=$ME_MIN_SNV \
			minDepthMeth=$ME_MIN_MET \
			p=$NT \
			bedOut=Y \
			context=ALL > /dev/null
	else
		echo -e "\e[031mSymbolic link ${My_bam/Bismark\/Alignment\/$Barcode/MethylExtract\/$Barcode\/$Chr} not found\e[0m\n"
		exit
	fi
fi # end of RUN_METHYL_EXT -eq 1 

###############################
## 13-2. MethylExtractBSCR
## Input: sam.gz 
## Output: MethylExtractBSCR.result.txt 
## time: 
###############################
if [ $RUN_METHYL_EXT_BSCR -eq 1 ]; then
	## convert bam to sam.gz 
	echo -e "samtools view -h $My_bam | gzip > ${My_bam%.bam}.sam.gz"
	time samtools view -h $My_bam | gzip > ${My_bam%.bam}.sam.gz
	if [ -s ${My_bam%.bam}.sam.gz ]; then
		echo -e "perl MethylExtractBSCR.pl $SLX.$Barcode.$Chr"
		time perl $ME_DIR/MethylExtractBSCR.pl \
			seqFile=$REF_SEQ \
			inFile=${My_bam%.bam}.sam.gz \
			flagW=$ME_FLAGW flagC=$ME_FLAGC > \
			$PROJECT_DIR/MethylExtract/$Barcode/$Chr/$SLX.$Barcode.$Chr.MethylExtractBSCR.result.txt
	else
		echo -e "\e[031m ${My_bam%.bam}.sorted.sam.gz not found\e[0m\n"
		exit
	fi
	time rm ${My_bam%.bam}.sorted.sam.gz
fi

#################################
## 13-3. MethylExtractBSPvalue
#################################
if [ $RUN_METHYL_EXT_BSP -eq 1 ]; then
	if [ -s $PROJECT_DIR/MethylExtract/$Barcode/$Chr/$SLX.$Barcode.$Chr.MethylExtractBSCR.result.txt ]; then
		time perl $ME_DIR/MethylExtractBSPvalue.pl \
		inFile= \
		BSCR= \
		outFile= \
		errorInterval= \
		FDR=
	else
		echo -e "\e[031m $PROJECT_DIR/MethylExtract/$Barcode/$Chr/$SLX.$Barcode.$Chr.MethylExtractBSCR.result.txt not found\e[0m\n"
		exit
	fi
fi

