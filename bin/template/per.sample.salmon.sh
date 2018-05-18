#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 27/Apr/2017
# Last modified: 27/Apr/2017
# Optimised and customised to run at the Darwin HPC

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh 
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir -p $PROJECT_DIR/scratch/$Barcode

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample=$SLX.$Barcode\e[0m\n"

<<fastq
SLX-9169.D704_D501.C6A20ANXX.s_1.r_1.fq.gz
SLX-9169.D704_D501.C6GNPANXX.s_5.r_1.fq.gz
SLX-9169.D704_D501.C6GNPANXX.s_6.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_1.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_2.r_1.fq.gz
SLX-9169.D704_D501.C6GW5ANXX.s_3.r_1.fq.gz
fastq

# assuming single-end
# [todo] what if paired end?
if [ $IS_SE -eq 1 ]; then
# per each barcode (sample) 
# do 1)FastQC 2)trim 3)tophat 4)bedtool, 5)htseq 
	FastQ_array=() #initialise
	Trimmed_FastQ_array=() #initialise
	FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`
	for FastQ_file in $FASTQ_FILES	# per fastq file (per lane)
	do
        # trimmed FASTQ
		Trimmed_FastQ=${FastQ_file/data\/fastq\/$SLX/results\/$PROJECT\/Trim\/$Barcode}
		Trimmed_FastQ=${Trimmed_FastQ%.fq.gz}_trimmed.fq.gz #SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1_trimmed.fq.gz 

		if [ -s $FastQ_file ]; then
			FastQ_array+=($FastQ_file) #push it
		else
			echo -e "\e[031m$FastQ_file not found\e[0m\n"
			exit
		fi
        # trimmed FASTQ
		#if [ -s $Trimmed_FastQ ]; then
		#	Trimmed_FastQ_array+=($Trimmed_FastQ) #push it
		#else
		#	echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
		#	exit
		#fi
	done # end of per FastQ_file 

	###########
	# SALMON 
	###########
	SALMON_OUT=$PROJECT_DIR/Salmon/$TR_PREFIX.$ENS_VER/$Barcode
	if [ $RUN_SALMON -eq 1 ]; then
	    mkdir -p $SALMON_OUT
		Merged_FastQ=$(printf " %s" "${FastQ_array[@]}")
		Merged_FastQ=${Merged_FastQ:1} # remove first space 

		echo -e "salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -r $Merged_FastQ -g $GTF -o $SALMON_OUT"
		time salmon quant \
			-p $NT \
			-i $SALMON_INDEX \
			-l $SALMON_LIBTYPE \
			-r $Merged_FastQ \
			-g $GTF \
			-o $SALMON_OUT
	fi # end of RUN_SALMON
else
	echo -e "\e[031m Paired-end mode\e[0m\n"
	Fw_FastQ_array=() #initialise
	Rv_FastQ_array=() #initialise

	Trimmed_Fw_FastQ_array=() #initialise
	Trimmed_Rv_FastQ_array=() #initialise

	FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*r_1.fq.gz | grep -v lost`
	Cell_lane_array=`for i in $FASTQ_FILES; do echo $i | cut -d/ -f7 | cut -d. -f3,4 ; done`
	## For each lane 
	for CELL_LANE in ${Cell_lane_array[*]} 
	do
		Cell=`echo $CELL_LANE | cut -d. -f1`
		Lane=`echo $CELL_LANE | cut -d. -f2`

		# Un-trimmed FastQ
		Fw_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_1.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_1.fq.gz
		Rv_FastQ=$HOME/data/fastq/$SLX/$SLX.$Barcode.$Cell.$Lane.r_2.fq.gz # SLX-8074.A001.C2N01ACXX.s_1.r_2.fq.gz
		if [ ! -s $Fw_FastQ ];then
			echo -e "\e[031m$Fw_FastQ not found\e[0m\n"
			exit
		else
			Fw_FastQ_array+=($Fw_FastQ) #push it
		fi

		if [ ! -s $Rv_FastQ ];then
			echo -e "\e[031m$Rv_FastQ not found\e[0m\n"
			exit
		else
			Rv_FastQ_array+=($Rv_FastQ) #push it
		fi

		# Trimmed FastQ
		Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1_val_1.fq.gz 
		Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2_val_2.fq.gz 
		#if [ -s $Trimmed_Fw_FastQ ];then
		#	Trimmed_Fw_FastQ_array+=($Trimmed_Fw_FastQ) #push it
		#else
		#	echo -e "\e[031m$Trimmed_Fw_FastQ not found\e[0m\n"
		#	exit
		#fi
        #
		#if [ -s $Trimmed_Rv_FastQ ];then
		#	Trimmed_Rv_FastQ_array+=($Trimmed_Rv_FastQ) #push it
		#else
		#	echo -e "\e[031m$Trimmed_Rv_FastQ not found\e[0m\n"
		#	exit
		#fi
	done # endof CELL_LANE

	############
	# Salmon 
	############
	#SALMON_OUT=$PROJECT_DIR/Salmon/$TR_PREFIX.$ENS_VER.novel.10/$Barcode
	SALMON_OUT=$PROJECT_DIR/Salmon/$TR_PREFIX.$ENS_VER/$Barcode
	if [ $RUN_SALMON -eq 1 ]; then
	    mkdir -p $SALMON_OUT
		Merged_Fw_FastQ=$(printf " %s" "${Fw_FastQ_array[@]}")
		Merged_Fw_FastQ=${Merged_Fw_FastQ:1} # remove first space 

		Merged_Rv_FastQ=$(printf " %s" "${Rv_FastQ_array[@]}")
		Merged_Rv_FastQ=${Merged_Rv_FastQ:1} # remove first space 

		echo -e "salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -1 $Merged_Fw_FastQ -2 $Merged_Rv_FastQ -g $GTF -o $SALMON_OUT"
		time salmon quant \
			-p $NT \
			-i $SALMON_INDEX \
			-l $SALMON_LIBTYPE \
			-1 $Merged_Fw_FastQ -2 $Merged_Rv_FastQ \
			-g $GTF \
			-o $SALMON_OUT
	fi # end of RUN_SALMON
fi


if [ -s $SALMON_OUT/quant.sf ]; then
    awk '/^ENST/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.sf | sort > $SALMON_OUT/$SLX.$Barcode.quant.$TR_PREFIX.$ENS_VER.enst.clean.txt
    #awk '/TCONS/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.sf | sort > $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.tcons.clean.txt
fi
if [ -s $SALMON_OUT/quant.genes.sf ]; then
    awk '/^ENSG/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.genes.sf | sort > $SALMON_OUT/$SLX.$Barcode.quant.$TR_PREFIX.$ENS_VER.ensg.clean.txt
    #awk '/XLOC/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.genes.sf | sort > $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.xloc.clean.txt 
fi
#join <(sort -k1,1 $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.tcons.clean.txt) <(sort -k1,1 /home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.hgnc_symbol.txt)  | sort -k 3,3 | awk 'BEGIN{OFS="\t"}{gene_count[$3]+=$2}END{for (gene in gene_count) print gene,gene_count[gene]}' | sort -k1,1 > $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.hgnc.clean.txt 
#join <(sort -k1,1 $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.tcons.clean.txt) <(sort -k1,1 /home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.ensg.txt)  | sort -k 3,3 | awk 'BEGIN{OFS="\t"}{gene_count[$3]+=$2}END{for (gene in gene_count) print gene,gene_count[gene]}' | sort -k1,1 > $SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.ensg.clean.txt 

echo -e "\e[32m$Barcode done for salmon\e[0m" 
