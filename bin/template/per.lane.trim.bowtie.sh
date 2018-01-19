#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 26/Feb/2015
# Last modified: 27/Feb/2015 
# Optimised and customised to run at the SBS and HPC machine

source $HOME/Pipelines/config/chip_seq.config 
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless', 'make_run_script'

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT_DIR=$RESULT_DIR/MY_SLX.MY_VERSION # e.g. SLX-8075.v1
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 
FASTQ_DIR=$HOME/data/fastq/$SLX

mkdir_unless $PROJECT_DIR	

echo -e `hostname`
echo -e "NO. of thread=$NT"
printf "\n\e[32mSample.Lane=$SLX.$Barcode.$Lane\e[0m\n"

#/disk1/ssg29/data/fastq/SLX-9338/SLX-9338.HAL21.000000000-ACU88.s_1.r_1.fq.gz
#/disk1/ssg29/data/fastq/SLX-9338/SLX-9338.HAL21.000000000-ACU88.s_1.r_2.fq.gz
FASTQ_FILES=`ls $FASTQ_DIR/$SLX.$Barcode.$Cell.$Lane.*.fq.gz | grep -v lost`

##################################################
# 2. run initial fastqc 
# output: results/SLX123/FastQC/A001/
##################################################
if [ $RUN_FASTQC -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/FastQC
	mkdir_unless $PROJECT_DIR/FastQC/$Barcode
	for FastQ_file in $FASTQ_FILES
	do
	echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
		time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
	done
fi

############################################################
# 3. Trim FastQ file per lane 	
# input: .fastq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1_val_1.fq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1_val_2.fq.gz
# output: results/SGA_v6/Trim/A001/SLX-8074.A001.C2N01ACXX.s_1.r_1.1.fastq.gz_trimming_report.txt
############################################################
if [ $RUN_TRIM -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/Trim
	mkdir_unless $PROJECT_DIR/Trim/$Barcode
	#for FastQ_file in $FASTQ_FILES
	#do
	#	FastQ_array+=($FastQ_file) #push it
	#done
	#dummy_fastq=$(printf " %s" "${FastQ_array[@]}")
	dummy_fastq=$(printf " %s" "${FASTQ_FILES[@]}")

	if [ $IS_PE -eq 1 ]; then
		echo -e "time trim_galore --paired --output_dir $PROJECT_DIR/Trim/$Barcode --fastqc $dummy_fastq\n"
		time trim_galore \
			--paired \
			--length $TR_LEN \
			--fastqc \
			--adapter $TR_ADAPTOR1 \
			--output_dir $PROJECT_DIR/Trim/$Barcode \
			$dummy_fastq
	else
		echo -e "time trim_galore --output_dir $PROJECT_DIR/Trim/$Barcode --fastqc $dummy_fastq\n"
		time trim_galore \
			--length $TR_LEN \
			--fastqc \
			--adapter $TR_ADAPTOR1 \
			--output_dir $PROJECT_DIR/Trim/$Barcode \
			$dummy_fastq
	fi
fi # end of run trim

if [ $IS_PE -eq 1 ]; then
	Trimmed_Fw_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1_val_1.fq.gz 
	Trimmed_Rv_FastQ=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_2_val_2.fq.gz 
	if [ ! -s $Trimmed_Fw_FastQ ];then
		echo -e "\e[031m$Trimmed_Fw_FastQ not found\e[0m\n"
		exit
	fi
	if [ ! -s $Trimmed_Rv_FastQ ];then
		echo -e "\e[031m$Trimmed_Rv_FastQ not found\e[0m\n"
		exit
	fi
	TRIMMED_FASTQ_FILES="$Trimmed_Fw_FastQ $Trimmed_Rv_FastQ"
else
	#SLX-8764.A001.C4EYHANXX.s_7.r_1_trimmed.fq.gz
	TRIMMED_FASTQ_FILES="$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.$Cell.$Lane.r_1"_trimmed.fq.gz
	if [ ! -s $TRIMMED_FASTQ_FILES ];then
		echo -e "\e[031m$TRIMMED_FASTQ_FILES not found\e[0m\n"
		exit
	fi
fi

if [ $RUN_BOWTIE -eq 1 ]; then
	MY_BOWTIE_OUT=$PROJECT_DIR/Bowtie/$Barcode/$SLX.$Barcode.$Cell.$Lane.bowtie.sam
	mkdir_unless $PROJECT_DIR/Bowtie
	mkdir_unless $PROJECT_DIR/Bowtie/$Barcode

	echo -e "bowtie2 --threads $NT -U $TRIMMED_FASTQ_FILES -x $BOWTIE2_INDEX_BASE -S $MY_BOWTIE_OUT 2> ${MY_BOWTIE_OUT%.sam}.log"
	time bowtie2 \
		--threads $NT \
		-U $TRIMMED_FASTQ_FILES \
		-x $BOWTIE2_INDEX_BASE \
		--rg-id $SLX \
		--rg PL:illumina \
		--rg PU:run	\
		--rg LB:$SLX \
		--rg SM:$Barcode \
		--rg CN:CamObsGynae \
		-S $MY_BOWTIE_OUT 2> ${MY_BOWTIE_OUT%.sam}.log

	if [ -s $MY_BOWTIE_OUT ];then
		echo -e "samtools view -bS $MY_BOWTIE_OUT  > ${MY_BOWTIE_OUT%.sam}.bam"
		time samtools view -bS $MY_BOWTIE_OUT  > ${MY_BOWTIE_OUT%.sam}.bam 

		if [ -s ${MY_BOWTIE_OUT%.sam}.bam ];then
			echo -e "rm -f $MY_BOWTIE_OUT \n"
			time rm -f $MY_BOWTIE_OUT 

			echo -e "samtools sort -m 25G ${MY_BOWTIE_OUT%.sam}.bam ${MY_BOWTIE_OUT%.sam}.sorted"
			time samtools sort -m 25G ${MY_BOWTIE_OUT%.sam}.bam ${MY_BOWTIE_OUT%.sam}.sorted
		else
			echo -e "\e[031m${MY_BOWTIE_OUT%.sam}.bam not found\e[0m\n"
			exit
		fi
		if [ -s ${MY_BOWTIE_OUT%.sam}.sorted.bam ];then
			echo -e "rm ${MY_BOWTIE_OUT%.sam}.bam"
			time rm ${MY_BOWTIE_OUT%.sam}.bam 

			echo -e "samtools index ${MY_BOWTIE_OUT%.sam}.sorted.bam"
			time samtools index ${MY_BOWTIE_OUT%.sam}.sorted.bam
		fi
	else
		echo -e "\e[031m$MY_BOWTIE_OUT not found\e[0m\n"
		exit
	fi
fi

# ~/results/SLX-8766.Homo_sapiens.v1/Bowtie/A005/SLX-8766.A005.C5UP5ANXX.s_2.bowtie.sorted.bam
My_bam=$PROJECT_DIR/Bowtie/$Barcode/$SLX.$Barcode.$Cell.$Lane.bowtie.sorted.bam
if [ $RUN_GENOMECOV -eq 1 ];then # run bedtools genomecov
	if [ -s $My_bam ];then
		echo -e "bedtools genomecov -ibam $My_bam > ${My_bam%.bam}.bedtools.genomecov.txt"
		time bedtools genomecov -ibam $My_bam > ${My_bam%.bam}.bedtools.genomecov.txt
	else
		echo -e "\e[031m$My_bam not found\e[0m\n"
		exit
	fi
fi

echo -e "\e[32m$Barcode done for perl.lane.trim.bowtie\e[0m" 
