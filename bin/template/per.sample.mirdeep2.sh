#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 27/Aug/2015
# Last modified: 27/Aug/2015
# Optimised and customised to run at the Darwin HPC
# based on miRDeep2

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.Homo_sapiens.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  
Lane="MY_LANE" # e.g. s_1 
Chunk="MY_CHUNK" # e.g. 1 

mkdir_unless $PROJECT_DIR	
mkdir_unless $PROJECT_DIR/miRDeep2

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
	Trimmed_FastQ_array=() #initialise
	FASTQ_FILES=`ls $HOME/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`
	for FastQ_file in $FASTQ_FILES	# per fastq file (per lane)
	do
		##################################################
		# 1. Run the inital fastqc 
		##################################################
		if [ $RUN_FASTQC -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/FastQC
			mkdir_unless $PROJECT_DIR/FastQC/$Barcode
			echo "fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode"
			time fastqc $FastQ_file -o $PROJECT_DIR/FastQC/$Barcode 
		fi
		#/home/ssg29/data/fastq/SLX-9176/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1_trimmed.fq.gz
		#/home/ssg29/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01.v1/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1.fq.gz_trimming_report.txt
		Trimmed_FastQ=${FastQ_file/data\/fastq\/$SLX/results\/$PROJECT\/Trim\/$Barcode}
		Trimmed_FastQ=${Trimmed_FastQ%.fq.gz}_trimmed.fq.gz #SLX-8546.D709_D501.C4FN0ANXX.s_6.r_1_trimmed.fq.gz 
		############################################################
		# 2. Trimming via cutadapt 
		############################################################
			#mkdir_unless $PROJECT_DIR/Trim
			#mkdir_unless $PROJECT_DIR/Trim/$Barcode
			#echo -e "\ncutadapt $FastQ_file"
			#time cutadapt \
			#	--trimmed-only \
			#	-q 20 \
			#	-O 8 \
			#	-m 15 \
			#	-a NEB3PrimeAdaptor=AGATCGGAAGAGCACACGTCT \
			#	-g NEB5PrimeAdaptor=GTTCAGAGTTCTACAGTCCGACGATC \
			#	-o $Trimmed_FastQ \
			#	$FastQ_file &> ${Trimmed_FastQ%_trimmed.fq.gz}.fq.gz_trimming_report.txt
			#if [ ! -s $Trimmed_FastQ ];then
			#	echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
			#	exit
			#else
			#	echo -e "\nfastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode"
			#	#time fastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode 
			#fi
		if [ $RUN_TRIM_NEB -eq 1 ]; then
			mkdir_unless $PROJECT_DIR/Trim
			mkdir_unless $PROJECT_DIR/Trim/$Barcode
			echo -e "\ncutadapt $FastQ_file"
			time cutadapt \
				--trimmed-only \
				-q $TR_QUAL \
				-O $TR_STRNCY \
				-m $TR_LEN \
				-a $TR_ADAPTOR1 \
				-g $TR_ADAPTOR2 \
				-o $Trimmed_FastQ \
				$FastQ_file &> ${Trimmed_FastQ%_trimmed.fq.gz}.fq.gz_trimming_report.txt
			if [ ! -s $Trimmed_FastQ ];then
				echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
				exit
			else
				#echo -e "\nfastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode"
				#time fastqc $Trimmed_FastQ -o $PROJECT_DIR/Trim/$Barcode 
				echo -e "\nfastqc skipped"
			fi
		fi # end of run trim

		if [ -s $Trimmed_FastQ ]; then
			Trimmed_FastQ_array+=($Trimmed_FastQ) #push it multi-laned fastq
		else
			echo -e "\e[031m$Trimmed_FastQ not found\e[0m\n"
			exit
		fi
	done # end of per FastQ_file 

	My_fastq=$PROJECT_DIR/Trim/$Barcode/$SLX.$Barcode.trim.merged.fq.gz
	My_mapped=$PROJECT_DIR/miRDeep2/mapper/$Barcode/$SLX.$Barcode.mirdeep2.mapped.arf
	if [ $RUN_MIRDEEP_MAPPER -eq 1 ];then
		#####################
		## 3. Merge FastQ ##
		#####################
		Merged_FastQ=$(printf " %s" "${Trimmed_FastQ_array[@]}")
		if [ ${#Trimmed_FastQ_array[@]} -eq 1 ]; then
			echo -e "ln -s $Merged_FastQ $My_fastq\n"
			time ln -s $Merged_FastQ $My_fastq
		else
		# only merge physically if there are at least two fastq files
			echo -e "zcat $Merged_FastQ | gzip -9 > $My_fastq\n"
			time zcat $Merged_FastQ | gzip -9 > $My_fastq
		fi
		#######################
		## 4. fastq to fasta ##
		#######################
		if [ ! -s ${My_fastq%.fq.gz}.fasta ];then
			#time fastqutils tofasta $My_fastq > ${My_fastq%.fq.gz}.fasta # this fails due to the space in header
			echo -e "zcat $My_fastq | awk \n"
			time zcat $My_fastq | awk "NR%4==1{print \">\"\$1}NR%4==2{print \$0}" > ${My_fastq%.fq.gz}.fasta
			#zcat ~/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1_trimmed.fq.gz | awk 'NR%4==1{print ">"$1}NR%4==2{print $0}' > ~/results/SLX-9176.Homo_sapiens.v1/Trim/NEBsmRNA01/SLX-9176.NEBsmRNA01.C7TWVANXX.s_6.r_1_trimmed.fasta
		fi

		#########################
		## 5. miRDeep2 Mapper  ##
		#########################
		if [ -s ${My_fastq%.fq.gz}.fasta ];then
			mkdir_unless $PROJECT_DIR/miRDeep2/mapper
			mkdir_unless $PROJECT_DIR/miRDeep2/mapper/$Barcode
			cd $PROJECT_DIR/miRDeep2/mapper/$Barcode
			echo -e "mapper.pl ${My_fastq%.fq.gz}.fasta\n"
			# -n: overwite
			# -c: input file is fasta format
			# -m: collapse reads
			# -j: remove all entries that have a sequence that contains letters
			#     other than a,c,g,t,u,n,A,C,G,T,U,N
			time mapper.pl \
				${My_fastq%.fq.gz}.fasta \
				-n -c -m -j \
				-o $NT \
				-p $BOWTIE_INDEX_BASE \
				-s ${My_fastq%.fq.gz}.collapsed.fasta \
				-t $My_mapped \
				-v &> ${My_mapped%.arf}.log
		else
			echo -e "\e[031m${My_fastq%.fq.gz}.fasta not found\e[0m\n"
			exit
		fi
	fi

	######################################################
	## 6. miRDeep2 Quantifier (reads must be collapsed) ##
	######################################################
	if [ $RUN_MIRDEEP_QUANTIFIER -eq 1 ];then
		if [ ! -s ${My_fastq%.fq.gz}.collapsed.fasta ];then
			# collapse reads via mapper.pl (collapsing only, not actual mapping)
			echo -e "mapper.pl ${My_fastq%.fq.gz}.fasta -c -m -s ${My_fastq%.fq.gz}.collapsed.fasta\n"
			time mapper.pl ${My_fastq%.fq.gz}.fasta -c -m -s ${My_fastq%.fq.gz}.collapsed.fasta
		fi

		if [ -s ${My_fastq%.fq.gz}.collapsed.fasta ];then
			mkdir_unless $PROJECT_DIR/miRDeep2/quantifier
			mkdir_unless $PROJECT_DIR/miRDeep2/quantifier/$Barcode
			cd $PROJECT_DIR/miRDeep2/quantifier/$Barcode
			time quantifier.pl \
				-t $MIRBASE_SPECIES \
				-r ${My_fastq%.fq.gz}.collapsed.fasta \
				-p $MIRBASE_PRECURSOR \
				-m $MIRBASE_MATURE \
				-k
		else
			echo -e "\e[031m${My_fastq%.fq.gz}.collapsed.fasta not found\e[0m\n"
			exit
		fi
	fi

	######################
	## 7. Core miRDeep2 ##
	######################
	if [ $RUN_MIRDEEP_CORE -eq 1 ];then
		if [ -s ${My_fastq%.fq.gz}.collapsed.fasta ] && [ -s $My_mapped ];then
			mkdir_unless $PROJECT_DIR/miRDeep2/core
			mkdir_unless $PROJECT_DIR/miRDeep2/core/$Barcode
			cd $PROJECT_DIR/miRDeep2/core/$Barcode
			echo -e "miRDeep2.pl ${My_fastq%.fq.gz}.collapsed.fasta \n"
			time miRDeep2.pl \
				${My_fastq%.fq.gz}.collapsed.fasta \
				$GENOME \
				$My_mapped \
				$MIRBASE_MATURE \
				${MIRBASE_MATURE%.fa}.homologues.fa \
				$MIRBASE_PRECURSOR \
				-t $MIRDEEP_SPECIES \
				-r $Barcode \
				-z $Barcode \
				-v &> $PROJECT_DIR/miRDeep2/core/$Barcode/$SLX.$Barcode.mirdeep2.core.log
		else
			echo -e "\e[031m${My_fastq%.fq.gz}.collapsed.fasta or $My_mapped not found\e[0m\n"
			exit
		fi

	fi

	#######################
	# 8. miRNA Enrichment # 
	#######################
	# miRNAs_expressed_all_samples_02_09_2015_t_23_08_11NEBsmRNA02.csv
	# merge count of 5p and 3p by the same mature-miRNA (pre-cursor)
	# to be used by DESeq2
	if [ $RUN_MIRNA_CNT -eq 1 ];then
		MIRNA_EXP=`ls $PROJECT_DIR/miRDeep2/core/$Barcode/miRNAs_expressed_all_samples_*$Barcode.csv`
		if [ -s ${MIRNA_EXP[0]} ];then
			awk 'BEGIN{OFS="\t"}!/^#/{Precursor[$3]+=$2}END{for(p in Precursor)print p,Precursor[p]}' ${MIRNA_EXP[0]} > $PROJECT_DIR/miRDeep2/core/$Barcode/$PROJECT.$Barcode.miRNA.cnt.txt 
			awk 'BEGIN{OFS="\t"}!/^#/{Mature[$1]+=$2}END{for(m in Mature)print m,Mature[m]}' ${MIRNA_EXP[0]} > $PROJECT_DIR/miRDeep2/core/$Barcode/$PROJECT.$Barcode.miRNA.mature.cnt.txt 
			#awk 'BEGIN{OFS="\t"}!/^#/{split($2,b,"."); print $1,b[1]}' ${MIRNA_EXP[0]} > $PROJECT_DIR/miRDeep2/core/$Barcode/$PROJECT.$Barcode.miRNA.mature.cnt.txt
		else
			echo -e "\e[031m${MIRNA_EXP[0]} not found\e[0m\n"
			exit
		fi
	fi

	#######################
	# 9. piRNA Enrichment # 
	#######################
	# to be used by DESeq2
	if [ $RUN_PIRNA_CNT -eq 1 ];then
		mkdir_unless $PROJECT_DIR/piRBase
		mkdir_unless $PROJECT_DIR/piRBase/$Barcode
		if [ -s $My_mapped ];then
			echo -e "parse_mappings.pl $My_mapped\n"
			time parse_mappings.pl \
				$My_mapped \
				-a 0 \
				-i 3 \
				| awk 'BEGIN{OFS="\t"}{start=$8-1; split($1,seq,"_x"); print "chr"$6,start,$9,$1,seq[2],$11}' \
				| intersectBed -a $PIRBASE_BED -b stdin -s -wao -f 0.3 \
				| sort -k1,1 -k2,2 -k3,3 -k4,4 -k6,6 \
				| groupBy -i stdin -grp 1,2,3,4,6 -c 11 -o sum \
				| awk 'BEGIN{OFS="\t"}{if($6<1)cnt=0;else cnt=$6; print $1,$2,$3,$4,cnt,$5}' \
				| sort -k4,4V > $PROJECT_DIR/piRBase/$Barcode/$PROJECT.$Barcode.piRBase.sorted.bed 
		else
			echo -e "\e[031m$My_mapped not found\e[0m\n"
			exit
		fi

		if [ -s $PROJECT_DIR/piRBase/$Barcode/$PROJECT.$Barcode.piRBase.sorted.bed ];then
			echo -e "groupBy -i $PROJECT_DIR/piRBase/$Barcode/$PROJECT.$Barcode.piRBase.sorted.bed\n"
			time groupBy \
				-i $PROJECT_DIR/piRBase/$Barcode/$PROJECT.$Barcode.piRBase.sorted.bed \
				-g 4 \
				-c 5 \
				-o sum > $PROJECT_DIR/piRBase/$Barcode/$PROJECT.$Barcode.piRBase.read.cnt.txt
		else
			echo -e "\e[031m$PROJECT_DIR/piRBase/$Barcode/$PROJECT.$Barcode.piRBase.sorted.bed not found\e[0m\n"
			exit
		fi
	fi
else
	echo -e "\e[031This script is for sigle-end read only\e[0m\n"
	exit
fi

echo -e "\e[32m$Barcode done for mapping & counting\e[0m" 
