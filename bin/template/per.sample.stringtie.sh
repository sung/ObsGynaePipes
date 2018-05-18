#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 21/Oct/2015
# Last modified: 21/Oct/2015
# Optimised and customised to run at the Darwin HPC
# based on StringTie <https://github.com/gpertea/stringtie>

source $HOME/Pipelines/config/rna_seq.config # export envrionment variables
source $HOME/Pipelines/lib/sung.sh #defines 'mkdir_unless'
source $HOME/config/sung.hpc.bash #defines PATH 

SLX="MY_SLX" # e.g. SLX-8080 
PROJECT="MY_SLX"."MY_VERSION" # e.g. SLX-8080.Homo_sapiens.v1
PROJECT_DIR=$RESULT_DIR/$PROJECT
Barcode="MY_BARCODE" # e.g. A001
Cell="MY_CELL" # e.g. C48CWACXX  

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

########
## FG ##
########
#for j in $BARCODE; do echo -e "$SLX $j\n"; samtools merge -f ~/results/$SLX.Homo_sapiens.SE125.v1/TopHat/$j/final_merged_accepted_hits.bam ~/results/$SLX.Homo_sapiens.SE125.v1/TopHat/$j/merged_accepted_hits.bam ~/results/$SLX.Homo_sapiens.v4/TopHat/$j/merged_accepted_hits.bam; samtools index ~/results/$SLX.Homo_sapiens.SE125.v1/TopHat/$j/final_merged_accepted_hits.bam; done

# assuming single-end
if [ $IS_SE -eq 1 ]; then
	My_bam=$PROJECT_DIR/TopHat/$Barcode/merged_accepted_hits.bam # SE125 (JD|FG)
	#My_bam=$PROJECT_DIR/TopHat/$Barcode/final_merged_accepted_hits.bam # SE125 + SE50 (FG)
	if [ ! -s $My_bam ];then 
		echo -e "\e[031m$My_bam not found. tophat failed? \e[0m\n"
		exit
	fi
	######################
	# 1. CoverageBed
	# Input: bam from tophat
	# (NB, to use 'bedtools genomecov', BAM _must_ be sorted by position, 
	# which is already done by tophat)
	# Output: $sample.genomecov.txt
	######################
	if [ $RUN_GENOMECOV_UNMAPPED -eq 1 ]; then
		mkdir_unless $PROJECT_DIR/Coverage/$Barcode
		# bam2coverage
		# echo -e "bedtools genomecov $My_bam\n"
		# time bedtools genomecov -split -ibam $My_bam > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.txt
		# bam2bedgraph
		echo -e "bedtools genomecov -bg -split -ibam $My_bam\n"
		time bedtools genomecov -bg -split -ibam $My_bam > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bedgraph
		# make chrom.sizes
		echo -e "samtools idxstats $My_bam\n"
		time samtools idxstats $My_bam | awk 'BEGIN{OFS="\t"}!/*/{print $1,$2}' > $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.samtool.idxstats.chr.size.txt
		# bedgraph2bw
		echo -e "bedGraphToBigWig $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bedgraph\n"
		time bedGraphToBigWig $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bedgraph $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.samtool.idxstats.chr.size.txt $PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.genomecov.split.unmap.rescued.bw
	fi
	###############
	## StringTie ## 
	## v1.2.3    ##
	###############
	# -f: minimum isoform fraction (default: 0.1)
	# -m: minimum assembled transcript length (default: 200)
	# -c: minimum reads per bp coverage to consider for transcript assembly (default: 2.5)
	# -j: minimum junction coverage (default: 1)
	# -M: fraction of bundle allowed to be covered by multi-hit reads (default:0.95)
	# -B: enable output of Ballgown table files which will be created in the
	#     same directory as the output GTF (requires -G, -o recommended)
	# -C: output a file with reference transcripts that are covered by reads
	# -A: gene abundance estimation output file
	if [ $RUN_STRINGTIE_ASM -eq 1 ];then 
		mkdir_unless $PROJECT_DIR/StringTie/$Barcode
		echo -e "stringtie $My_bam\n"
		time stringtie \
			$My_bam \
			-p $NT \
			-G $GTF \
			-l $SLX.$Barcode.STRG \
			-c $STR_MIN_READ \
			-j $STR_MIN_JNC_COV \
			-B \
			-o $PROJECT_DIR/StringTie/$Barcode/$PROJECT.$Barcode.stringtie.$TR_PREFIX.$ENS_VER.assembly.gtf \
			-C $PROJECT_DIR/StringTie/$Barcode/$PROJECT.$Barcode.stringtie.$TR_PREFIX.$ENS_VER.coverage.txt \
			-A $PROJECT_DIR/StringTie/$Barcode/$PROJECT.$Barcode.stringtie.$TR_PREFIX.$ENS_VER.gene.cnt.txt
	fi
else
	echo -e "\e[031This script is for signle-end only\e[0m\n"
	exit
fi

echo -e "\e[32m$Barcode done for StringTie & Bam2BigWig\e[0m" 
