#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 22/Apr/2014
# Last modified: 22/Apr/2014 

# define a function
function mkdir_unless(){
	if [ ! -d $1 ]; then
		mkdir $1
	fi	
}

BowtieIndexBase="/whale-data/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
Genome="$BowtieIndexBase.fa" 
GTF='/whale-data/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf'
#GTF_ASM='/whale-data/ssg29/RNA-Seq/data/Pilot1.gtf_assemblies.txt' # R3-R10
GTF_ASM='/whale-data/ssg29/RNA-Seq/data/Pilot1.R3-R8.gtf_assemblies.txt' # R3-R8
SAMPLE_SHEET='/whale-data/ssg29/RNA-Seq/data/Pilot1.cuffdiff.sample_sheet.txt'
PROJECT='Pilot1_v5'
CONTRAST='/whale-data/ssg29/RNA-Seq/data/Pilot1.contrasts.txt'
CUFFMERGE_DIR="/whale-data/ssg29/RNA-Seq/Results/$PROJECT/Cuffmerge"
CUFFDIFF_DIR="/whale-data/ssg29/RNA-Seq/Results/$PROJECT/Cuffdiff"
LIB_TYPE='fr-firststrand' # see http://cufflinks.cbcb.umd.edu/manual.html#library
NT=8

mkdir_unless /whale-data/ssg29/RNA-Seq/Results/$PROJECT
mkdir_unless $CUFFMERGE_DIR
mkdir_unless $CUFFDIFF_DIR

if [ -s $GTF_ASM ]; then
	# real    13m33.081s
	echo -e "\ncuffmerge -o $CUFFMERGE_DIR -g $GTF -s $Genome -p $NT $GTF_ASM"
	#time cuffmerge -o $CUFFMERGE_DIR -g $GTF -s $Genome -p $NT $GTF_ASM
else
	echo -e "\e[031m$GTF_ASM not found. cufflinks failed? \e[0m\n"
	exit
fi

if [ -s $CUFFMERGE_DIR/merged.gtf ]; then
	# Warning: No conditions are replicated, switching to 'blind' dispersion method

	# 1. contrast & sample sheet: seg fault
	#echo -e "\ntime cuffdiff -o $CUFFDIFF_DIR --library-type $LIB_TYPE -b $Genome -p $NT -u --contrast-file $CONTRAST  --use-sample-sheet $CUFFMERGE_DIR/merged.gtf $SAMPLE_SHEET"
	#time cuffdiff -o $CUFFDIFF_DIR --library-type $LIB_TYPE -b $Genome -p $NT -u --contrast-file $CONTRAST  --use-sample-sheet $CUFFMERGE_DIR/merged.gtf $SAMPLE_SHEET

	# 2. contrast only: seg fault
	#echo -e "\ntime cuffdiff -o $CUFFDIFF_DIR --library-type $LIB_TYPE -b $Genome -p $NT -L R3,R4,R5,R6,R7,R8,R9,R10 --contrast-file $CONTRAST -u $CUFFMERGE_DIR/merged.gtf /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R3/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R4/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R5/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R6/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R7/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R8/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R9/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R10/accepted_hits.bam "
	#time cuffdiff -o $CUFFDIFF_DIR --library-type $LIB_TYPE -b $Genome -p $NT -L R3,R4,R5,R6,R7,R8,R9,R10 --contrast-file $CONTRAST -u $CUFFMERGE_DIR/merged.gtf /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R3/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R4/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R5/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R6/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R7/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R8/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R9/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R10/accepted_hits.bam 

	# 3. no contrast & no sample sheet (cufflink_matched_case_control.log.2)
	# real    1137m36.573s
	echo -e "\ntime cuffdiff -o $CUFFDIFF_DIR --library-type $LIB_TYPE -b $Genome -p $NT -L R3,R4,R5,R6,R7,R8 -u $CUFFMERGE_DIR/merged.gtf /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R3/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R4/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R5/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R6/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R7/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R8/accepted_hits.bam "
	time cuffdiff -o $CUFFDIFF_DIR --library-type $LIB_TYPE -b $Genome -p $NT -L R3,R4,R5,R6,R7,R8 -u $CUFFMERGE_DIR/merged.gtf /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R3/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R4/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R5/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R6/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R7/accepted_hits.bam /whale-data/ssg29/RNA-Seq/Results/Pilot1_v4/TopHat/R8/accepted_hits.bam 
else
	echo -e "\e[031m$CUFFMERGE_DIR/merged.gtf not found. cuffmerge failed? \e[0m\n"
	exit
fi
