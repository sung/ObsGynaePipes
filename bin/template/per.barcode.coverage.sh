#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 30/Jan/2015
# Last modified: 30/Jan/2015 
# Optimised and customised to run at the SBS machine

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/lib/sung.sh

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

# SLX-7634.HALO1.000000000-A5NK2.bwa.bam
My_bam=$PROJECT_DIR/BWA/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.bam 
if [ ! -s $My_bam ];then
	echo -e "\e[031m$My_bam is expected but not found\e[0m\n"
	exit
fi

My_coverageBed=$PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.bam
if [ $RUN_COVBED -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/Coverage
	mkdir_unless $PROJECT_DIR/Coverage/$Barcode
	### CoverageStas [ PCR/optical duplicateReads were removed]
	echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET > $My_coverageBed.CoverageStats \n"
	time samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET > $My_coverageBed.CoverageStats
	### CoveragePerBase [ PCR/optical duplicateReads were removed]
	echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET -d > $My_coverageBed.CoveragePerBase \n"
	time samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET -d > $My_coverageBed.CoveragePerBase
	### CoverageHistogram [ PCR/optical duplicateReads were removed]
	echo -e "samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET -hist > $My_coverageBed.CoverageHist \n"
	time samtools view -uF 0x400 $My_bam | coverageBed -abam stdin -b $TARGET -hist > $My_coverageBed.CoverageHist

	echo -e "\e[032mcoverageBed done for $My_bam\e[0m\n"
fi


My_CallableLoci_Out=$PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.bam.bases.callable
My_CallableLoci_Summary=$PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.bam.bases.callable.summary
if [ $RUN_GATK_CALLABLE -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/Coverage
	mkdir_unless $PROJECT_DIR/Coverage/$Barcode
	###  No of callable Bases
	echo -e "java -Xmx$GATK_MEM -jar $GATK/GenomeAnalysisTK.jar -T CallableLoci\n"
	time java -Xmx$GATK_MEM -jar $GATK/GenomeAnalysisTK.jar \
		-T CallableLoci \
		-R $GENOME \
		-I $My_bam \
		-L $TARGET \
		-l INFO \
		--minMappingQuality $GATK_MIN_MAPPING_QUAL \
		--minBaseQuality $GATK_MIN_BASE_QUAL \
		-format BED \
		-summary $My_CallableLoci_Summary \
		-o $My_CallableLoci_Out
	echo -e "\e[032mGATK CallableLoci done for $My_bam\e[0m\n"

	if [ -s $My_CallableLoci_Out ];then
		### Generate Bases callable by target file
		echo -e "awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $My_CallableLoci_Out | intersectBed -a $TARGET -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $My_CallableLoci_Out.byTarget \n"
		time awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $My_CallableLoci_Out | intersectBed -a $TARGET -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $My_CallableLoci_Out.byTarget
	fi

	if [ -s $My_CallableLoci_Out.byTarget ];then
		#### Perl script to generate this format
		##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
		##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
		echo -e "perl $REARRANGE_COV $My_CallableLoci_Out.byTarget > $My_CallableLoci_Out.byTarget.rearranged \n"
		time perl $REARRANGE_COV $My_CallableLoci_Out.byTarget > $My_CallableLoci_Out.byTarget.rearranged
	fi
fi

My_DOV_Out=$PROJECT_DIR/Coverage/$Barcode/$SLX.$Barcode.$Cell.bwa.sorted.realigned.recalibrated.OnTarget.q$MAPQ.bam.depth.of.cov.txt
if [ $RUN_GATK_DOV -eq 1 ]; then
	mkdir_unless $PROJECT_DIR/Coverage
	mkdir_unless $PROJECT_DIR/Coverage/$Barcode
	### Depth of coverage
	echo -e "java -Xmx$GATK_MEM -jar $GATK/GenomeAnalysisTK.jar -T DepthOfCoverage\n"
	time java -Xmx$GATK_MEM -jar $GATK/GenomeAnalysisTK.jar \
		-T DepthOfCoverage \
		-R $GENOME \
		-I $My_bam \
		-L $TARGET \
		--minMappingQuality $GATK_MIN_MAPPING_QUAL \
		--minBaseQuality $GATK_MIN_BASE_QUAL \
		--includeDeletions \
		--summaryCoverageThreshold \
		--outputFormat table \
		--summaryCoverageThreshold \
		--omitDepthOutputAtEachBase \
		-o $My_DOV_Out
	echo -e "\e[032mGATK DepthOfCoverage done for $My_bam\e[0m\n"
fi

if [ $RUN_PICARD_SUMMARY -eq 1 ];then
	### CollectAlignmentSummaryMetrics
	echo -e "java -Xmx$GATK_MEM -jar $PICARD_HOME/CollectAlignmentSummaryMetrics.jar\n"
	time java -Xmx$GATK_MEM -jar $PICARD_HOME/CollectAlignmentSummaryMetrics.jar \
		INPUT=$My_bam \
		OUTPUT=$My_bam.picard.aln.summary.txt \
		REFERENCE_SEQUENCE=$GENOME \
		ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT
	echo -e "\e[032mPicard CollectAlignmentSummaryMetrics done for $My_bam\e[0m\n"
fi

echo -e "\e[32m$Barcode done for per.barcode.coverage\e[0m" 
