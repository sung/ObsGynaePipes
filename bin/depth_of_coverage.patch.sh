#!/bin/bash

#Barcodes=(A001 A002 A003 A004 A005 A010)
Barcodes=(A012)
Depth=(1 3 5 7 10)

My_bam=/whale-data/ssg29/Methyl-Seq/results/SGA_v3/Bismark/Alignment/A003/A003.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.sorted.bam
genomeSize=`samtools view -H $My_bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'`

for Barcode in ${Barcodes[*]}
do
	echo -e "$Barcode"	
	CpG_report=/whale-data/ssg29/Methyl-Seq/results/SGA_v3/Bismark/MethylCall/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.q1.deduplicated.CpG_report.txt
	#time cat $CpG_report | awk '{sum+=($4+$5)}END{print "CpG Depth of Cov=",sum/NR}'
	Genome_Cov=/whale-data/ssg29/Methyl-Seq/results/SGA_v3/Coverage/$Barcode/$Barcode.bedtools.genomecov.txt

	#time awk 'BEGIN{print "chrBase chr base strand coverage freqC freqT"}{str="F"}{if($3=="-"){str="R"}}{if($4+$5>0){cov=$4+$5;print $1"."$2,$1,$2,str,cov,$4/cov*100,$5/cov*100}}' $CpG_report > ${CpG_report%.txt}.methylkit.txt

	time cat $CpG_report | awk '{sum+=$4+$5}END{print "CpG Depth of Cov=",sum/NR}' > /whale-data/ssg29/Methyl-Seq/results/SGA_v3/Coverage/$Barcode/$Barcode.CpG.coverage.txt
	time grep -P '^genome' $Genome_Cov | awk "{sum+=\$2*\$3} END { print \"depth of coverage= \",sum/$genomeSize}" > /whale-data/ssg29/Methyl-Seq/results/SGA_v3/Coverage/$Barcode/$Barcode.coverage.txt

	for d in ${Depth[*]}
	do
		echo "calculating x$d coverage..."

		# CpG Coverage
		time cat $CpG_report | awk "{if(\$4+\$5>=$d){sum++}}END{print \"x$d CpG BoC=\",sum/NR*100}" >> /whale-data/ssg29/Methyl-Seq/results/SGA_v3/Coverage/$Barcode/$Barcode.CpG.coverage.txt

		#Genome Coverage
		time grep -P '^genome' $Genome_Cov | awk "{if(\$2>=$d){sum+=\$3}} END {print \"x$d BoC=\",sum/$genomeSize*100}" >> /whale-data/ssg29/Methyl-Seq/results/SGA_v3/Coverage/$Barcode/$Barcode.coverage.txt
	done

done
