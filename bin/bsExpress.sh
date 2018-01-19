#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# Reference: http://code.google.com/p/oxbs-sequencing-qc/wiki/bsExpressDoc
# Called by bin/template/trim_bismark_template.v{n}.sh
# First created: 3/Jul/2014
# Last modified: 3/Jul/2014 

		mkdir_unless $PROJECT_DIR/BsExpress/$Barcode
		# mapping against known control 
		echo -e "bismark $Trimmed_Fw_FastQ $Trimmed_Rv_FastQ (for bsExpress)"
		# output: $Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam 
		time bismark \
			--bam \
			--bowtie2 \
			--maxins $BM_MAXINS \
			-p $((NT/2)) \
			--temp_dir $PROJECT_DIR/scratch/$Barcode \
			--output_dir $PROJECT_DIR/BsExpress/$Barcode \
			$BS_EXP_REF_DIR \
			-1 $Trimmed_Fw_FastQ \
			-2 $Trimmed_Rv_FastQ

<<bsExpress1
		time bsExpress --verbose \
			--skip_fastqc \
			--skip_trim \
			--skip_shorten \
			--maxlen $BM_MAXINS \
			--input $Trimmed_Fw_FastQ $Trimmed_Rv_FastQ \
			--outdir $PROJECT_DIR/BsExpress/$Barcode \
			--prefix $Barcode \
			--ref $BS_EXP_REF
		#this failed to run, don't know why
bsExpress1

		# sort bam file
		# SLX-10404.A001.C854JANXX.s_4.r_1_val_1.fq.gz_bismark_bt2_pe.bam
		if [ -s $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam ]; then
			echo -e "samtools sort $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.sorted"
			time samtools sort $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.sorted 
		else
			echo -e "\e[031m$PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam not found\e[0m\n"
			exit
		fi

		# index bam and then run bsExpress
		if [ -s $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.sorted.bam ]; then
			echo -e "samtools index $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.sorted.bam"
			time samtools index $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.sorted.bam 

			echo -e "bsExpress $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.bam"
			time bsExpress --verbose \
				--skip_fastqc \
				--skip_trim \
				--skip_shorten \
				--skip_aln \
				--maxlen $BM_MAXINS \
				--input $PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.sorted.bam \
				--outdir $PROJECT_DIR/BsExpress/$Barcode \
				--prefix $Barcode \
				--ref $BS_EXP_REF
		else
			echo -e "\e[031m$PROJECT_DIR/BsExpress/$Barcode/$Barcode.r_1_val_1.fq.gz_bismark_bt2_pe.sorted.bam not found\e[0m\n"
			exit
		fi
