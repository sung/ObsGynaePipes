#!/bin/bash
# Author: Sung Gong <ssg29@cam.ac.uk>
# First created: 30/Jan/2015
# Last modified: 19/Mar/2018
# A wrapper script to run BWA and GATK for WGS (e.g. HiSeq-Xten)
# Optimised and customised to run at the UIS CSD3 (eg. peta4 skylake) machines

source $HOME/Pipelines/config/target.seq.sbs.config # calls global.config, slurm.config, samples.config 
source $HOME/Pipelines/lib/sung.sh # defines user defined functions (e.g. make_run_script)
# 'make_slurm_script' defined within config/slurm.config # defines 'make_slurm_script'

mkdir -p $RESULT_DIR	
mkdir -p $TOP/logs
mkdir -p $BIN_TOP/script
mkdir -p $BIN_TOP/slurm

#run_pipeline "per.run.count.read.base.fq" 0 # per SLX

#run_pipeline "per.lane.trim.bwa" 0 # per lane
#run_pipeline "per.barcode.trim.bwa" 1 # per barcode (sample) 
#run_pipeline "per.barcode.gatk" 1 # based on GATK3
#run_pipeline "per.barcode.gatk4" 1 # based on GATK4 via Spark (beta; not perfect)
#run_pipeline "per.run.GenotypeGVCFs" 1 # based on GATK3
#run_pipeline "per.run.VQSR" 1 # based on GATK3

#run_pipeline "per.barcode.coverage" 0 # per barcode (sample)
#run_pipeline "per.barcode.gatk.somatic" 0 # per barcode (sample)
#run_pipeline "per.barcode.mutect" 0 # per barcode (sample)
#run_pipeline "per.barcode.somatic.snp" 0 # per barcode (sample)
