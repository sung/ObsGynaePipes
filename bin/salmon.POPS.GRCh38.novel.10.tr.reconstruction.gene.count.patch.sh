#!/bin/bash 

IS_SE=1 # 1:single-end, 0:paired-end
SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' 
TR_PREFIX=GRCh38 # GRCh37|GRCh38
ENS_VER=82 # 82 (for manual downlaod)
SALMON_INDEX=$HOME/data/Salmon/transcriptome_index/$SPECIES/Ensembl/$TR_PREFIX/POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.salmon.idx
SALMON_LIBTYPE=A # automatically infer the library type
GTF=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.gtf
GTF2=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.xloc.txt
GTF3=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.hgnc_symbol.txt
GTF4=/home/ssg29/results/RNA-Seq/Placentome/Cuffcompare/POPS/POPS.GRCh38.novel.10.tr.reconstruction.ensg.txt
NT=12

######################
## Placenta Samples ##
######################
# 1. FG
#SLX_IDS="SLX-9168 SLX-9169" 
# 2. JD
SLX_IDS="SLX-10281 SLX-10402 SLX-9792 SLX-10284 SLX-10285 SLX-10283 SLX-10287"
# 3. Plasma ultra-deep PE151
#SLX_IDS="SLX-11368 SLX-11369" 

####################
## Plasma Samples ##
####################
for SLX in $SLX_IDS; do
	#PROJECT_DIR=$HOME/results/$SLX.$SPECIES.SE125.v2 # FG-placenta
	PROJECT_DIR=$HOME/results/$SLX.$SPECIES.v1 # JD-placenta or Plasma
    ALL_BARCODES=($(for i in `ls $HOME/rcs/rcs-ssg29-obsgynae/POPS/data/fastq/$SLX/$SLX*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f10 | cut -d'.' -f2 ; done | uniq))
	for Barcode in ${ALL_BARCODES[*]}; do
		echo -e "$SLX:$Barcode"
		SALMON_OUT=$PROJECT_DIR/Salmon/GRCh38.82.novel.10/$Barcode; mkdir -p $SALMON_OUT
		TARGET=$SALMON_OUT/$SLX.$Barcode.quant.POPS.$TR_PREFIX.$ENS_VER.novel.10.tr.reconstruction.clean.txt

		if [ -s $TARGET ]; then
			#join <(sort -k1,1 $TARGET) <(sort -k1,1 $GTF3) | sort -k 3,3 | awk 'BEGIN{OFS="\t"}{gene_count[$3]+=$2}END{for (gene in gene_count) print gene,gene_count[gene]}' | sort -k1,1 > ${TARGET%.clean.txt}.hgnc.clean.txt
			join <(sort -k1,1 $TARGET) <(sort -k1,1 $GTF4) | sort -k 3,3 | awk 'BEGIN{OFS="\t"}{gene_count[$3]+=$2}END{for (gene in gene_count) print gene,gene_count[gene]}' | sort -k1,1 > ${TARGET%.clean.txt}.ensg.clean.txt
		else
            echo -e "$TARGET not found\n"
#			if [ $IS_SE -eq 1 ]; then
#				FastQ_array=() #initialise
#				for FastQ_file in `ls /scratch/obsgynae/POPS/data/fastq/$SLX/$SLX.$Barcode*.fq.gz | grep -v lost`; do FastQ_array+=($FastQ_file); done
#				Merged_FastQ=$(printf " %s" "${FastQ_array[@]}")
#				Merged_FastQ=${Merged_FastQ:1} # remove the first space 
#				echo -e "salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -r $Merged_FastQ -g $GTF2 -o $SALMON_OUT"
#				time salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -r $Merged_FastQ -g $GTF2 -o $SALMON_OUT
#			else
#				echo -e "\e[031m Paired-end mode\e[0m\n"
#				Fw_FastQ_array=() #initialise
#				Rv_FastQ_array=() #initialise
#				# Un-trimmed FastQ
#				for Fw_FastQ_file in `ls /scratch/obsgynae/POPS/data/fastq/$SLX/$SLX.$Barcode*.r_1.fq.gz | grep -v lost`; do Fw_FastQ_array+=($Fw_FastQ_file); done
#				for Rv_FastQ_file in `ls /scratch/obsgynae/POPS/data/fastq/$SLX/$SLX.$Barcode*.r_2.fq.gz | grep -v lost`; do Rv_FastQ_array+=($Rv_FastQ_file); done
#				Merged_Fw_FastQ=$(printf " %s" "${Fw_FastQ_array[@]}")
#				Merged_Fw_FastQ=${Merged_Fw_FastQ:1} # remove first space 
#
#				Merged_Rv_FastQ=$(printf " %s" "${Rv_FastQ_array[@]}")
#				Merged_Rv_FastQ=${Merged_Rv_FastQ:1} # remove first space 
#
#				echo -e "salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -1 $Merged_Fw_FastQ -2 $Merged_Rv_FastQ -g $GTF2 -o $SALMON_OUT"
#				time salmon quant -p $NT -i $SALMON_INDEX -l $SALMON_LIBTYPE -1 $Merged_Fw_FastQ -2 $Merged_Rv_FastQ -g $GTF2 -o $SALMON_OUT
#			fi
#
#			awk '/TCONS/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.sf > $TARGET
#			awk '/XLOC/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.genes.sf > ${TARGET%.clean.txt}.xloc.clean.txt
#
#			join <(sort -k1,1 $TARGET) <(sort -k1,1 $GTF3) | sort -k 3,3 | awk 'BEGIN{OFS="\t"}{gene_count[$3]+=$2}END{for (gene in gene_count) print gene,gene_count[gene]}' | sort -k1,1 > ${TARGET%.clean.txt}.hgnc.clean.txt
#			join <(sort -k1,1 $TARGET) <(sort -k1,1 $GTF4) | sort -k 3,3 | awk 'BEGIN{OFS="\t"}{gene_count[$3]+=$2}END{for (gene in gene_count) print gene,gene_count[gene]}' | sort -k1,1 > ${TARGET%.clean.txt}.ensg.clean.txt
#			#echo -e "\e[031m$TARGET not found\e[0m\n"
#			#exit
		fi
	done
done
