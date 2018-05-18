#!/bin/bash 

SPECIES='Homo_sapiens' # 'Homo_sapiens' or 'Mus_musculus' 
TR_PREFIX=GRCh38 # GRCh37|GRCh38
ENS_VER=82 # 82 (for manual downlaod)
SALMON_INDEX=$HOME/data/Salmon/transcriptome_index/$SPECIES/Ensembl/$TR_PREFIX/$SPECIES.$TR_PREFIX.$ENS_VER.cdna.ncrna.salmon.idx
SALMON_LIBTYPE=A # automatically infer the library type
GENOME=$HOME/data/genome/$SPECIES/Ensembl/$TR_PREFIX/Sequence/WholeGenomeFasta/genome.fa
GTF=$HOME/data/genome/$SPECIES/Ensembl/$TR_PREFIX/Annotation/Genes/$SPECIES.$TR_PREFIX.$ENS_VER.gtf
NT=16

#time gffread $GTF -g $GENOME -w $FASTA 
#time salmon index -t $GTF -i $SALMON_INDEX

SLX_IDS="SLX-ILM-Plasma2017"
for SLX in $SLX_IDS; do
	PROJECT_DIR=$HOME/results/$SLX.$SPECIES.v1
    ALL_BARCODES=($(for i in `ls $HOME/data/fastq/$SLX/$SLX.*.fq.gz | grep -v lost`; do echo $i | cut -d'/' -f7 | cut -d'.' -f2 ; done | uniq))
	for Barcode in ${ALL_BARCODES[*]}; do
		echo -e "$SLX:$Barcode"
		SALMON_OUT=$PROJECT_DIR/Salmon/$TR_PREFIX.$ENS_VER/$Barcode; mkdir -p $SALMON_OUT

        Fw_FastQ_array=() #initialise
        Rv_FastQ_array=() #initialise

		for Fw_FastQ_file in `ls $HOME/data/fastq/$SLX/$SLX.$Barcode*r_1.fq.gz | grep -v lost`; do Fw_FastQ_array+=($Fw_FastQ_file); done
		for Rv_FastQ_file in `ls $HOME/data/fastq/$SLX/$SLX.$Barcode*r_2.fq.gz | grep -v lost`; do Rv_FastQ_array+=($Rv_FastQ_file); done

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
        if [ -s $SALMON_OUT/quant.sf ]; then
            awk '/^ENST/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.sf | sort > $SALMON_OUT/$SLX.$Barcode.quant.$TR_PREFIX.$ENS_VER.enst.clean.txt
        fi
        if [ -s $SALMON_OUT/quant.genes.sf ]; then
            awk '/^ENSG/{printf "%s\t%.f\n", $1,$5}' $SALMON_OUT/quant.genes.sf | sort > $SALMON_OUT/$SLX.$Barcode.quant.$TR_PREFIX.$ENS_VER.ensg.clean.txt
        fi
        echo -e "\e[32m$Barcode done for salmon -quant\e[0m" 
	done
done
