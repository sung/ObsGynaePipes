#!/usr/bin/bash

NT1=16
NT2=$((NT1/2))
NT3=$(($NT1*2))
NT4=$(($NT1%2))

echo $NT1
echo $NT2
echo $NT3
echo $NT4

BM_BUFFER_SIZE=60G 
echo $((${BM_BUFFER_SIZE%G}/20))G # 3

sum=0
for i in {1..10}; 
do
	let "sum = sum + 1"	
done
echo -e "sum=$sum"

TEMP_CHR=(chr7 chr12 chr18 chr17 chr3 chr13 chr10 chr4 chr2 chr15 chrX chr21 chr6 chr22 chr16)
for Chr in ${TEMP_CHR[*]}
do 
	echo -e "chr=$Chr"
done

source $HOME/Pipelines/test.config
#CASECTR_PAIRS=(D705-D501,D706-D505 D705-D507,D707-D503 D703-D505)

## initialise 'ALL_BARCODES'
for MY_CASECTR in ${CASECTR_PAIRS[*]}
do
	dummy_array=(${MY_CASECTR//,/ }) # (D709_D501,D709_D502,D101_D102) => (D709_D501 D709_D502 D101_D102)
	for dummy_element in ${dummy_array[*]}
	do
		ALL_BARCODES+=($dummy_element) #push it
	done
done

for MY_BARCODE in ${ALL_BARCODES[*]}
do
	printf "\e[32mdoing $MY_BARCODE\e[0m\n"
done


SLX="SLX-8549"
FLOWCELL="C4G2JANXX"
file="SLX-8549.D707_D503.C4G2JANXX.s_1.r_1.fq.gz"
if [[ "$file" =~ "$SLX" ]] && [[ "$file" =~ "$FLOWCELL" ]]; then
	echo yes
else
	echo no
fi
