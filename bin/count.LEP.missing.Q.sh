#########################
## All based on GRCh37 ##
#########################
# 1. BED files from Schultz and RoadMap
for i in `ls /home/ssg29/data/RoadMap/Schultz2015/RNA-Seq/*bed.gz`; do echo $i; zcat $i | awk '$1=="chr7" && ($2==127894456 || $2==127894459){print $2}' | sort | uniq -c; done
for i in `ls /home/ssg29/data/RoadMap/Consortium/RNA-Seq/Somatic/*bed.gz`; do echo $i; zcat $i | awk '$1=="chr7" && ($2==127894456 || $2==127894459){print $2}' | sort | uniq -c; done
for i in `ls /home/ssg29/data/RoadMap/Consortium/RNA-Seq/Placenta/*bed.gz`; do echo $i; zcat $i | awk '$1=="chr7" && ($2==127894456 || $2==127894459){print $2}' | sort | uniq -c; done
for i in `ls /home/ssg29/data/RoadMap/Consortium/RNA-Seq/Fetal/*bed.gz`; do echo $i; zcat $i | awk '$1=="chr7" && ($2==127894456 || $2==127894459){print $2}' | sort | uniq -c; done

# 2. BAM files from illumina plasma-2017
# q > 30 : uniqutely aligned reads
# with itself and mate mapped: 0x1 bit set and neither 0x4 nor 0x8 bits set
# no secondary alignment : 0x100
for i in `ls ~/results/Illumina/Plasma-RNA-2017/Bam/PL*bam`;do echo $i; samtools view -q 30 -f 0x1 -F 0x4 -F 0x8 -F 0x100 $i | awk '$3=="chr7" && ($4==127894457 || $4==127894460){print $4}' | sort | uniq -c; done # NONE 
for i in `ls ~/results/Illumina/Plasma-RNA-2017/Bam/RNA*bam`;do echo $i; samtools view -q 30 -f 0x1 -F 0x4 -F 0x8 -F 0x100 $i | awk '$3=="chr7" && ($4==127894457 || $4==127894460){print $4}' | sort | uniq -c; done # SEE BELOW

/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-100_S96_accepted_hits.bam
1 127894457
2 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-102_S97_accepted_hits.bam
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-103_S98_accepted_hits.bam
12 127894457
11 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-105_S99_accepted_hits.bam
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-106_S100_accepted_hits.bam
10 127894457
4 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-107_S101_accepted_hits.bam
1 127894457
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-108_S102_accepted_hits.bam
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-88_S84_accepted_hits.bam
1 127894457
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-89_S85_accepted_hits.bam
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-90_S86_accepted_hits.bam
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-91_S87_accepted_hits.bam
11 127894457
25 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-92_S88_accepted_hits.bam
3 127894457
5 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-93_S89_accepted_hits.bam
1 127894457
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-94_S90_accepted_hits.bam
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-95_S91_accepted_hits.bam
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-96_S92_accepted_hits.bam
2 127894457
1 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-97_S93_accepted_hits.bam
9 127894457
2 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-98_S94_accepted_hits.bam
51 127894457
40 127894460
/home/ssg29/results/Illumina/Plasma-RNA-2017/Bam/RNA-99_S95_accepted_hits.bam

