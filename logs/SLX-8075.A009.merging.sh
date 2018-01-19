ls: cannot access /home/ssg29/scratch/results/SLX-8075.v1/Bismark/Alignment/A009/SLX-8075.A009*q1.bam: No such file or directory
samtools merge /home/ssg29/scratch/results/SLX-8075.v1/Bismark/Alignment/A009/SLX-8075.A009.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam   

Usage:   samtools merge [-nr] [-h inh.sam] <out.bam> <in1.bam> <in2.bam> [...]

Options: -n       sort by read names
         -r       attach RG tag (inferred from file names)
         -u       uncompressed BAM output
         -f       overwrite the output BAM if exist
         -1       compress level 1
         -l INT   compression level, from 0 to 9 [-1]
         -@ INT   number of BAM compression threads [0]
         -R STR   merge file in the specified region STR [all]
         -h FILE  copy the header in FILE to <out.bam> [in1.bam]

Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users
      must provide the correct header with -h, or uses Picard which properly maintains
      the header dictionary in merging.


real	0m0.023s
user	0m0.001s
sys	0m0.002s
[031m /home/ssg29/scratch/results/SLX-8075.v1/Bismark/Alignment/A009/SLX-8075.A009.r_1_val_1.fq.gz_bismark_bt2_pe.q1.bam not found[0m

