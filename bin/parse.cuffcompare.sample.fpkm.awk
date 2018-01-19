# http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/#tracking-transfrags-through-multiple-samples-outprefixtracking
# NB. multi-transcript per tcon:
# q2:SLX-10283.D711_D502.STRG.23877|SLX-10283.D711_D502.STRG.23877.2|100|1.438221|0.000000|0.000000|17.084448|1571,SLX-10283.D711_D502.STRG.23877|SLX-10283.D711_D502.STRG.23877.1|100|1.034284|0.000000|0.000000|12.286129|3420
BEGIN{FS="\t";OFS=","}
{
	track="";
	for (i = 5; i <= NF; i++){
		if($i=="-"){
			#track=track",NA";
		}else{
			split($i,multi,","); # multi-trascript annotation for this sample
			for(m in multi){
				split(multi[m],info,"|");  # info[4]: FPKM, info[7]:coverage, info[8]:length
				split(info[2],sample,"."); # info[1]: q2:SLX-10283.D711_D502.STRG.23877, info[2]: SLX-10283.D711_D502.STRG.23877.2
				if(info[8]=="-")info[8]="NA";
				print $1,sample[1],sample[2],sample[3]"."sample[4]"."sample[5],info[4],info[7],info[8]
			}
		}
	}
}
# OUTPUT
# TCONS_00323031,SLX-10283,D711_D502,STRG.23877.2,1.438221,17.084448,1571
# TCONS_00323031,SLX-10283.D711_D502,STRG.23877.1,1.034284,12.286129,3420
