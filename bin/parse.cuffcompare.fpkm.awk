#http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/#tracking-transfrags-through-multiple-samples-outprefixtracking
BEGIN{FS="\t"}
{
	track="";
	for (i = 5; i <= NF; i++){
		if($i=="-"){
			track=track",NA";
		}else{
			cnt=0;sum=0;split($i,multi,",");
			for(m in multi){
				split(multi[m],info,"|");
				if(info[4]!="-"){
					sum+=info[4];
				}
				cnt++
			}
			track=track","sum/cnt;
		}
	}
	print $1 track
}
