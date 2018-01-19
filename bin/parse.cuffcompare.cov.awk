BEGIN{FS="\t"}
{
	track="";
	for (i = 5; i <= NF; i++){
		if($i=="-"){
			track=track",NA";
		}else{
			cnt=0;sum=0;split($i,multi,","); # multi-trascript annotation for this sample
			for(m in multi){
				split(multi[m],info,"|");
				if(info[7]!="-"){ # info[7]: coverage of this transcript
					sum+=info[7];
				}
				cnt++
			}
			track=track","sum/cnt;
		}
	}
	print $1 track
}
