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
				if(info[8]!="-"){
					sum+=info[8];
				}
				cnt++
			}
			track=track","sum/cnt;
		}
	}
	print $1 track
}
