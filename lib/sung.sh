#!/bin/bash

function mkdir_unless(){
	if [ ! -d $1 ]; then
		mkdir -p $1
	fi	
}

function make_run_script(){
    MY_TEMPLATE=$1 # e.g. template/tophat.persample or template/edgeR 
    MY_APP=$2 # e.g. SLX-8547.D704_D503.tophat.persample.sh

    MY_BARCODE=$3 # e.g. D704_D503 (optional)
    MY_CELL=$4 # e.g. C48CWACXX (optional)
    MY_LANE=$5 # e.g. s_1 (optional)
    MY_CHUNK=$6 # e.g. 1 (optional)

    cp $MY_TEMPLATE $BIN_TOP/script/$MY_APP

    sed -i "s/MY_SLX/$SLX/" $BIN_TOP/script/$MY_APP
    sed -i "s/MY_COHORT/$COHORT/" $BIN_TOP/script/$MY_APP
    sed -i "s/MY_VERSION/$VERSION/" $BIN_TOP/script/$MY_APP
	if [ -n "$MY_BARCODE" ];then #string is not null.
    	sed -i "s/MY_BARCODE/$MY_BARCODE/" $BIN_TOP/script/$MY_APP
	fi
	if [ -n "$MY_CELL" ];then #string is not null.
    	sed -i "s/MY_CELL/$MY_CELL/" $BIN_TOP/script/$MY_APP
	fi
	if [ -n "$MY_LANE" ];then #string is not null.
    	sed -i "s/MY_LANE/$MY_LANE/" $BIN_TOP/script/$MY_APP
	fi
	if [ -n "$MY_CHUNK" ];then #string is not null.
    	sed -i "s/MY_CHUNK/$MY_CHUNK/" $BIN_TOP/script/$MY_APP
	fi
}

#http://stackoverflow.com/questions/12147040/division-in-script-and-floating-point/24431665#24431665
div ()  # Arguments: dividend and divisor
{
	if [ $2 -eq 0 ]; then echo division by 0; exit; fi
	local p=12                            # precision
	local c=${c:-0}                       # precision counter
	local d=.                             # decimal separator
	local r=$(($1/$2)); echo -n $r        # result of division
	local m=$(($r*$2))
	[ $c -eq 0 ] && [ $m -ne $1 ] && echo -n $d
	[ $1 -eq $m ] || [ $c -eq $p ] && return
	local e=$(($1-$m))
	let c=c+1
	div $(($e*10)) $2
} 

# from Stuart Rankin to monitor job list
jl(){
	if [[ `hostname` =~ "login-sand" ]]; then
		squeue -p sandybridge --format="%.20i %.20A  %.20a %.20F %.9P %.8j %.8u %.8T %.9M %.9l %.6D %20R %Q" | egrep -ve "QOSResource|AssocGrpCPUMins" | less
	else
		squeue -p skylake --format="%.20i %.20A  %.20a %.20F %.9P %.8j %.8u %.8T %.9M %.9l %.6D %20R %Q" | egrep -ve "QOSResource|AssocGrpCPUMins" | less
	fi
}
