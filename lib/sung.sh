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

function run_pipeline(){
	MY_CALLER=$1 # 
	USE_SLURM=$2 # 1 (use slurm) or 0 (do not use slurm)

	TEMPLATE=$BIN_TOP/template/$MY_CALLER.sh # template/tophat.persample.sh # per sample 
	if [ ! -f $TEMPLATE ];then
		echo "$TEMPLATE not found"
		exit
	fi

	printf "\e[32mdoing $MY_CALLER for $SLX\e[0m\n"
	#################################
	# per SLX run (count read.base) #
	#################################
	if [[ $MY_CALLER =~ ^per.run ]]; then
		MY_APP=$SLX.$VERSION.$COHORT.$MY_CALLER.sh
		make_run_script $TEMPLATE $MY_APP

		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$VERSION.$COHORT.$MY_CALLER.sh
		make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER

		# submit this job to the scheduler
		if [ $USE_SLURM -eq 1 ]; then
			cd $TOP/logs
			echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
			sbatch $SLURM_OPT $MY_SLURM_SCRIPT
		else
			# run the script locally
			echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.$COHORT.log"
			time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.$COHORT.log &
		fi
	###########################
	# per sample (or barcode) #
	###########################
	elif [[ $MY_CALLER =~ ^per.sample ]] || [[ $MY_CALLER =~ ^per.barcode ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\e[32mdoing $SLX.$VERSION.$MY_BARCODE\e[0m\n"
			FASTQ_FILES=`ls $HOME/data/fastq/$SLX/${SLX}.${MY_BARCODE}.*r_1.fq.gz | grep -v lost`
			MY_CELL=`echo $FASTQ_FILES[0] | cut -d/ -f7 | cut -d. -f3`
            #MY_CELL="XXX"
			MY_APP=$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.$MY_CELL.sh
			make_run_script $TEMPLATE $MY_APP $MY_BARCODE $MY_CELL

			# submit this job to the scheduler
			if [ $USE_SLURM -eq 1 ]; then
				MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.$MY_CELL.sh
				make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $MY_CELL
				cd $TOP/logs
				echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
				sbatch $SLURM_OPT $MY_SLURM_SCRIPT
			# run the script locally
			else
				echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.log"
				time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_BARCODE.$MY_CALLER.log &
			fi
		done # end of MY_BARCODE
	############
	# per lane #
	############
	elif [[ $MY_CALLER =~ '^per.lane' ]]; then
		for MY_BARCODE in ${ALL_BARCODES[*]}
		do
			printf "\n\e[32m$MY_CALLER for $SLX.$MY_BARCODE\e[0m\n"
			FASTQ_FILES=`ls $FASTQ_DIR/${SLX}.${MY_BARCODE}.*r_1.fq.gz | grep -v lost`
			Cell_lane_array=`for i in $FASTQ_FILES; do echo $i | cut -d/ -f7 | cut -d. -f3,4 ; done`
			for MY_CELL_LANE in ${Cell_lane_array[*]} 
			do
				MY_CELL=`echo $MY_CELL_LANE | cut -d. -f1`
				MY_LANE=`echo $MY_CELL_LANE | cut -d. -f2`
				MY_APP=$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.sh
				make_run_script $TEMPLATE $MY_APP $MY_BARCODE $MY_CELL $MY_LANE

				# submit this job to the scheduler
				if [ $USE_SLURM -eq 1 ]; then
					MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.sh
					make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER $MY_BARCODE $MY_CELL $MY_LANE 1
					cd $TOP/logs
					echo -e "sbatch $MY_SLURM_SCRIPT"
					sbatch $MY_SLURM_SCRIPT
				# run the script locally
				else
					echo -e "bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.log"
					nice -n 10 time bash $BIN_TOP/script/$MY_APP >& $TOP/logs/$SLX.$MY_BARCODE.$MY_CELL.$MY_LANE.$MY_CALLER.log &
				fi # end of USE_SLURM
			done # end of MY_CELL_LANE
		done
	#####################
	# per batch (edgeR) #
	#####################
	elif [[ $MY_CALLER =~ ^edgeR ]]; then
		MY_APP=$SLX.$VERSION.$MY_CALLER.sh
		make_run_script $TEMPLATE $MY_APP

		MY_SLURM_SCRIPT=$BIN_TOP/slurm/slurm.$SLX.$VERSION.$MY_CALLER.sh
		make_slurm_script $MY_SLURM_SCRIPT $MY_APP $MY_CALLER

		# submit this job to the scheduler
		if [ $USE_SLURM -eq 1 ]; then
			cd $TOP/logs
			echo -e "sbatch $SLURM_OPT $MY_SLURM_SCRIPT"
			sbatch $SLURM_OPT $MY_SLURM_SCRIPT
		# run the script locally
		else
			echo -e "bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.log"
			time bash $BIN_TOP/script/$MY_APP &> $TOP/logs/$SLX.$VERSION.$MY_CALLER.log &
		fi
	else
		echo -e "\e[031m$MY_CALLER not found\e[0m\n"
		exit
	fi # end of if MY_CALLER
} # end of run_pipeline



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
