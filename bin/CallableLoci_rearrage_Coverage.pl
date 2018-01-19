# awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' 14so00002DN.q8.OnTarget.bases.callable | intersectBed -a Fluidigm_Target.bed -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse | perl callable.pl 
# Input
#chr21   35742756        35743170        KCNE2_00001043939_1     LOW_COVERAGE=5,CALLABLE=370,LOW_COVERAGE=36
#chr21   35821481        35821977        KCNE1_00001537342_1     CALLABLE=494,NO_COVERAGE=1
#chr3    38591760        38593059        SCN5A_00001672933_1;SCN5A_00001672933_2;SCN5A_00001672933_3     CALLABLE=1298

# Output
#chr21	35742756	35743170	KCNE2_00001043939_1	371	42	0
#chr21	35821481	35821977	KCNE1_00001537342_1	495	0	2
#chr3	38591760	38593059	SCN5A_00001672933_1;SCN5A_00001672933_2;SCN5A_00001672933_3	1299	0	0

#!/usr/bin/perl -w
use strict;
my %data=();
my $name;
my $val;
my @key=();
my ($input) = @ARGV;
open(INPUT,"$input") or die "Can't create log file";

print "Chr\tStart\tEnd\tGene\tCALLABLE\tLOW_COVERAGE\tNO_COVERAGE\n";
while(<INPUT>)
	{
	chomp();
	my @col_name = split ("\t",$_);
	my @info = split ("\,",$col_name[4]);

	#print "\e[032m",$col_name[4],"\e[0m\n";

	my %data;
	foreach my $entry (@info){
		my ($key,$val)=split(/=/,$entry);	
		if(defined $data{$key}){
			#remove $val=$val+1
			$val=$val;
			$data{$key}=$data{$key}+$val;
		}else{
			$data{$key}=$val;
		}
	}

	##Take out increment +1 ++$data{CALLABLE}:0; 
	$data{CALLABLE} = exists $data{CALLABLE}? $data{CALLABLE}:0;
	$data{LOW_COVERAGE} = exists $data{LOW_COVERAGE}? $data{LOW_COVERAGE}:0;
	$data{NO_COVERAGE} = exists $data{NO_COVERAGE}? $data{NO_COVERAGE}:0;
	print "$col_name[0]\t$col_name[1]\t$col_name[2]\t$col_name[3]\t$data{CALLABLE}\t$data{LOW_COVERAGE}\t$data{NO_COVERAGE}\n";

}

