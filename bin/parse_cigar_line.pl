#!/usr/bin/perl -w 
#===============================================================================
#
#         FILE:  parse_cigar_line.pl
#
#        USAGE:  ./parse_cigar_line.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Sungsam Gong (sung), sung@bio.cc
#      VERSION:  1.0
#      CREATED:  19/06/14 17:43:46
#     REVISION:  ---
#===============================================================================

while (<STDIN>){
	@a=split/\s+/;
	#$a[5]=~/\*|([0-9]+[MIDNSHPX=])+/;

	@b=split(/([0-9]+[MIDNSHPX=])/,$a[5]);
	$bad=0;
	foreach (@b){
		if(/([0-9]+)[IDX]/){ #insertion, deletion, or mismatch
			$bad=$bad+1;
		}
	}
	print "$a[4]\t$bad\n"; # 'MAPQ' 'NO of I|D|X'

}

