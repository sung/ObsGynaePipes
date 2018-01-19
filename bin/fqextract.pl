#!/usr/bin/perl -w 
# impression from https://www.biostars.org/p/10353/
# and https://github.com/lh3/readfq/blob/master/readfq.pl
# coded by Sung Gong <ssg29@cam.ac.uk>
# crated in 9/Jun/2014

die("Usage: cat in.fq | $0 <max_length_of_read>\n") unless (scalar(@ARGV) == 1);
my $max_len= shift(@ARGV);

sub readfq {
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if (!defined(@$aux));
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
		while (<$fh>) {
			chomp;
			if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
				$aux->[0] = $_;
				last;
			}
		}
		if (!defined($aux->[0])) {
			$aux->[1] = 1;
			return;
		}
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
		chomp;
		$c = substr($_, 0, 1);
		last if ($c eq '>' || $c eq '@' || $c eq '+');
		$seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return ($name, $seq, $qual);
		}
	}
	$aux->[1] = 1;
	return ($name, $seq);
}

my @aux = undef;
my $n=0;
my ($name, $seq, $qual);
while (($name, $seq, $qual) = readfq(\*STDIN, \@aux)) {
	next if length($seq)>=$max_len;
	print "\@$name\n$seq\n+\n$qual\n";
}
