#!/usr/bin/perl

# subsamp_sites.pl

# subsamples sites in even increments

use warnings;
use strict;

my $nsamp = 0;
my $start = 0;

die(qq/
subsamp_sites.pl <start coordinate (0-based)> <number sites to sample> <sites file>

Assumes input file has a header.
\n/) if (!@ARGV || scalar @ARGV < 3);

my $fname = $ARGV[2];
open(my $ifh, '<', $fname) or die("Unable to open file of sites $fname: $!\n");

$start = $ARGV[0];
$nsamp = $ARGV[1];

chomp(my $nsites = `wc -l $fname | cut -f1 -d ' '`);
$nsites = $nsites-$start;
die ("Start coordinate exceeds number of sites in input file: $!\n") if ($nsites < 0);

die("Invalid number of sites to subsample: $nsamp\n") if ($nsamp < 0);

my $inc = ($nsamp == 0) ? 0:int($nsites/$nsamp);
print STDERR "Sampling interval: $inc entries\n";

my $header = <$ifh>;
print STDOUT "$header";

my $k = 0;
for ($k = 0; $k < $start; $k++)  {
	<$ifh>;
}

$k = 0;
while (<$ifh>) {
	$k++;
	if ($k == $inc) {
		print STDOUT $_;
		$k = 0;
	}
}

close $ifh;

exit(0);
