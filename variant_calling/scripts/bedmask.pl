#!/usr/bin/perl

# bedmask.pl

# updated 2023-03-13 to work on input with multiple chromosomes

use warnings;
use strict;
use Getopt::Long;

my $mindepth = 0;
my $maxdepth = "inf";
my $minmq = 35;
my $maxmq0 = 0.2;
my $minbq = 20;
my $maxbq0 = 1;

die(qq/
** version 0.2 **

bedmask.pl [options] <input bamstats file> <output file prefix>

Options:
--mindepth   INT   Minimum total depth [$mindepth]
--maxdepth   INT   Maximum total depth [$maxdepth]
--minmq      FLOAT Minimum RMS map quality [$minmq]
--maxmq0     FLOAT Maximum fraction of map quality zero reads [$maxmq0]
--minbq	     FLOAT Minimum RMS base quality [$minbq]
--maxbq0     FLOAT Maximum fraction of base quality zero reads [$maxbq0]
\n/) if (!@ARGV || scalar(@ARGV) < 2);


my $outprefix = pop @ARGV;
my $infile = pop @ARGV;
die ("Input bamstats not a regular file: $infile\n") if (!-f $infile);

GetOptions('mindepth=i' => \$mindepth, 'maxdepth=i' => \$maxdepth, 'minmq=f' => \$minmq, 
'maxmq0=f' => \$maxmq0, 'minbq=f' => \$minbq, 'maxbq0=f' => \$maxbq0);

my $outpass = $outprefix . "_pass.bed";
my $outfail = $outprefix . "_fail.bed";

open(my $infh, '<', $infile) or die("Unable to open input bamstats file: $infile\n");
open(my $passfh, '>', $outpass) or die("Unable to open output file of PASS sites: $outpass\n");
open(my $failfh, '>', $outfail) or die("Unable to open output file of FAIL sites: $outfail\n");

my $userargs = "\nInput file: $infile\nPassed sites bed: $outpass\nFailed sites bed: $outfail
Min Depth: $mindepth\nMax Depth: $maxdepth\nMinimum RMS map quality: $minmq\nMax fraction MQ0 reads: $maxmq0
Minimum RMS BQ: $minbq\nMaximum fraction BQ0 reads: $maxbq0\n";

print STDOUT "$userargs\n";

my @fnames = ('NoCov', 'LowCov', 'ExcessCov', 'LowBQ', 'ExcessBQ0', 'LowMQ', 'ExcessMQ0');
my $nsites = 0;
my $nfail = 0;
my %counts;
foreach (@fnames) {
	$counts{$_} = 0;
}
my %flags;

my ($chr, $start, $end, $flast, $fstr) = ('', 0, 0, '', '');
my @l = split(/\t/, <$infh>);
while (<$infh>) {
	chomp;
	@l = split(/\t/, $_);
	$nsites++;

	# check for new chromosome or new region
	if (!$chr || ($chr && $l[0] ne $chr) || ($chr && $l[1]-1 > $end)) {
		if ($chr) {
			if ($flast eq "PASS") {
				print $passfh "$chr\t$start\t$end\n";
			} else {
				print $failfh "$chr\t$start\t$end\t$flast\n";
			}
		}
		$flast = '';
	}

	my @filters;
	if ($l[2] > 0) {
		if ($l[2] < $mindepth) {
			push @filters, 'LowCov';
			$counts{LowCov}++;
		}
		if ($l[2] > $maxdepth) {
			push @filters, 'ExcessCov';
			$counts{ExcessCov}++;
		}
		if ($l[4] < $minbq) {
			push @filters, 'LowBQ';
			$counts{LowBQ}++;
		}
		if ($l[5] > $maxbq0) {
			push @filters, 'ExcessBQ0';
			$counts{ExcessBQ0}++;
		}
		if ($l[7] < $minmq) {
			push @filters, 'LowMQ';
			$counts{LowMQ}++;
		}
		if ($l[8] > $maxmq0) {
			push @filters, 'ExcessMQ0';
			$counts{ExcessMQ0}++;
		}
	} else {
		push @filters, 'NoCov';
		$counts{NoCov}++;
	}

	if (@filters) {
		$fstr = join(';', @filters);
		$nfail++;
		$flags{$fstr} = 0 if (!exists $flags{$fstr});
		$flags{$fstr}++;
	} else {
		$fstr = "PASS";
	}
	
	if (!$flast) {
		# this case just happens at the start of a new chromosome or region
		($chr, $start, $end, $flast) = ($l[0], $l[1]-1, $l[1], $fstr);
	} elsif ($fstr eq $flast) {
		$end++;
	} elsif ($fstr ne $flast) {
		if ($flast eq "PASS") {
			print $passfh "$chr\t$start\t$end\n";
		} else {
			print $failfh "$chr\t$start\t$end\t$flast\n";
		}
		($chr, $start, $end, $flast) = ($l[0], $l[1]-1, $l[1], $fstr);
	}

}

# print last region
if ($flast eq "PASS") {
	print $passfh "$chr\t$start\t$end\n";
} else {
	print $failfh "$chr\t$start\t$end\t$flast\n";
}

close $infh;
close $passfh;
close $failfh;

print STDOUT "===== REPORT =====\nNumber sites analyzed: $nsites\nNumber FAIL sites: $nfail\n";
print STDOUT "\n===== FAIL TYPES =====\n";
foreach (keys %flags) {
        print STDOUT "$_: $flags{$_}\n";
}
print STDOUT "\n===== FILTER COUNTS (among fail types) =====\n";
foreach my $id (@fnames) {
	print STDOUT "$id: $counts{$id}\n";
}

exit;
