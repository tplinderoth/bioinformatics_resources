#!/usr/bin/perl

# splitFastqLanes.pl <library name, e.g. 'T_090_CSFP210021209-1a'> <directory of reads> <output directory>

use warnings;
use strict;
use IO::Zlib;

die(qq/
splitFastqLanes.pl <library name, e.g. 'T_090_CSFP210021209-1a'> <directory of reads> <output directory>

This script works for reads with headers conforming to the casava 1.8 format.
\n/) if (!@ARGV || scalar @ARGV < 3);

my $lib = $ARGV[0];
my $fqdir = $ARGV[1];
my $outdir = $ARGV[2];

$fqdir =~ s/\/$//;
$outdir =~ s/\/$//;

my $r1fq = "${fqdir}/${lib}_R1.fastq.gz";
my $r2fq = "${fqdir}/${lib}_R2.fastq.gz";

# find all of the unique flowcell and lane combinations in fastq file

my %id;
my $fqfh = fileopen($r1fq);
die("-->exiting\n") if (!$fqfh);

while (<$fqfh>) {
	if ($_ =~ /^@/) {
		my $idx = $1 if ($_ =~ /^([^:]+:[^:]+:[^:]+:[^:]+):/);
		$id{$idx} = 1 if (!exists $id{$idx});
	}
}

close $fqfh;

# extract out sets of reads from different lanes

foreach my $lane (keys %id) {
	my $n = $1 if ($lane =~ /^[^:]+:[^:]+:[^:]+:(\d+)$/);
	my $out1 = "${outdir}/${lib}_${n}_R1.fastq.gz";
	my $out2 = "${outdir}/${lib}_${n}_R2.fastq.gz";
	my $cmd1 = "zgrep -A3 $lane $r1fq | gzip > $out1";
	my $cmd2 = "zgrep -A3 $lane $r2fq | gzip > $out2";
	print STDOUT "Extracting $lane reads from $r1fq\n";
	system($cmd1);
	print STDOUT "Extracting $lane reads from $r2fq\n";
	system($cmd2);
}

exit;

sub fileopen {
	my $fh;
	if(open($fh, '<', $_[0])) {
 		sysread $fh, my $magic, 2;
		close $fh;
		if ($magic eq "\x1f\x8b") {
			$fh = new IO::Zlib;
			$fh->open($_[0], "rb");
		} else {
			open($fh, '<', $_[0]);
		}
	} else {
		print STDERR "ERROR: Unable to open file $_[0]\n";
		return 0;
	}
	return $fh;
}

