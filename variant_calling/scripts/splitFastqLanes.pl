#!/usr/bin/perl

# splitFastqLanes.pl [string to append to library name] <library name, e.g. 'T_090_CSFP210021209-1a'> <directory of reads> <output directory>

use warnings;
use strict;
use IO::Zlib;

die(qq/
splitFastqLanes.pl [string to append to library name] <library name, e.g. 'T_090_CSFP210021209-1a'> <directory of reads> <output directory>

This script works for reads with headers conforming to the casava 1.8 format.

The output fastq files will be named <library><name modifier if provided>_<batch number>_R[1|2].fastq.gz
\n/) if (!@ARGV || scalar @ARGV < 3);

my $outdir = pop @ARGV;
my $fqdir = pop @ARGV;
my $lib = pop @ARGV;

my $lib_modifier = @ARGV ? $ARGV[0] : '';

$fqdir =~ s/\/$//;
$outdir =~ s/\/$//;

my $r1fq = "${fqdir}/${lib}_R1.fastq.gz";
my $r2fq = "${fqdir}/${lib}_R2.fastq.gz";

# find all of the unique flowcell and lane combinations in fastq file

my %id;
my $fqfh = fileopen($r1fq);
die("-->exiting\n") if (!$fqfh);

do {
	my $fqline = <$fqfh>;
	if ($fqline =~ /^@/) {
		my $idx = ($fqline =~ /^([^:]+:[^:]+:[^:]+:[^:]+):/) ? $1 : '';
		die ("Error: Unrecognized fastq header format:\n$idx\n") if (!$idx);
		$id{$idx} = 1 if (!exists $id{$idx});
	} else {
		die("Error: Expected read header did not start with \'@\':\n$fqline\n");
	}
	for (my $i=0; $i<3; $i++) {
		$fqline = <$fqfh>;
	}
} while (!eof $fqfh);

close $fqfh;

# extract out sets of reads from different lanes

my $n = 1;
foreach my $lane (keys %id) {
	#my $n = $1 if ($lane =~ /^[^:]+:[^:]+:[^:]+:(\d+)$/);
	my $out1 = "${outdir}/${lib}${lib_modifier}_${n}_R1.fastq.gz";
	my $out2 = "${outdir}/${lib}${lib_modifier}_${n}_R2.fastq.gz";
	my $cmd1 = "zgrep -A3 $lane $r1fq | gzip > $out1";
	my $cmd2 = "zgrep -A3 $lane $r2fq | gzip > $out2";
	print STDOUT "Extracting $lane reads from $r1fq\n";
	system($cmd1);
	print STDOUT "Extracting $lane reads from $r2fq\n";
	system($cmd2);
	$n++;
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

