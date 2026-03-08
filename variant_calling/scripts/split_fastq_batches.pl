#!/usr/bin/perl

# split_fastq_batches.pl

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my $version = "0.1.0";
my @pe;
my $outprefix = undef;
my $se = undef;
my $batches = undef;

die(qq/
split_fastq_batches.pl <inputs>

Required input:
--pe   <fastq1 fastq2>   Name of foward read fastq file followed by name of reverse read fastq file.
--outprefix  <STR>   Output file prefix to which batch numbers and file suffixes will be appended.

Optional input:
--single   <fastq>   Name of fastq file containing unpaired reads.
--batches  <file>    List of batch names to extract.

Notes:
- Output forward will be named <out>_<batch number>_R1.fastq.gz
- Output reverse reads will be named <out>_<batch number>_R2.fastq.gz
- Output unpaired reads will be named <out>_<batch number>_U.fastq.gz
- If using --batches, batch names must be formatted as <instrumenet name>:<run ID>:<flowcell ID>:<lane number> as provided in the fastq header.

version $version
\n/) if (!@ARGV || scalar @ARGV < 4);

my $res = GetOptions('pe=s{2}' => \@pe, 'outprefix=s' => \$outprefix, 'single=s' => \$se, 'batches=s' => \$batches);
die("--> exiting prematurely\n") unless ($res);

print STDOUT "\nFiles containing reads:\n";
print STDOUT "Forward reads: $pe[0]\nReverse reads: $pe[1]\n" if @pe;
print STDOUT "Unpaired reads: $se\n" if $se;

# Find all of the unique run, flowcell and lane combinations in input fastq files or user-supplied list of batch names.
# Note that for paired-end reads the batch searching could be sped up by assuming that forward and reverse files
# are synched (only look in one of the two files), but avoiding that assumption for safety.

my @fqfiles;
@fqfiles = @pe if @pe;
push @fqfiles, $se if $se;
my @suffix;
if (scalar @fqfiles == 1) {
	@suffix = ("U.fastq.gz");
} elsif (scalar @fqfiles == 2) {
	@suffix = ("R1.fastq.gz", "R2.fastq.gz");
} elsif (scalar @fqfiles == 3) {
	@suffix = ("R1.fastq.gz", "R2.fastq.gz", "U.fastq.gz");
} else {
	die("Error: Invalid number of input fastq files\n");
}

my %batch_hash;
if ($batches) {
	open(my $bfh, '<', $batches) or die("Error: unable to open batch ID list $batches\n");
	while(<$bfh>) {
		chomp;
		$batch_hash{$_} = 1;
	}
	close $bfh;
} else {
	foreach my $infile (@fqfiles) {
		my $zip;
		my $fh = fileopen($infile, \$zip);
		die("--> exiting prematurely\n") if (!$fh);

		while(1) {
			my $header = <$fh>;
			if ($header) {
				my $batch_id = $1 if $header =~ /^@([^:]+:[^:]+:[^:]+:[^:]):/;
				$batch_hash{$batch_id} = 1;
				for (my $i = 0; $i < 3; $i++) {
					<$fh>;
				}
			} else {
				last if eof($fh);
				die ("Error: $infile format appears incorrect\n");
			}
		}

		close $fh;
	}
}

my @batch_arr = sort {$a cmp $b} keys %batch_hash;
print STDOUT "\nDiscovered batches:\n";
foreach (@batch_arr) {
	print STDOUT "$_\n";
}
print STDOUT "\n";

# extract batches from input fastq files

my $n = 1;
my $iszip;
my $grepcmd;
my $testfh;
foreach my $bid (@batch_arr) {
	my $c = 0;
	foreach my $fq (@fqfiles) {
		$testfh = fileopen($fq, \$iszip);
		die("--> exiting prematurely\n") if (!$testfh);
		close $testfh;
		$grepcmd = $iszip ? "zgrep" : "grep";
		my $out = "${outprefix}_${n}_${suffix[$c]}";
		my $cmd = "$grepcmd -A3 $bid $fq | gzip > $out";
		print STDOUT "Extracting $bid reads from $fq\n";
		system($cmd);
		$c++;
	}	
	$n++;
}

exit;

sub fileopen {
	my ($f, $iszip) = @_;
	my $fh = undef;
	if(open($fh, '<', $f)) {
		sysread $fh, my $magic, 2;
		close $fh;
		if ($magic eq "\x1f\x8b") {
			$fh = new IO::Zlib;
			$fh->open($f, "rb");
			$$iszip = 1;
		} else {
			open($fh, '<', $f);
			$$iszip = 0;
		}
	} else {
		print STDERR "Error: Unable to open file $_[0]\n";
		$$iszip = -1;
		return 0;
	}
	return $fh;
}

