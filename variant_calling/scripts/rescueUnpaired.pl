#!/usr/bin/perl

#rescueUnpaired.pol

use warnings;
use strict;
use Getopt::Long;
use IO::Zlib;
use IO::Compress::Gzip qw($GzipError);

my $fq1 = undef;
my $fq2 = undef;
my $outname = undef;
my $minlen = 70;

die(qq/
rescueUnpaired.pl [options]

Options
--fq1       FILE     Forward read input fastq file.
--fq2       FILE     Reverse read input fastq file.
--outname   STRING   Output file name \(output is compressed with gzip\).
--minlen    INT      Minimum sequence length to retain read. [$minlen]
\n/) if (!@ARGV || scalar @ARGV < 2);

GetOptions('fq1=s' => \$fq1, 'fq2=s' => \$fq2, 'outname=s' => \$outname, 'minlen=i' => \$minlen);

my @fqfiles;
if  ($fq1) {
	if (-f $fq1) {
		push @fqfiles, $fq1;
	} else {
		die("Error: Read 1 file, $fq1, does not exist\n");
	}
}

if ($fq2) {
	if (-f $fq2) {
		push @fqfiles, $fq2;
	} else {
		die("Error: Read 2 file, $fq2, does not exist\n");
	}
}

die("Error: Invalid minimum read length, $minlen\n") if ($minlen < 1);

my $outfh = new IO::Compress::Gzip($outname) or die("Error: Unable to open output file -> $GzipError\n");

foreach my $infile (@fqfiles) {
	my $fh = fileopen($infile);
	die("Error: Unable to open input fastq $infile\n") if (!$fh);

	do {
		chomp(my $header = <$fh>);
		chomp(my $seq = <$fh>);
		chomp(my $opt = <$fh>);
		chomp(my $qual = <$fh>);
		if (length($seq) >= $minlen) {
			print $outfh "$header\n$seq\n$opt\n$qual\n";
		}
	} while (!eof $fh);

	close $fh;
}

close $outfh;

### subroutines ###

sub fileopen {
        my $fh = undef;
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
                print STDERR "Error: Unable to open file $_[0]\n";
                return 0;
        }
        return $fh;
}
