#!/usr/bin/perl

#rescueUnpaired.pl v. 0.0.2

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
	my $zip;
	my $fh = fileopen($infile, \$zip);
	die("Error: Unable to open input fastq $infile\n") if (!$fh);

	do {
		my $header = <$fh>;
		if ($header) {
			chomp($header);
			chomp(my $seq = <$fh>);
			chomp(my $opt = <$fh>);
			chomp(my $qual = <$fh>);
			if (length($seq) >= $minlen) {
				print $outfh "$header\n$seq\n$opt\n$qual\n";
			}
		} else {
			last if eof($fh);
			# If the sequence/base quality string length is zero at the end of a gzipped file, the eof method is unreliable, 
			# therefore, test for any empty header string, which should always be nonempty in a correctly formatted fastq, and then
			# exit normally if eof is true. If a nonempty header is found before the end of the fastq file, somethig must be
			# wrong with the fastq.
			die ("Error: $infile format appears incorrect\n");
		}
	} while (1);

	close $fh;
}

close $outfh;

### subroutines ###

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
