#!/usr/bin/perl

#find_sequencing_batches.pl v. 0.0.1

use warnings;
use strict;
use IO::Zlib;
use IO::Compress::Gzip qw($GzipError);

die(qq/
find_sequencing_batches.pl <fastq 1> [fastq 2 ... fastq N]
\n/) if (!@ARGV);

my @fqfiles = @ARGV;
my %batches;

foreach my $infile (@fqfiles) {
	my $zip;
	my $fh = fileopen($infile, \$zip);
	die("--> exiting prematurely\n") if (!$fh);

	while(1) {
		my $header = <$fh>;
		if ($header) {
			my $batch_id = $1 if $header =~ /^@([^:]+:[^:]+:[^:]+:[^:]):/;
			$batches{$batch_id} = 1;
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

foreach (sort {$a cmp $b} keys %batches) {
	print STDOUT "$_\n";
}

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
