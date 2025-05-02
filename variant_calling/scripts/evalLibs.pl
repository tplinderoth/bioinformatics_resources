#!/usr/bin/perl

# evalLibs.pl

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my @raw;
my @clean;
my $libname;
my $maxlen = 299; # max read length

die(qq/
evalLibs.pl [options] --raw <raw.fastq1 raw.fastq2 .. raw.fastqn> --clean <clean.fastq1 clean.fastq2 .. clean.fastqn>

Options
--library  STR  Name of library
--maxlen   INT  Maximum read length (bp) [$maxlen]
\n/) if !@ARGV;

GetOptions('raw=s{1,}' => \@raw, 'clean=s{1,}' => \@clean, 'library=s' => \$libname, 'maxlen=i' => \$maxlen);

# initialize length distribution hash
my %sizedist;
for (my $i=0; $i <= $maxlen; $i++) {
	$sizedist{$i} = 0;
}

print STDOUT "## LIBRARY $libname ##\n\n" if ($libname);

print STDOUT "--RAW FASTQ--\n";
my ($raw_reads, $raw_bases, $raw_n) = (0,0,0);
my $perc_n;
foreach my $rawfq (@raw) {
	print STDOUT "$rawfq counts\n";
	my $fqfh = fileopen($rawfq);
	die("-->exiting\n") if (!$fqfh);
	my ($libreads, $libbases, $libn) = fastqCounts($fqfh, undef);
	$perc_n = sprintf '%.2f', $libn/$libbases * 100;
	print STDOUT "number reads: $libreads\nnumber bases: $libbases\nnumber N bases: $libn ($perc_n%)\n";
	$raw_reads += $libreads;
	$raw_bases += $libbases;
	$raw_n += $libn;
	close $fqfh;
}

$perc_n = sprintf '%.2f', $raw_n/$raw_bases * 100;
print STDOUT "\nRaw library totals\n";
print STDOUT "raw number reads: $raw_reads\nraw number bases: $raw_bases\nraw number N bases: $raw_n ($perc_n%)\n";

print STDOUT "\n--CLEAN FASTQ--\n";
my ($clean_reads, $clean_bases, $clean_n) = (0,0,0);
foreach my $cleanfq (@clean) {
	print STDOUT "$cleanfq counts\n";
	my $fqfh = fileopen($cleanfq);
	die("-->exiting\n") if (!$fqfh);
	my ($libreads, $libbases, $libn) = fastqCounts($fqfh, \%sizedist);
	$perc_n = sprintf '%.2f', $libn/$libbases * 100;
	print STDOUT "number reads: $libreads\nnumber bases: $libbases\nnumber N bases: $libn ($perc_n%)\n";
	$clean_reads += $libreads;
        $clean_bases += $libbases;
        $clean_n += $libn;
	close $fqfh;
}

$perc_n = sprintf '%.2f', $clean_n/$clean_bases * 100;
print STDOUT "\nClean library totals\n";
print STDOUT "clean number reads: $clean_reads\nclean number bases: $clean_bases\nclean number N bases: $clean_n ($perc_n%)\n";

print STDOUT "\nData removed by quality control\n";
my $rmv_reads = $raw_reads - $clean_reads;
my $read_perc = $raw_reads > 0 ? sprintf '%.2f', $rmv_reads/$raw_reads * 100 : 0;
my $rmv_bases = $raw_bases - $clean_bases;
my $base_perc = $raw_bases > 0 ? sprintf '%.2f', $rmv_bases/$raw_bases * 100 : 0;
print STDOUT "filtered reads: $rmv_reads ($read_perc%)\nfiltered bases: $rmv_bases ($base_perc%)\n";

# print cleaned read length distribution
print STDOUT "\n--CLEAN READ LENGTHS--\n";
print STDOUT "LENGTH\tCOUNT\n";
my @readlen = sort { $a <=> $b } keys %sizedist;
for (my $i = 0; $i <= $readlen[$#readlen]; $i++) {
	my $c = exists $sizedist{$i} ? $sizedist{$i} : 0;
	print STDOUT "$i\t$c\n";
}

exit;

sub tallySize {
	my ($counts, $len) = @_;
	if (exists $$counts{$len}) {
		$$counts{$len}++;
	} else {
		$$counts{$len} = 1;
	}
}

sub fastqCounts {
	my ($fh, $lcounts) = @_;
	my ($nreads, $nbases, $nn) = (0,0,0);
	<$fh>; # skip first line
	chomp(my $line = <$fh>);
	$nreads++;
	$nbases += length($line);
	$nn += () = $line =~ /N/ig; # count number of 'N' bases
	tallySize($lcounts, $nbases) if ($lcounts);
	while (<$fh>) {
		my $n = 0;
		while ($n < 3 && ($line = <$fh>)) { $n++; } # faster than modulo (% 4 == 0) test
		if ($line) {
			chomp($line);
			$nreads++;
			my $len = length($line);
			$nbases += $len;
			$nn += () = $line =~ /N/ig;
			tallySize($lcounts, $len) if ($lcounts);
		}
	#print STDERR "$nreads\t$nbases\n" if ($nreads % 1000 == 0);
	}
	return $nreads, $nbases, $nn;
}

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
