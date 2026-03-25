#!/usr/bin/perl

# annotate_vcf.pl

# updated 2022-03-14 to work on input with multiple chromosomes

use warnings;
use strict;
use Getopt::Long;
use IO::Zlib;
use List::MoreUtils qw(uniq);
#use warnings FATAL => 'all';

my $version = '1.3.3'; # fixed identification of sites overlapping deletions
my $alfile = undef;
my $alfields = undef;
my $anctype = 'parsimony';
my $dpbounds = undef;
my $hetbound = undef;
my $vcf = undef;
my $rmfilter = undef;
my $fstart = 4;
my $bed = undef;
my $bedlist = undef;
my $genorep = undef;
my $overwrite;
my $dohet = undef;
my $maxpoly = undef;
my $mingq = 0;
my $devsnp = undef; # ngsParalog calcLR output file
my $devsnp_maxlr = undef; # maximum ngsParalog calcLR LR
my $anno = undef; # file with custom annotations

my %default_flags = (
'REPGQ' => "Number of individuals from different classes that are genotyped with a minimum degree of confidence",
'REPDP' => "Number of genotyped individuals from different classes with a minimum DP",
'HET' => "Proportion of heterozygous individuals",
'DEVLRT' => "Likelihood ratio test statistic from ngsParalog test of a deviant variant",
'LowDP' => "Total site DP below excessively low cutoff of INTEGER",
'HighDP' => "Total site DP above excessively high cutoff of INTEGER",
'NoCov' => "No mapped reads",
'LowCov' => "Excessively low total site depth",
'ExcessCov' => "Excessively high total site depth",
'HighHet' => "P-value for test of excess heterozygosity is below FLOAT",
'LowMQ' => "Low average mapping quality",
'ExcessMQ0' => "Excessively high proportion of reads with mapping quality of zero",
'LowBQ' => "Low average base quality",
'ExcessBQ0' => "Excessively high proportion of reads with base quality of zero",
'Poly' => "Proportion of polymorphic samples above excessively high cutoff of FLOAT",
'DeviantVar' => "Deviant variant based on ngsParalog LRT statistic above cutoff of FLOAT");

die(qq/
annotate_vcf.pl version $version

annotate_vcf.pl [options] <vcf file>
OR
cat <vcf file> | annotate_vcf.pl [options]

--anno       TSV file of annotations for the VCF header, use --anno 'default' to see default descriptions and more info
--alfile     File of outgroup alleles in pafAlleles format
--alfields   Column numbers of alleles to use from alleles file, 1-indexed and ','-separated [$fstart...]
--anctype    Infer ancestral allele using 'parsimony' or as the 'major' allele among outgroups [$anctype]
--dpbounds   LowDP and HighDP bounds to annotate FILTER field, format '<INT lower bound>,<INT upper bound>', triggers LowDP, HighDP flags
--hetbound   Lower p-value bound for excess heterozygosity test \(uses INFO\/ExcHet field\), triggers HighHet flag
--rmfilter   ','-delimited list of FILTER annotations to remove
--bed        Bed format file with FILTER annotations to add in the 4th column, see --flags 'default' for accepted flags
--bedlist    ','-delimited list of FILTER annotations to ignore from bed
--overwrite  Replace all existing filters, otherwise new filters are appended
--genorep    File of sample subset indexes for appending genotype representation to the INFO field, triggers REPGQ, REP flags
--dohet      Calculate proportion of heterozygous genotypes, adds INFO\/HET
--mingq      Minimum genotype quality when counting genotypes [$mingq]
--maxpoly    Upper heterozygosity bound \(uses INFO\/HET field\)
--devlr_file ngsParalog calcLR file of deviant SNP LRT test statistics.
--max_devlr  Maximum ngsParalog calcLR deviant SNP LRT statistic for passing sites.

Notes:
*Assumes the same contig order in VCF and allele files.

*For inferring ancestral alleles with parsimony it is assumed that outgroups in the pafAlleles file 
 are ordered with increasing divergence from the ingroup (greater column number = greater divergence).
 Alternatively, the outgroup columns can be passed in order of their increasing divergence from the
 ingroup using the --alfields option.

*Parsimony currently only implemented for 3 outgroup case.

*Particular FILTER annotations will be updated if the function that applies them (--dpbounds, --bed)
 is called.

*Bed must be sorted in same order as VCF.

*Assumes input VCF indels are left-aligned and normalized (bcftools norm -f <fasta>)

*ngsParalog calcLR input file \(--devlr_file\) must be sorted in the same order as the VCF.
*Assumes ngsParalog file contains a header.

*--genorep file must have columns <group name> <comma-delimited, 0-based indexes of samples in VCF> <min GQ> <min DP>.
 Lines with '#' are ignored (so preceed a header with # for example).
\n/) if (-t STDIN && !@ARGV);

GetOptions('alfile=s' => \$alfile, 'alfields=s' => \$alfields, 'dpbounds=s' => \$dpbounds, 'hetbound=f' => \$hetbound, 'rmfilter=s' => \$rmfilter, 
'anctype=s' => \$anctype, 'bed=s' => \$bed, 'bedlist=s' => \$bedlist, 'overwrite' => \$overwrite, 'genorep=s' => \$genorep, 'dohet' => \$dohet, 
'mingq=i' => \$mingq, 'maxpoly=f' => \$maxpoly, 'devlr_file=s' => \$devsnp, 'max_devlr=f' => \$devsnp_maxlr, 'anno=s' => \$anno);

# this is for printing the command to the VCF header
my %userargs = ('--anno' => $anno, '--alfile' => $alfile, '--alfields' => $alfields, '--dpbounds' => $dpbounds, '--hetbound' => $hetbound, '--rmfilter' => $rmfilter, 
'--anctype' => $anctype, '--bed' => $bed, '--bedlist' => $bedlist, '--overwrite' => $overwrite, '--genorep' => $genorep, '--dohet' => $dohet, '--mingq' => $mingq, 
'--maxpoly' => $maxpoly, '--devlr_file' => $devsnp, '--max_devlr' => $devsnp_maxlr);
# disable printing commands if not applicable
$userargs{"--anctype"} = undef if (!$alfile);
$userargs{"--bedlist"} = undef if (!$bedlist);
$userargs{"--genorep"} = undef if (!$genorep);
$userargs{"--anno"} = undef if (!$anno);

# print information about annotation flags if requested

die(qq/
Optional VCF annotation flags, the default descriptions, and the argument that triggers the flag for inclusion.

To change the default descriptions provide a tab-delimited file to --anno with columns [1] flag ID \(case sensitive\), [2] description \(should not contain tabs\).
Each row in the --anno file is a different flag. Do not enclose descriptions in double quotes.

INFO annotations:
REPGQ  \"$default_flags{REPGQ}\" \(--genorep\)
REPDP  \"$default_flags{REPDP}\" \(--genorep\)
HET    \"$default_flags{HET}\" \(--dohet\)
DEVLRT \"$default_flags{DEVLRT}\" (--devlr_file)

FILTER annotations:
LowDP      \"$default_flags{LowDP}\" \(--dpbounds\)
HighDP     \"$default_flags{HighDP}\" \(--dpbounds\)
HighHet    \"$default_flags{HighHet}\" \(--hetbound\)
Poly       \"$default_flags{Poly}\" \(--maxpoly\)
DeviantVar \"$default_flags{DeviantVar}\" \(--max_devlr\)
NoCov      \"$default_flags{NoCov}\" \(--bed\)
LowCov     \"$default_flags{LowCov}\" \(--bed\)
ExcessCov  \"$default_flags{ExcessCov}\" \(--bed\)
LowMQ      \"$default_flags{LowMQ}\" \(--bed\)
ExcessMQ0  \"$default_flags{ExcessMQ0}\" \(--bed\)
LowBQ      \"$default_flags{LowBQ}\" \(--bed\)
ExcessBQ0  \"$default_flags{ExcessBQ0}\" \(--bed\)
\n/) if ($anno && ! -f $anno && $anno eq "default");

# check arguments
die("Error: --mingq $mingq invalid\n") if ($mingq < 0);

# read in annotation file
if ($anno) {
	open(my $annofh, '<', $anno) or die("Unable to open file of flag descriptions $anno\n");
	while (<$annofh>) {
		chomp;
		my @tok = split(/\t/, $_);
		if (exists $default_flags{$tok[0]}) {
			$default_flags{$tok[0]} = $tok[1];
		} else {
			die("Unrecognized flag $tok[0] in flag description file\n");
		}
	}
	close $annofh;
}

# read in VCF

$vcf = pop @ARGV;
my $vcffh;
if ($vcf) {
	open($vcffh, '<', $vcf) or die("Unable to open VCF $vcf: $!\n");
	sysread $vcffh, my $magic, 2;
	close $vcffh;

	if ($magic eq "\x1f\x8b") {
		# gzipped file
		$vcffh = new IO::Zlib;
		$vcffh->open($vcf, "rb");
	} else {
		open($vcffh, '<', $vcf);
	}
} else {
	die("No input VCF\n") if (-t STDIN);
	$vcffh = \*STDIN;
}


# initialize ancestral allele inference

my ($alfh, $nalleles, $aamethod, @outal, @col, @parstable, %tree_counts);
my ($l1, $l2, $chr, $nsites) = (undef, undef, "", 0);

if ($alfile) {
	open($alfh, '<', $alfile) or die("Could not open allele file $alfile: $!\n");

	if ($alfields) {
		@col = map {die("Argument to --alfields must be ','-delimited integers > 0") if $_ =~ /\D/; $_ - 1} split(/,/, $alfields);
                $nalleles = scalar @col;
	} else {
		$nalleles = split(/\s+/, <$alfh>)-3;
		@col = ($fstart-1) .. ($fstart + $nalleles - 2);
	}

        $l1 = <$alfh>; # first site of new chr
        die("No sites in allele file $alfile\n") if (!$l1);

	if (lc($anctype) eq 'parsimony') {
		$aamethod = 1;
		@parstable = makeLookup(\%tree_counts);
		treeCounter(\%tree_counts, $nalleles, 0);
	} elsif (lc($anctype) eq 'major') {
		$aamethod = 2;
	} else {
		die("Invalid --anctype argument: $anctype\n");
	}
} else {
	$anctype = undef;
}

# initialize ngsParalog annotation

my ($devfh, %devidx);
if ($devsnp) {
	open($devfh, '<', $devsnp) or die("Unable to open --devlr_file $devsnp: $!\n");
	print STDERR "Indexing ngsParalog file\n";
	IndexScaff($devfh, \%devidx, 1); # should make the header a variable because raw calcLR output does not have a header
}

# initialize genotype representation annotations

my (%grp, %grpseen, @grpid, %grpcounts, $grpfh);

if ($genorep) {
	open($grpfh, '<', $genorep) or die("Unable to open sample group index file $genorep: $!\n");
	while(<$grpfh>) {
		if ($_ !~ /^#/) {
			chomp;
			my @tok = split(/\s+/, $_);
			if (! exists $grpseen{$tok[0]}) {
				push @grpid, $tok[0];
				$grpseen{$tok[0]} = 1;
				$grpcounts{$tok[0]}{gq} = 0;
				$grpcounts{$tok[0]}{dp} = 0;
			}
			die("Group info file must have four columns <group ID> <VCF sample indexes> <min GQ> <min DP>: $!\n") if scalar @tok != 4;
			foreach my $idx (split(/,/,$tok[1])) {
				my @idx2 = $idx =~ /(\d+)-(\d+)/ ? ($1..$2) : $idx;
				foreach my $i (@idx2) {
					$grp{$i}{id} = $tok[0];
					$grp{$i}{gq} = $tok[2];
					$grp{$i}{dp} = $tok[3];
				}
			}
		}
	}
close $grpfh;
}

# initialize filter annotation

my (@dpcutoff, %delfilter, $bedfh);

if ($dpbounds) {	
	@dpcutoff = split(/,/, $dpbounds);
	if (scalar(@dpcutoff) != 2) {
		die("ERROR: --dpbounds requires 2 values, a lower and an upper coverage threshold\n");
	} elsif ($dpcutoff[0] >= $dpcutoff[1]) {
		die("ERROR: --dpbounds upper bound must be greater than the lower bound value\n");
	}

	$default_flags{LowDP} =~ s/INTEGER$/$dpcutoff[0]/;
	$default_flags{HighDP} =~ s/INTEGER$/$dpcutoff[1]/;
	
	$delfilter{HighDP} = 1;
	$delfilter{LowDP} = 1;
}

if ($hetbound) {
	die ("--hetbound outside [0,1] interval\n") if ($hetbound < 0 || $hetbound > 1);
	$default_flags{HighHet} =~ s/FLOAT$/$hetbound/;
	$delfilter{HighHet} = 1;
}

if ($devsnp_maxlr) {
	die("--max_devlr must be >= 0\n") if ($devsnp_maxlr < 0);
	$default_flags{DeviantVar} =~ s/FLOAT$/$devsnp_maxlr/;
	$delfilter{DeviantVar} = 1;
}

if (defined $maxpoly && ($maxpoly < 0 || $maxpoly > 1)) {
	die("Error: --maxpoly out of range [0,1]\n");
	$default_flags{Poly} =~ s/FLOAT$/$maxpoly/;
	$delfilter{Poly} = 1;
}

my %bedfilter = (LowCov => 1, ExcessCov => 1, LowMQ => 1, ExcessMQ0 => 1, NoCov => 1, LowBQ => 1, ExcessBQ0 => 1);
my %bedidx;
if ($bed) {
	if ($bedlist) {
		$bedlist .= "," if ($bedlist !~ /,$/);
		$bedlist =~ s/,/\|/g;
	}

	open($bedfh, '<', $bed) or die("Unable to open bed file $bed: $!\n");
	print STDERR "Indexing bed file\n";
	IndexScaff($bedfh, \%bedidx, 0);

	# debug
	#foreach my $k (keys %bedidx) {
	#	seek $bedfh, $bedidx{$k}, 0;
	#	while(<$bedfh>) {
	#		print STDOUT "$k: $_";
	#		last;
	#	}
	#}
	#exit;

	foreach (keys %bedfilter) {
		if ($bedlist && $bedlist =~ /$_\|/) {
			delete $bedfilter{$_};
			next;
		}
		$delfilter{$_} = 1;
	}
}

if ($rmfilter) {
	foreach my $filter_anno (split(/,/,$rmfilter)) {
		$delfilter{$filter_anno} = 1;
	}
}

# print pretty VCF header

my $datestr = localtime();
my $command = "##annotate_vcfCommand=<ID=annotate_vcf.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
	$command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

# filtering parameters.

my ($grp_gq_string, $grp_dp_string);

my $aa_string = "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n";
if ($genorep) {
	$grp_gq_string = "##INFO=<ID=REPGQ,Number=" . scalar(@grpid) . ",Type=Integer,Description=\"$default_flags{REPGQ}\">\n";
	$grp_dp_string = "##INFO=<ID=REPDP,Number=" . scalar(@grpid) . ",Type=Integer,Description=\"$default_flags{REPDP}\">\n";
}
my $het_string="##INFO=<ID=HET,Number=1,Type=Float,Description=\"$default_flags{HET}\">\n";
my $vt_string = "##INFO=<ID=VT,Number=.,Type=String,Description=\"Variant type\">\n";
my $devvar_lrt_string="##INFO=<ID=DEVLRT,Number=1,Type=Float,Description=\"$default_flags{DEVLRT}\">\n";

my $lowdp_string = "##FILTER=<ID=LowDP,Description=\"$default_flags{LowDP}\">\n";
my $highdp_string = "##FILTER=<ID=HighDP,Description=\"$default_flags{HighDP}\">\n";
my $highhet_string = "##FILTER=<ID=HighHet,Description=\"$default_flags{HighHet}\">\n";
my $poly_string = "##FILTER=<ID=Poly,Description=\"$default_flags{Poly}\">\n";
my $lowcov_string = "##FILTER=<ID=LowCov,Description=\"$default_flags{LowCov}\">\n";
my $excesscov_string = "##FILTER=<ID=ExcessCov,Description=\"$default_flags{ExcessCov}\">\n";
my $lowmq_string = "##FILTER=<ID=LowMQ,Description=\"$default_flags{LowMQ}\">\n";
my $mqzero_string = "##FILTER=<ID=ExcessMQ0,Description=\"$default_flags{ExcessMQ0}\">\n";
my $lowbq_string = "##FILTER=<ID=LowBQ,Description=\"$default_flags{LowBQ}\">\n";
my $bqzero_string = "##FILTER=<ID=ExcessBQ0,Description=\"$default_flags{ExcessBQ0}\">\n";
my $nocov_string = "##FILTER=<ID=NoCov,Description=\"$default_flags{NoCov}\">\n";
my $devvar_string = "##FILTER=<ID=DeviantVar,Description=\"$default_flags{DeviantVar}\">\n";

my @headorder = ('fileformat', 'reference', 'contig', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'other'); # header order

my %header = (fileformat => undef, ALT => undef, FILTER => undef, INFO => undef, FORMAT => undef, reference => undef, contig => undef, other => undef);

my %seen = (aa => 0, lowdp => 0, highdp => 0, highhet => 0, het => 0, poly => 0, repgq => 0, repdp => 0, vt => 0, devlrt => 0, devfilter => 0);

my $hline;
my $filter;
while (($hline = <$vcffh>) =~ /^##/) {
	if ($hline =~ /^##([^=]+)/i) {
		my $annotation = $1;

		if (($rmfilter || $bed) && $hline =~ /^##FILTER=<ID=([^,]+)/) {
			next if exists $delfilter{$1}
		}

		if ($alfile && $hline =~ /INFO=<ID=AA/) {
			# update AA annotation to ensure it is correct
			push @{$header{$annotation}}, $aa_string;
			++$seen{aa};
		} elsif ($genorep && $hline =~ /INFO=<ID=REPGQ>/) {
			# update REPGQ annotation
			push @{$header{$annotation}}, $grp_gq_string;
			++$seen{repgq};
		} elsif ($genorep && $hline =~ /INFO=<ID=REPDP>/) {
			# update REPDP annotation
			push @{$header{$annotation}}, $grp_dp_string;
			++$seen{repdp};
		} elsif ($dohet && $hline =~ /INFO=<ID=HET>/) {
			# update HET annotation
			push @{$header{$annotation}}, $het_string;
			++$seen{het};
		} elsif ($devsnp && $hline =~ /INFO=<ID=DEVLRT>/) {
			# update DEVLRT annotation
			push @{$header{$annotation}}, $devvar_lrt_string;
			++$seen{devlrt};
		} elsif ($hline =~ /INFO=<ID=VT>/) {
			# update VT annotation
			push @{$header{$annotation}}, $vt_string;
			++$seen{vt};
		} elsif ($dpbounds && $hline =~ /FILTER=<ID=LowDP/) {
			# update LowDP annotation
			push @{$header{$annotation}}, $lowdp_string;
			++$seen{lowdp};
		} elsif ($dpbounds && $hline =~ /FILTER=<ID=HighDP/) {
			# update HighDP annotation
			push @{$header{$annotation}}, $highdp_string;
			++$seen{highdp};
		} elsif ($hetbound && $hline =~ /FILTER=<ID=HighHet>/) {
			# update HighHet annotation
			push @{$header{$annotation}}, $highhet_string;
			++$seen{highhet};
		} elsif ($devsnp_maxlr && $hline =~ /FILTER=<ID=DeviantVar>/) {
			# update DeviantVar annotation
			push @{$header{$annotation}}, $devvar_string;
			++$seen{devfilter};
		} elsif (defined $maxpoly && $hline =~ /FILTER=<ID=Poly>/) {
			# update max number of polymorphic samples string
			push @{$header{$annotation}}, $poly_string;
			++$seen{poly}
		} elsif ($annotation =~ /^fileformat|ALT|FILTER|INFO|FORMAT|reference|contig/) {
			push @{$header{$annotation}}, $hline;
		} else {
			push @{$header{other}}, $hline;
		}
	}
}

push @{$header{INFO}}, $aa_string if ($alfile && !$seen{aa});
push @{$header{INFO}}, $grp_gq_string if ($genorep && !$seen{repgq});
push @{$header{INFO}}, $grp_dp_string if ($genorep && !$seen{repdp});
push @{$header{INFO}}, $het_string if ($dohet && !$seen{het});
push @{$header{INFO}}, $vt_string if (!$seen{vt});
push @{$header{INFO}}, $devvar_lrt_string if ($devsnp && !$seen{devlrt});
push @{$header{FILTER}}, $lowdp_string if ($dpbounds && !$seen{lowdp});
push @{$header{FILTER}}, $highdp_string if ($dpbounds && !$seen{highdp});
push @{$header{FILTER}}, $highhet_string if ($hetbound && !$seen{highhet});
push @{$header{FILTER}}, $devvar_string if ($devsnp_maxlr && !$seen{devfilter});
push @{$header{FILTER}}, $poly_string if (defined $maxpoly && !$seen{poly});
push @{$header{FILTER}}, $nocov_string if ($bedfilter{NoCov});
push @{$header{FILTER}}, $lowcov_string if ($bedfilter{LowCov});
push @{$header{FILTER}}, $excesscov_string if ($bedfilter{ExcessCov});
push @{$header{FILTER}}, $lowmq_string if ($bedfilter{LowMQ});
push @{$header{FILTER}}, $mqzero_string if ($bedfilter{ExcessMQ0});
push @{$header{FILTER}}, $lowbq_string if ($bedfilter{LowBQ});
push @{$header{FILTER}}, $bqzero_string if ($bedfilter{ExcessBQ0});
push @{$header{other}}, $command;


# make matrix for quick genorep logging
my @grpmat;
if ($genorep) {
	my $fieldidx = 0;
	foreach (split(/\s+/, $hline)) {
		push @grpmat, [0,0,0,0];
		if ($fieldidx > 8) {
			my $sampidx = $fieldidx - 9;
			if (exists $grp{$sampidx}) {
				$grpmat[$fieldidx][0] = 1;
				$grpmat[$fieldidx][1] = $grp{$sampidx}{id};
				$grpmat[$fieldidx][2] = $grp{$sampidx}{gq};
				$grpmat[$fieldidx][3] = $grp{$sampidx}{dp};
			}
		}
		$fieldidx++;
	}
}

# check that info required for deviant SNP filtering is supplied
if ($devsnp_maxlr) {
	die("Deviant variant flagging with --max_devlr requires that INFO/DEVLRT exists or a file with ngsParalog LRT statistics is supplied with --devlr_file\n") unless ($seen{devlrt} || $devsnp);
}

# print header info

foreach my $tag (@headorder) {
	foreach (@{$header{$tag}}) {
		print STDOUT $_;
	}
}

print STDOUT $hline;

# process VCF sites
print STDERR "Annotating VCF\n";

my %filtercounts;
map { $filtercounts{$1} = 0 if $_ =~ /ID=([^,]+)/ } @{$header{FILTER}};
my %filterconfigs;

my $vcfchr;
my @bedtok = ('',0,0);
my $site_count = 0;
my $fail_counts = 0;
my $aa_count = 0;
my @tok;
my %missing_chr;
my %bedlast;
my @indel = (0,0,''); # start and stop coordinates for indel
my @devtok = ('',0,0);
my %devlast;

# variables to hold variant counts
my $nsnps = 0; # total number SNPs
my $nsnp_sites; # total number sites containing SNPs
my $nsnp_widel = 0; # number of SNP sites within a deletion
my $nbisnps = 0; # number sites with biallelic snps
my $nbisnps_widel = 0; # number sites with biallelic snps within a deletion
my $ninsert = 0; # number of insertions
my $ndelete = 0; # number of deletions

$" = "\t";

while (<$vcffh>) {
	@tok = split(/\s+/, $_);
	$vcfchr = uc($tok[0]);
	$site_count++;

	## Determine variant type
	@indel = (0,0,'') if ($vcfchr ne $indel[2]); # reset indel for new chromosome
	die("Multiple reference alleles found at $tok[0] $tok[1]") if ($tok[3] =~ /,/);
	die("Unexpecte allele in reference at $tok[0] $tok[1]") if ($tok[4] =~ /[^ACGT]i/);
	my $reflen = length($tok[3]);
	if ($reflen > 1) {
	# left-alignment means this is a deletion wrt reference
		# figure out if site is within a new deletion
		my $varend = $tok[1] + ($reflen-1);
		if ($varend > $indel[1]) {
			my $start = $tok[1]+1;
			$indel[0] = $tok[1]+1 if $start > $indel[1]; # this is the first base actually deleted
			$indel[1] = $varend; # this is the last base actually deleted
			$indel[2] = uc($tok[0]); # chromosome
		}
	}
	my $widel = ($tok[1] >= $indel[0] && $tok[1] <= $indel[1] && $vcfchr eq $indel[2]) ? 1 : 0;
	
	my ($issnp, $isins, $isdel) = (0,0,0);
	my @snp_alleles = ();
	if ($tok[4] ne '.') {
		my @altarr = split(/,/, $tok[4]);
		my $altn;
		my $refbase = substr $tok[3], 0, 1;
		my $var_count = 1;

		foreach my $a (@altarr) {
			if ($a eq '*') {
				# do not record star alleles
				$var_count++;
				next;
			}
			my $altlen = length($a);
			if ($altlen == $reflen) {
				# SNP
				$issnp = 1;
				my $altbase = substr $a, 0, 1;
				if ($altbase ne $refbase) {
					$altn = $var_count;
					push @snp_alleles, $altbase;
					$nsnps++;
				}
			} else {
				# indel
				if ($altlen > $reflen) {
					$isins = 1;
					$ninsert++;
				} else {
					$isdel = 1;
					$ndelete++;
				}
			}
			$var_count++;
		}

	}

	$nsnp_sites++ if ($issnp);
	$nsnp_widel++ if ($issnp && $widel);
	if (scalar(@snp_alleles) == 1) {
		$nbisnps++;
		$nbisnps_widel++ if ($widel);
	}

	my $vt = '';
	if ($tok[7] =~ /VT=([^\s|;|=]+)/) {
		$vt = $1;
		my $vt_original = $vt;
		$vt .= ',SNP' if ($issnp && $vt !~ /snp/i);
		$vt .= ',DEL' if ($isdel && $vt !~ /deletion/i);
		$vt .= ',INS' if ($isins && $vt !~ /insertion/i);
		$vt .= ',OVLDEL' if ($widel && $issnp && $vt !~ /OVLDEL/i); # SNP Overlaps Deletion
		if ($tok[4] eq ".") {
			# monomorphic site
			$vt = $widel ? 'OVLDEL' : '.'; # site Overlaps Deletion
		}
		$tok[7] =~ s/VT=$vt_original/VT=$vt/;
	} else {
		$vt .= ',SNP' if ($issnp);
		$vt .= ',DEL' if ($isdel);
		$vt .= ',INS' if ($isins);
		$vt .= ',OVLDEL' if ($widel && $issnp);
		$vt =~ s/^,//;
		if ($tok[4] eq ".") {
			# monomorphic site
			$vt = $widel ? 'OVLDEL' : '.';
		}
		$tok[7] .= ";VT=$vt";
	}


	## update bedgraph region - assumes VCF and bed are sorted in the same order
	if ($bed && !exists $bedlast{$vcfchr}) {
		# update chromosome
		if ($bedtok[0] ne $vcfchr && exists $bedidx{$vcfchr}) {
			seek $bedfh, $bedidx{$vcfchr}, 0;
			bedsplit(\@bedtok, $bedfh);
		}
		
		#update position
		while ($bedtok[0] eq $vcfchr && $bedtok[2] < $tok[1] && !eof($bedfh)) { 
			bedsplit(\@bedtok, $bedfh);
			if ($bedtok[0] ne $vcfchr) {
				$bedlast{$vcfchr} = 1;
				last;
			}
		}
	}

	## update deviant SNP file site - assumes VCF and ngsParalog calcLR file are sorted in the same order
	if ($devsnp && !exists $devlast{$vcfchr}) {
		# update chromosome
		if ($devtok[0] ne $vcfchr && exists $devidx{$vcfchr}) {
			seek $devfh, $devidx{$vcfchr}, 0;
			bedsplit(\@devtok, $devfh);
		}
			
		# update position
		while ($devtok[0] eq $vcfchr && $devtok[1] < $tok[1] && !eof($devfh)) {
			bedsplit(\@devtok, $devfh);
			if ($devtok[0] ne $vcfchr) {
				$devlast{$vcfchr} = 1;
				last;
			}
		}
	}
	
	## ancestral allele
	if ($alfile) {
		if ($chr ne $tok[0]) {
		# store ancestral alleles in massive matrix
			@outal = ();
			$chr = $1 if ($l1 =~ /^(\S+)/);

			while ($chr ne $tok[0]) {
				my $l1 = <$alfh>;
				die("$tok[0] not found in allele file $alfile\n") if (!$l1);
				$chr = $1 if ($l1 =~ /^(\S+)/);
			}
			$l2 = $l1;
			my $c;

			do {
				push @outal, [(split(/\s+/, $l2))[@col]];
				if (defined($l2 = <$alfh>)) {
					$c = $1 if ($l2 =~ /^(\S+)/);
				}
			} while (defined($l2) && $c eq $chr);

			$l1 = $l2;
			$nsites = scalar(@outal);

			# debug
			#print "$nsites\n";
			#foreach (@outal) {
			#	print "@$_\n";
			#}
			#exit;
		}

		# collect ingroup alleles
		die("More than one reference allele at $tok[0] $tok[1]: $tok[3]\n") if ($tok[3] =~ /,/);
		my $reflen = length($tok[3]);
		my @ingroup = ($tok[3], split(/,/, $tok[4]));

		# collect alleles for each outgroup over the length of the ingroup ref allele
		my @outgroup;
		my $start = $tok[1]-1;
		my $end = $start + $reflen;
		# loop over each outgroup assembly
		for (my $assem=0; $assem<$nalleles; $assem++) {

			if (length($outal[$start]->[$assem]) > $reflen) {
				# insert wrt to ingroup reference
				push @outgroup, $outal[$start]->[$assem];
			} else {
				# non-insertion wrt to ingroup reference
				$outgroup[$assem] = "";
				for (my $s=$start; $s<$end; $s++) {
					$outgroup[$assem] .= $outal[$s]->[$assem] if ($outal[$s]->[$assem] ne '*');
				}

				if (!$outgroup[$assem]) {
					# check if all aligned sites were deleted
					$outgroup[$assem] = '*';
				} else {
					# check for deletions longer that ingroup deletions
					if ($outal[$end-1]->[$assem] eq '*' && $end <= $#outal && $outal[$end]->[$assem] eq '*') {
						$outgroup[$assem] .= '*'; # denotes that next site beyond last reference allele site is deleted
					}
				}
			}
		}

		# find ancestral allele
		my $anc_allele;

		if ($aamethod == 1) {
			# use parsimony to find ancestral allele

			# Find the outgroup matching configuration for each ingroup allele in the lookup table
			my @outsp;
			my $tree;

			foreach my $inallele (@ingroup) {

				my $tabref = \@parstable;
				$tree = "";

				foreach my $outallele (@outgroup) {
					if ($outallele =~ /N/ ) {
						# outgroup allele is missing
						$tabref = ${$tabref}[2];
						$tree .= 2;
					} else {
						# check if ingroup allele matches outgroup
						my $match = $inallele eq $outallele ? 1 : 0;
						$tabref = ${$tabref}[$match];
						$tree .= $match;
					}

				}

				push @outsp, $tabref;
				$tree_counts{$tree}++;
				$tree .= ',';
			}

			# resolve exceptional cases for most parsimonious solution
			if ($nalleles == 3) {
				# 3 outgroup case
				if ($tree =~ /100/ && $tree !~ /011/) {
					push @outsp, 0;
				}
			}

			# pick ancestral allele from the most related outgroup species
			my $spid = (sort {$a <=> $b} @outsp)[0];
			$anc_allele = $spid != 9 ? $outgroup[$spid] : 'N';

			# debug
			#print "$tok[0]\t$tok[1]\t@ingroup;\t@outgroup;\t@outsp;\t$anc_allele\n"; next;			

		} elsif ($aamethod == 2) {
			# use major allele among outgroups as ancestral allele

			my %allele_counts = ();
			foreach my $inallele (@ingroup) {
				$allele_counts{$inallele} = 0;
			}

			# count occurances of outgroup alleles among the ingroup
			my $aa_seen = 0; # records whether at least one ref or alt allele is present among outgroups
			foreach my $outallele (@outgroup) {
				if (exists $allele_counts{$outallele}) {
					$allele_counts{$outallele}++;
					$aa_seen = 1;				
				}
			}

			if ($aa_seen) {
				my $maxcount = 0;
				my $nmax = 0;
					foreach my $allele (@ingroup) {
					if ($allele_counts{$allele} >= $maxcount) {
						$nmax = $allele_counts{$allele} == $maxcount ? $nmax + 1:1;
						$maxcount = $allele_counts{$allele};
						$anc_allele = $allele;
					}
				}
				$anc_allele = "N" if ($nmax > 1 || $maxcount < 1);
		
			} else {
				$anc_allele="N";
			}
		}

		# add ancestral allele to INFO
		if (!($tok[7] =~ s/AA=[^;|\s]+/AA=$anc_allele/)) {
			$tok[7] .= ";AA=$anc_allele";
		}
		$aa_count++ if ($anc_allele ne 'N');
	}

	## find FORMAT subfield indices
	my ($gtidx, $dpidx, $gqidx, $iter) = (-9, -9, -9, 0);
	if ($genorep || $dohet) {
		foreach my $fmtid (split(':', $tok[8])) {
			if ($fmtid eq 'GT') {
				$gtidx = $iter;
			} elsif ($fmtid eq 'GQ') {
				$gqidx = $iter;
			} elsif ($fmtid eq 'DP') {
				$dpidx = $iter;
			}
			last if ($gtidx > -9 && $gqidx > -9 && $dpidx > -9);
			$iter++;
		}
		die("FMT/GT missing at $tok[0] $tok[1]\n") if ($dohet && $gtidx == -9);
		die("FMT/GQ missing at $tok[0] $tok[1]\n") if ($dohet && $mingq > 0 && $tok[4] ne "." && $gqidx == -9);
		die ("FMT/DP missing at $tok[0] $tok[1]\n") if ($genorep && $dpidx == -9);
	}

	## genotype respresentation info
	if ($genorep) {
		# reset counts
		foreach my $id (@grpid) {
			$grpcounts{$id}{gq} = 0;
			$grpcounts{$id}{dp} = 0;
		}

		# count genotypes
		for (my $i = 9; $i <= $#tok; $i++) {
			if ($grpmat[$i][0]) {
				my @geno = split(':', $tok[$i]);
				$grpcounts{$grpmat[$i][1]}{dp}++ if ($dpidx > -9 && $geno[0] !~ /\./ && $geno[$dpidx] >= $grpmat[$i][3]);
				$grpcounts{$grpmat[$i][1]}{gq}++ if ($gqidx > -9 && $geno[0] !~ /\./ && $geno[$gqidx] >= $grpmat[$i][2]);
			}
		}

		# record genotype representation info
		my $repdp_str = "REPDP=";
		my $repgq_str = "REPGQ=";
		$iter = 0;
		foreach my $id (@grpid) {
			$repdp_str.="$grpcounts{$id}{dp}" if ($dpidx > -9);
			$repgq_str.="$grpcounts{$id}{gq}" if ($gqidx > -9);
			if ($iter < $#grpid) {
				$repdp_str .= ",";
				$repgq_str .= ",";
			}
			$iter++;
		}
		$repdp_str = "" if ($dpidx == -9);
		$repgq_str = "" if ($gqidx == -9);
		
		addInfo(\$tok[7], "REPDP=", \$repdp_str);
		addInfo(\$tok[7], "REPGQ=", \$repgq_str);
		$tok[7] =~ s/;{2,}/;/g;
		$tok[7] =~ s/;$//;
	}

	## heterozygote frequency
	if ($dohet && $tok[4] ne ".") {
		my ($ngt, $nhet) = (0,0);
		for (my $i = 9; $i <= $#tok; $i++) {
			my @indarr = split(':',$tok[$i]);
			if ($indarr[$gtidx] =~ /\./) {
				# missing genotype
				next;
			} else {
				if ($mingq == 0 || ($gqidx != -9 && $indarr[$gqidx] >= $mingq)) {
					$ngt++;
					my @a = split(/[\/\|]/,$indarr[$gtidx]);
					$nhet++ if (keys %{{ map{$_,1} @a }} > 1);
				}
			}

		}
		my $h = ($ngt == 0) ? '.' : $nhet/$ngt;
		$h = sprintf("%.7f",$h) if ($h ne '.' && $h != 0);
		addInfo(\$tok[7], "HET=", \"HET=${h}");
	}

	## ngsParalog INFO annotation
	if ($devsnp) {
		my $deviant_lrt = '.';
		$deviant_lrt = $devtok[4] if ($devtok[0] eq $vcfchr && $devtok[1] eq $tok[1]);
		addInfo(\$tok[7],"DEVLRT=",\"DEVLRT=${deviant_lrt}");
	}

	## remove filter annotations
	if ($rmfilter) {
		if ($tok[6] ne '.') {
			my $filter_str = "";
			foreach my $fid (split(/;/, $tok[6])) {
				$filter_str .= "$fid;" unless exists $delfilter{$fid};
			}
			$filter_str =~ s/;*$//;
			$tok[6] = (!$filter_str) ? '.' : $filter_str; 
		}
	}

	## apply filters
	if ($dpbounds || $hetbound || $bed) {
		my @filter_annotations;
	
		# coverage filters
		if ($dpbounds) {
			my $dp = $1 if ($tok[7] =~ /DP=(\d+)/);
			if (defined $dp) {
				if ($dp < $dpcutoff[0]) {
					push @filter_annotations, "LowDP";
				} elsif ($dp > $dpcutoff[1]) {
					push @filter_annotations, "HighDP";
				}
			} else {
				print STDERR "No DP info for $tok[0] $tok[1], INFO=<$tok[7]>\n";
			}
		}

		# excess heterozygosity filter
		# Flags a site based on the lowest p-value for multiallelic sites
		if ($hetbound && $tok[7] =~ /ExcHet=([^;]+)/) {
			foreach my $p (split(/,/, $1)) {
				push @filter_annotations, "HighHet" if ($p < $hetbound);
			}
		}

		# excess polymorphic (i.e. heterozygous) sample filter based on just the number of polymorphic samples
		# This is useful for flagging polymorphic clones.
		if (defined $maxpoly && $tok[4] ne '.') {
			if ($tok[7] =~ /HET=([^;]+)/) {
				push @filter_annotations, "Poly" if ($1 > $maxpoly);
			} else {
				print STDERR "Warning: INFO/HET not found at $tok[0] $tok[1] for use with --maxpoly. Try running with --dohet.\n";
			}
		}

		# ngsParalog deviant variant filter
		if (defined $devsnp_maxlr && $tok[7] =~ /DEVLRT=([^;]+)/) {
			my $lrtval=$1;
			push @filter_annotations, "DeviantVar" if ($lrtval =~ /^\d/ && $lrtval > $devsnp_maxlr);
		}

		# bed filters
		if ($bed && $vcfchr eq $bedtok[0] && $tok[1] > $bedtok[1] && $tok[1] <= $bedtok[2]) {
			map {push @filter_annotations, $_ if exists $bedfilter{$_}} split(';',$bedtok[3]); 
		}

		# update FILTER field with new annotations
		if (@filter_annotations) {
			if ($tok[6] eq "." || $tok[6] eq "PASS") {
				$tok[6] = join(';',@filter_annotations);
			} else {
				if ($overwrite) {
					$tok[6] .= join(';',@filter_annotations);
				} else {
					push @filter_annotations, split(';', $tok[6]);
					$tok[6] = join(';', uniq @filter_annotations);
				}
			}
		} else {
			$tok[6] = "PASS" if ($tok[6] eq '.');
		}
	}
	
	if ($tok[6] ne '.' && uc($tok[6]) ne 'PASS') {
		$fail_counts++;
		map {$filtercounts{$_}++} split(';', $tok[6]);
	}
	if (exists $filterconfigs{$tok[6]}) {
		$filterconfigs{$tok[6]}++;
	} else {
		$filterconfigs{$tok[6]} = 1;
	}

	print STDOUT "@tok\n";
}

my $nvariants = $nsnps + $ndelete + $ninsert; # total number of variants

close $vcffh;
close $alfh if ($alfile);
close $bedfh if ($bed);
close $devfh if ($devfh);

# print counts
print STDERR "\n===== REPORT =====\n";
print STDERR "Last contig processed: $tok[0]\n";
print STDERR "Number sites in VCF: $site_count\n";

print STDERR "\n===== VARIANT COUNTS =====\n";
print STDERR "Number variants: $nvariants\n";
print STDERR "Number SNPs: $nsnps\n";
print STDERR "Number insertions: $ninsert\n";
print STDERR "Number deletions: $ndelete\n";
print STDERR "Number SNP sites: $nsnp_sites\n";
print STDERR "Number bialleleic SNP sites: $nbisnps\n";
print STDERR "Number SNP sites within a deletion: $nsnp_widel\n";
print STDERR "Number biallelic SNP sites within a deletion: $nbisnps_widel\n";

if ($alfile) {
	print STDERR "\n===== ANCESTRAL ALLELE REPORT =====\n";
	print STDERR "Number of sites with ancestral alleles: $aa_count\n";
	if ($aamethod == 1) {
		print STDERR "Ingroup allele matching configurations\n";
		treeCounter(\%tree_counts, $nalleles, 1);
	}
}

print STDERR "\n===== FILTER SUMMARY =====\n";
print STDERR "Number of FAIL sites: $fail_counts\n";
print STDERR "\n===== FILTER CONFIGURATION COUNTS =====\n";
foreach (keys %filterconfigs) {
	print STDERR "$_: $filterconfigs{$_}\n" unless ($_ eq 'PASS' || $_ eq '.');
}

print STDERR "\n===== FILTER FLAG COUNTS =====\n";
foreach (@{$header{FILTER}}) {
	if ($_ =~ /ID=([^,]+)/) {
		print STDERR "$1: $filtercounts{$1}\n" unless ($1 eq 'PASS' || $1 eq '.');
	}
}

exit;

sub bedsplit {
	my ($tok, $fh) = @_;
	@$tok = split(/\s+/, readline($$fh));
	$$tok[0] = uc($$tok[0]);
}

sub IndexScaff {
	# indexes scaffold position in any file where the first field is scaffold ID
	my ($fh, $scaffhash, $header) = @_;
	my $idx = 0;
	chomp(my $scaff = uc((split(/\s+/,readline($$fh)))[0]));
	$$scaffhash{$scaff} = $idx; # this might be the file header if there is one, but that's okay
	my $lastscaff = $scaff;
	$idx = tell $$fh;
	while(my $line = readline($$fh)) {
		chomp($line);
		my @tok = split(/\s+/, $line);
		$scaff = uc($tok[0]);
		$$scaffhash{$scaff} = $idx if ($scaff ne $lastscaff && !exists $$scaffhash{$scaff});
		$idx = tell $$fh;
		$lastscaff = $scaff;
	}
	if ($header) {
		my @scaffkeys = grep {$$scaffhash{$_} == 0} keys %{$scaffhash};
		delete $$scaffhash{$scaffkeys[0]};
	}
}


sub makeLookup {
### returns species representing the ancestral allele ###

### species in descending relatedness to Astatotilapia calliptera
# Simochromis diagramma, Tropheini tribe (0)
# Neolamprologus brichardi, Lamprologini tribe (1)
# Oreochromis niloticus (2)
# not able to assign (9)

### codes corresponding to elements of of the lookup table
# [S. diagramma]->[N. brichardi]->[O. niloticus]
# 0: mismatch to A. calliptera allele
# 1: match A. calliptera allele
# 2: missing data

        my @anc;

        $anc[0]->[0]->[0] = 9;
        $anc[0]->[0]->[1] = 9;
        $anc[0]->[0]->[2] = 9;
        $anc[0]->[1]->[0] = 9;
        $anc[0]->[1]->[1] = 1;
        $anc[0]->[1]->[2] = 1;
        $anc[0]->[2]->[0] = 9;
        $anc[0]->[2]->[1] = 2;
        $anc[0]->[2]->[2] = 9;
        $anc[1]->[0]->[0] = 9;
        $anc[1]->[0]->[1] = 0;
        $anc[1]->[0]->[2] = 0;
        $anc[1]->[1]->[0] = 0;
        $anc[1]->[1]->[1] = 0;
        $anc[1]->[1]->[2] = 0;
        $anc[1]->[2]->[0] = 0;
        $anc[1]->[2]->[1] = 0;
        $anc[1]->[2]->[2] = 0;
        $anc[2]->[0]->[0] = 9;
        $anc[2]->[0]->[1] = 2;
        $anc[2]->[0]->[2] = 9;
        $anc[2]->[1]->[0] = 1;
        $anc[2]->[1]->[1] = 1;
        $anc[2]->[1]->[2] = 1;
        $anc[2]->[2]->[0] = 9;
        $anc[2]->[2]->[1] = 2;
        $anc[2]->[2]->[2] = 9;

        if (@_) {
                my $config_counts = $_[0];
                for (my $i=0; $i<3; $i++) {
                        for (my $j=0; $j<3; $j++) {
                                for (my $k=0; $k<3; $k++) {
                                        $$config_counts{"$i$j$k"} = 0;
                                }
                        }
                }
        }

	return @anc;
}

sub treeCounter {
# fun = 0: initalize count hash
# fun = 1: print count hash

	my ($counts, $nout, $fun) = @_;

	# generate array of configurations
	my @config = ("");

	for (my $i=0; $i<$nout; $i++) {
		my @arr;
		foreach my $str (@config) {
			foreach my $state ((0,1,2)) {
				push @arr, "$str$state";
			}
		}
		@config = @arr;
	}

	# process configurations
	foreach (@config) {
		if ($fun == 0) {
			$$counts{$_} = 0;
		} elsif ($fun == 1) {
			print STDERR "$_: $$counts{$_}\n";
		}
	}
}

sub addInfo {
	my ($info, $oldstr, $newstr) = @_;
	if ($$info =~ /$oldstr/) {
		$$info =~ s/${oldstr}[^;]+/$$newstr/;
	} else {
		$$info .= ";$$newstr";
	}
}
