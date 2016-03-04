#!/usr/bin/perl -w

# ***************************************************************************
# renameChr.pl (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
# Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
# All rights reserved.
# ---------------------------------------------------------------------------
# Last modified: 05 August 2011 (DL)
# ---------------------------------------------------------------------------
# Filter to rename chromosome names into UCSC/IGV compliant names for hg18/19
# ***************************************************************************

use strict;
use Getopt::Long;

# Hardcoded mapping of bowtie to UCSC chromosomes names
my %mapping = (
	"1" => "chr1",
	"2" => "chr2",
	"3" => "chr3",
	"4" => "chr4",
	"5" => "chr5",
	"6" => "chr6",
	"7" => "chr7",
	"8" => "chr8",
	"9" => "chr9",
	"10" => "chr10",
	"11" => "chr11",
	"12" => "chr12",
	"13" => "chr13",
	"14" => "chr14",
	"15" => "chr15",
	"16" => "chr16",
	"17" => "chr17",
	"18" => "chr18",
	"19" => "chr19",
	"20" => "chr20",
	"21" => "chr21",
	"22" => "chr22",
	"23" => "chrX",
	"24" => "chrY",
	"X" => "chrX",
	"Y" => "chrY",
	"M" => "chrM",
	"MT" => "chrM"
);

my $sizes = 0;

GetOptions(
	"sizes" => \$sizes,
) or die "Unknown option, only valid open: -s for processing .sizes files\n";

while(<STDIN>) {
	if ($sizes) {
		if ($_ =~ /(\S+)\t(\d+)/) {
			my $chr = $mapping{$1};
			if (defined($chr)) {
				$_ =~ s/$1/$chr/;
				print;
			} else {
				print;
			}
		} else {
			die "Invalid file format\n";
		}
	} else {
		if ($_ =~ /^fixedStep chrom=(\S+)/) {
			my $chr = $mapping{$1};
			if (defined($chr)) {
				$_ =~ s/chrom=$1/chrom=$chr/;
				print;
			} else {
				print;
			}
		} else {
			print;
		}
	}
}
