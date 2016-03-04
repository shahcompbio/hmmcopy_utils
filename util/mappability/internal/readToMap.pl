#!/usr/bin/perl -w

# ***************************************************************************
# readToMap.pl (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
# Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
# All rights reserved.
# ---------------------------------------------------------------------------
# Last modified: 03 August 2011 (DL)
# ---------------------------------------------------------------------------
# Parses bowtie output and calculates a mappability value for each position
# ***************************************************************************

use strict;
use Getopt::Long;

my $MAXHITS = 4; # Minimum non-zero value == 1 / max_hits
my $PRECISION = 2; # Maximum number of digits after decimal point

GetOptions(
	"max_hits=i" => \$MAXHITS,
	"precision=i" => \$PRECISION
) or die "Unknown option";

$PRECISION += 2;

my $last_chr = "";
my $expected_pos = 0;

while(<STDIN>) {
	# print $_;
	chomp;
	next if (/^$/);
	my ($pos, undef, undef, undef, undef, undef, $hits) = split(/\t/);
	(my $chr, $pos) = split(/:/, $pos);

	if ($last_chr ne $chr) {
		print "fixedStep chrom=$chr start=1 step=1\n";
		$last_chr = $chr;
		$expected_pos = 0;
	}
	while ($pos != $expected_pos) {
		print "0\n";
		++$expected_pos;
	}

	++$hits;
	my $value = 0.0;
	if ($hits <= $MAXHITS) {
		$value = 1.0 / $hits;
	}
	printf("%.$PRECISION" . "s\n", $value);
	++$expected_pos;
}
