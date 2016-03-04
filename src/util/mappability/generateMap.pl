#!/usr/bin/perl -w

# ***************************************************************************
# generateMap.pl (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
# Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
# All rights reserved.
# ---------------------------------------------------------------------------
# Last modified: 03 August 2011 (DL)
# ---------------------------------------------------------------------------
#
# Takes in a FASTA reference and outputs and BigWig file of mappability values
#
# In reality just a wrapper script around a four program pipe, namely:
#     fastaToRead | bowtie | readToMap | wigToBigWig
#
# If more fine control is required in the process, users can call the individual
# programs manually.  In practice with a 3.06 GHz CPU and 8GB RAM, E.coli:O157
# (NC_002655.2, 5.5MB) could be processed in about 2 minutes, and C. elegans
# (ce6, 100MB) in 40 minutes on a single core, with all default settings (35-mers).
#
# Hidden options:
#     -l, --list              List all chromosomes in FASTA reference file
#     -c, --chr <string>      Specify the entries and order of sequences to analyze [ALL]
#                             the <string> should be a comma-delimited list (NO spaces)
#                             NOTE: BigWig files are only valid with all chromsomes,
#                             Only reason to use this would be to create WIG files for each
#                             chromosome, cat all resultant files, then run wigToBigWig
#
#     -m, --maxhits <int>     Changes the minimum mappability value, where
#                             minimum non-zero mappability == 1 / maxhits
#     -u, --uncompress        Outputs WIG to stdout instead of BigWig to file,
#                             for debugging, or manually running wigToBigWig.
#     -r, --rename            Renames chromosomes to be UCSC Genome Browser/IGV compliant,
#                             only applicable to human.
# ***************************************************************************

use strict;
use Getopt::Long;
use File::Basename;

sub usage() {
	print "Usage: $0 [options] <FASTA reference>";
	print "\n";
	print "Options:\n";
	print "    -o, --output <string>   Output (also takes 'stdout') [default: <FASTA reference>.map.bw]\n";
	print "    -w, --window <int>      Specify the fragment size to calculate mappability values [35]\n";
	print "\n";
	print "    -i, --index <string>    Location of a ready built bowtie index of the FASTA input\n";
	print "    -b, --build             Build bowtie index for given FASTA reference, take caution\n";
	print "                            when using on large genomes, 'bowtie-build' for details\n";
	print "\n";
	print "    -h, --help              Prints this message\n";
	print "\n";
	print "Example:\n";
	print "    $0 hg18.fa\n";
	print "\n";
	print "Author: Daniel Lai <jujubix\@cs.ubc.ca>\n";
}

my $window = 35;
my $list;
my $build;
my $sequence = "";
my $index;
my $maxhits = 4;
my $output;
my $uncompress;
my $rename;

GetOptions(
	"rename" => \$rename,
	"output=s" => \$output,
	"uncompress" => \$uncompress,
	"maxhits=i" => \$maxhits,
	"window=i" => \$window,
	"index=s" => \$index,
	"list" => \$list,
	"build" => \$build,
	"chr=s" => \$sequence,
	"help" => sub{ usage(); exit(0); }
) or die "Unknown option or missing option argument, -h for help\n";

my ($file, $dir) = fileparse ($0);
my $fastaToRead = "$dir" . "internal/fastaToRead";
my $readToMap = "$dir" . "internal/readToMap.pl";
my $renameChr = "$dir" . "../renameChr.pl";
my $wigToBigWig = "$dir" . "../bigwig/wigToBigWig";

unless(-e "$fastaToRead") {
	die "$fastaToRead missing, try recompiling...\n";
}

unless(-e "$readToMap") {
	die "readToMap missing, try recompiling...\n";
}

unless(-e "$renameChr") {
	die "$renameChr missing, try recompiling...\n";
}

unless(-e "$wigToBigWig") {
	die "$wigToBigWig missing, try recompiling...\n";
}

my $status = `bowtie --version`;
if (!defined($status)) {
	die "bowtie not found, obtain from http://bowtie-bio.sourceforge.net/index.shtml and/or add to \$PATH\n";
}

if (!defined($ARGV[0])) {
	warn "Please provide a FASTA reference file\n";
	usage();
	exit(0);
}

my $fasta = $ARGV[0];
if (!defined($index)) {
	$index = $fasta;
}

my $sizes = "$fasta.sizes";
unless (-e $sizes) {
	system("$fastaToRead -l $fasta > $fasta\.sizes");
}

if ($list) {
	system("$fastaToRead -l $fasta");
	exit(0);
}

if ($build) {
	system("bowtie-build $fasta $fasta");
	exit(0);
}

if ($sequence ne "") {
	$sequence = "-c $sequence";
}

if (!defined($output)) {
	$output = "$fasta.map.bw";
	warn "No output file specified, writing BigWig file to: $output\n";
}

$rename = ($rename) ? "| $renameChr " : "";
my $first = "$fastaToRead -w $window $sequence $fasta";
my $second = "bowtie $index -f - -v 0 --quiet";
my $third = "$readToMap -m $maxhits $rename";
my $fourth = "$wigToBigWig stdin $sizes $output";

if ($uncompress) {
	$status = system("$first | $second | $third");
} else {
	$status = system("$first | $second | $third | $fourth");
}

exit($status);
