# HMM Copy Utils

Author: Daniel Lai <dalai@bccrc.ca> Department of Molecular Oncology, BC Cancer Research Agency
Date:   August 02, 2011

## Overview

This repo provides tools for extracting read counts and gc and mappability statistics in preparation for running HMMCopy.

## Installation

### Using CMake

Requires cmake installed and on PATH (http://www.cmake.org/).
In the source directory, install as follows:

	cmake .
	make

Should result in main binaries in /bin and useful binaries in /util

## Usage

### mapCounter

*fast* average mappability counter using BigWig files

Overview:

mapCounter is a small program for calculating average [mappability] for non-overlapping
windows of fixed width across all sequences (i.e. chromosomes) present in a [BigWig file].
It is built mainly on top of an independent subset of files obtained from the
[UCSC Genome Browser source code] (i.e. kent library) made available by Jim Kent.
Generating average mappability files in 1000 base windows (default) on hg18 took 160
seconds on a 3.06 GHz Intel Core 2 Duo with 8GB RAM.

[mappability]: http://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeMapability
[BigWig file]: http://genome.ucsc.edu/goldenPath/help/bigWig.html
[UCSC Genome Browser source code]: http://genome.ucsc.edu/admin/git.html

Features:

- Average mappability calculation for fixed non-overlapping windows

CLI:

	Usage: ./mapCounter [options] <BigWig file>

	Options:
		-w, --window <int>       Specify the size of non-overlapping windows [1000]
		-l, --list               List all chromosomes in BigWig file
		-s, --sequence <string>  Specify the entries and order of sequences to analyze [ALL],
								 the <string> should be a comma-delimited list (NO spaces)
	Example:
		./mapCounter -w 100000 -s 1,3,5,X hg18.bw > hg18.map.seg

### gcCounter

Overview:

gcCounter is a small program for calculating GC content for non-overlapping
windows of fixed width across all sequences (i.e. chromosomes) present in a FASTA
reference file.  It is built mainly on top of the fastahack source code made
available by Erik Garrison <erik.garrison@bc.edu> of the Marth Lab at Boson College
at GitHub (https://github.com/ekg/fastahack).  Calculating the GC content in windows
of 1000 on hg18 took 54 seconds on a 3.06 GHz Intel Core 2 Duo with 8GB RAM.

Features:

- FASTA index (.fai) generation for FASTA files (~30 seconds for hg18 on 3GHz Core 2 Duo)
- GC content calculation for fixed non-overlapping windows
- N content calculation (anything not ACGT) for fixed non-overlapping windows

CLI:

	Usage: ./gcCounter [options] <FASTA reference>

	Options:
		-w, --window <int>       Specify the size of non-overlapping windows [1000]
		-l, --list               List all chromosomes in FASTA reference file
		-s, --sequence <string>  Specify the entries and order of sequences to analyze [ALL],
								 the <string> should be a comma-delimited list (NO spaces)
	Example:
		./gcCounter -w 100000 -s 1,3,5,X hg18.fasta > hg18.gc.seg

### readCounter

Fast and memory efficient counting of reads in BAM files

Overview:

mapCounter is a small program for counting the number of reads in non-overlapping
windows of fixed width directly from BAM files.  It is build on top of the [BamTools API]
made available by Derek Barnett on GitHub (https://github.com/pezmaster31/bamtools).

[BamTools API]: http://bioinformatics.oxfordjournals.org/content/27/12/1691

Features:

- Average mappability calculation for fixed non-overlapping windows

CLI:

	Usage: ./readCounter [options] <BAM file>

	Options:
	    -w, --window <int>       Specify the size of non-overlapping windows [1000]
	    -q, --quality <int>      Specify the mapping quality value below which reads are ignored

	    -l, --list               List all chromosomes in BAM reference file
	    -s, --sequence <string>  Specify the entries and order of sequences to analyze [ALL],
	                             the <string> should be a comma-delimited list (NO spaces)

	    -b, --build              Build BAM index for file (same index format as SAMtools)
	Example:
	    ./readCounter -w 100 -s 1,3,5,X aligned_reads.bam > readcounts.seg

## Generating the Genome Reference Mappability File

If you would like to generate your own mappability file, we provide some scripts to do so.

### TLDR:
The process is:
1) Obtain FASTA genome file of interest
2) Run the FASTA genome through `util/mappability/generateMap.pl`, which generates a 1-basepair resolution BigWig (.bw) file.  Note this is READLENGTH-specific. (have `--window` parameter match your expected read length)
3) Run mapCounter with the above .bw file to generate BIN SIZE specific .wig file (have `--window` parameter match your desired bin width)

### Example
Given `refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa`:

Note the key parameter here is `--window` which should roughly match the expected read length of the sequencing library (e.g. 125 for modern Illumina sequencers): 

	$HMMCOPY_DIR/util/mappability/generateMap.pl refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -o refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw

This may fail if you have not built a bowtie index of your reference yet. To do this, make sure first that bowtie is in your $PATH variable. Once this is done, you can go:

	$HMMCOPY_DIR/util/mappability/generateMap.pl -b refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -o refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw

This will build the bowtie index. Now you can re-run the first command. In the end, you should get a BigWig file which you then need to convert into a wig file with `mapCounter`, the key parameter here is `--window` which should match your desired bin width. You can do this by going:

	$HMMCOPY_DIR/bin/mapCounter -w 1000 refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw > refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.ws_1000.wig

#### On bin widths

Typically a good bin width is one that allows you at least 200 reads per bin.  For 30x this means roughly bins of 1000, and for single cell sequencing (at 0.05x) this is closer to 500,000.

Also, if your bin width EXCEEDS the size of any chromosome, you will end up with a chromosome with a single bin, which will cause the following error:

	Initialization
		    EM iteration: 1 Log likelihood: -Inf
		    Expectation
		    Error: INTEGER() can only be applied to a 'integer', not a 'NULL'

Some recommended solutions are:
1) Decrease the bin size when running readCounter, and generating the GC and mappability wig files
2) Remove the offending chromosome line from the three input files (reads, GC, mappability), simply open the file and delete the last two lines in each that correspond to the MT chromosome

### Generating the Genome Reference GC Content File

You can similarly create the GC wig file with gcCounter using the same reference.  Ensure your `--window` parameter here matches that of the other readCounter and mapCounter tools.

	$HMMCOPY_DIR/bin/gcCounter Homo_sapiens.GRCh37.75.dna.primary_assembly.fa > Homo_sapiens.GRCh37.75.dna.primary_assembly.gc.wig
