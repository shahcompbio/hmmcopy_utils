correctReadcount suite - various tools for correcting GC content bias in readcounts

Author: Daniel Lai <jujubix@cs.ubc.ca> Department of Molecular Oncology, BC Cancer Research Agency
Date:   August 02, 2011

Installation:

	Requires cmake installed and on PATH (http://www.cmake.org/).
	In the source directory, install as follows:
		cmake .
		make
	Should result in main binaries in /bin and useful binaries in /util

/********** mapCounter **********/

mapCounter --- *fast* average mappability counter using BigWig files

Overview:

	mapCounter is a small program for calculating average mappability[1] for non-overlapping
	windows of fixed width across all sequences (i.e. chromosomes) present in a BigWig file[2].
	It is built mainly on top of an independent subset of files obtained from the UCSC Genome
	Browser source code (i.e. kent library)[3] made available by Jim Kent.  Generating average
	mappability files in 1000 base windows (default) on hg18 took 160 seconds on a 3.06 GHz Intel
	Core 2 Duo with 8GB RAM.

	[1] http://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeMapability
	[2] http://genome.ucsc.edu/goldenPath/help/bigWig.html
	[3] http://genome.ucsc.edu/admin/git.html

Features:

	- Average mappability calculation for fixed non-overlapping windows

Usage: ./mapCounter [options] <BigWig file>

Options:
    -w, --window <int>       Specify the size of non-overlapping windows [1000]
    -l, --list               List all chromosomes in BigWig file
    -s, --sequence <string>  Specify the entries and order of sequences to analyze [ALL],
                             the <string> should be a comma-delimited list (NO spaces)
Example:
    ./mapCounter -w 100000 -s 1,3,5,X hg18.bw > hg18.map.seg

/********** gcCounter **********/

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

Usage: ./gcCounter [options] <FASTA reference>

Options:
    -w, --window <int>       Specify the size of non-overlapping windows [1000]
    -l, --list               List all chromosomes in FASTA reference file
    -s, --sequence <string>  Specify the entries and order of sequences to analyze [ALL],
                             the <string> should be a comma-delimited list (NO spaces)
Example:
    ./gcCounter -w 100000 -s 1,3,5,X hg18.fasta > hg18.gc.seg

/********** readCounter **********/

readCounter --- fast and memory efficient counting of reads in BAM files

Overview:

	mapCounter is a small program for counting the number of reads in non-overlapping
	windows of fixed width directly from BAM files.  It is build on top of the BamTools
	API[1] made available by Derek Barnett on GitHub (https://github.com/pezmaster31/bamtools).

	[1] http://bioinformatics.oxfordjournals.org/content/27/12/1691

Features:

	- Average mappability calculation for fixed non-overlapping windows

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