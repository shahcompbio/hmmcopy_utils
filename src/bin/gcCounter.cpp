// ***************************************************************************
// gcCounter.cpp (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
// Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 05 August 2011 (DL)
// ---------------------------------------------------------------------------

#include "Fasta.h"
#include <stdlib.h>
#include <getopt.h>

void analyze(string sequence, double &gc_content, double &n_content, long &length) ;
void usage(const char* app);

bool forgiving = 0;

int main (int argc, char** argv) {

// NOTE: 0-based positions internally

	string fasta_file;
	long window_size = 1000;
	string sequence_list = "";
    bool memory_map = false;
    int list = 0;
	char sep = '\t';
	bool seg = 0;
	bool ncontent = 0;

	int c;
	while (true) {
		static struct option long_options[] = {
			{"help",		no_argument, 0, 'h'},
			{"list",		no_argument, 0, 'l'},
			{"seg",			no_argument, 0, 's'},
			{"n",			no_argument, 0, 'n'}, // Hidden option, outputs N content instead of GC
			{"forgiving",	no_argument, 0, 'f'}, // Hidden option, allows GC content for bins with N
			{"window",		required_argument, 0, 'w'},
			{"chromosome",	required_argument, 0, 'c'},
			{0, 0, 0, 0}
		};

		int option_index = 0;
		c = getopt_long (argc, argv, "hlnsfw:c:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c) {
			case 's':
				seg = 1;
				break;

			case 'n':
				ncontent = 1;
				break;

			case 'f':
				forgiving = 1;
				break;

			case 'w':
				window_size = atoi(optarg);
				break;

			case 'c':
				sequence_list = optarg;
				break;

			case 'h':
				usage(argv[0]);
				exit(0);
				break;

			case 'l':
				list = 1;
				break;

			default:
				usage(argv[0]);
				exit(1);
		}
	}

    if (optind < argc) {
        fasta_file = argv[optind];
    } else {
        cerr << "Please specify a FASTA reference file." << endl;
		usage(argv[0]);
        exit(1);
    }

    FastaReference fr;
    fr.open(fasta_file, memory_map);

	if (list) {
		cerr << "Listing chromosome names and lengths in file: " << fasta_file << endl;
		for (size_t i = 0; i < fr.index->sequenceNames.size(); ++i) {
			string name = fr.index->sequenceNames[i];
			FastaIndexEntry seq = fr.index->entry(name);
			cout << name << "\t" << seq.length << endl;
		}
		exit(0);
	}

	vector<string> chromosomes = fr.index->sequenceNames;
	if (!sequence_list.empty()) {
		bool quit = 0;
		chromosomes = split(sequence_list, ",");

		for (size_t i = 0; i < chromosomes.size(); ++i) {
			string name = chromosomes[i];
			FastaIndexEntry seq = fr.index->entry(name);
			if (seq.length == 0) {
				quit = 1;
			}
		}
		if (quit) {
			cerr << "Use --list for list of valid chromosome names." << endl;
			exit(-1);
		}
	}

	if (seg) cout << "data" << sep << "chr" << sep << "start" << sep << "end" << sep << "n" << sep << "gc" << endl;

	for (size_t i = 0; i < chromosomes.size(); ++i) {
		string name = chromosomes[i];
		FastaIndexEntry seq = fr.index->entry(name);
		long pos = 0;
		long window = window_size;

		if (window_size > seq.length) {
			cerr << "Window size greater than chromosome length of " << name << ", adjusting to chromosome length: " << seq.length << endl;
			window = seq.length;
		}

		// NOTE: WIG has 1-based positions
		if (!seg) cout << "fixedStep chrom=" << name << " start=1 step=" << window << " span=" << window << endl;

		while(pos < seq.length) {
			string sequence = fr.getSubSequence(name, pos, window);
			// cout << name << "\t" << pos << "\t" <<  sequence << endl;

			double gc = 0;
			double n = 0;
			long length = 0;
			analyze(sequence, gc, n, length);


			if (seg) {
				// NOTE: SEG uses 0-based positions
				cout << "gc" << sep << name << sep << pos << sep << pos + length << sep << n << sep << gc << endl;
			} else {
				if (ncontent) {
					cout << n << endl;
				} else {
					cout << gc << endl;
				}
			}
			pos += window;
		}
	}

	return 0;
}

void analyze(string sequence, double &gc_content, double &n_content, long &length) {
	long gc_count = 0;
	long n_count = 0;
	length = 0;

	for (size_t i = 0; i < sequence.length(); ++i) {
		switch(sequence.at(i)) {
			case 'G':
			case 'C':
			case 'g':
			case 'c':
				++gc_count;
				++length;
				break;
			case 'A':
			case 'T':
			case 'a':
			case 't':
				++length;
				break;
			case '>':
				goto END;
			default:
				++n_count;
				++length;
				break;
		}
	}

	END:

	bool invalid = (n_count > 0) ? 1 : 0;
	if (forgiving && n_count < length) {
		invalid = 0;
	}

	gc_content = (invalid) ? -1 : (double)gc_count / (length - n_count);
	n_content = (double)n_count / length;

	return;
}

void usage(const char* app) {
	cerr << "Usage: " << app << " [options] <FASTA reference>" << endl
		<< endl
		<< "Options:" << endl
		<< "    -s, --seg                 Outputs in SEG format" << endl
		<< "    -w, --window <int>        Specify the size of non-overlapping windows [1000]" << endl
		<< "    -l, --list                List all chromosomes in FASTA reference file" << endl
		<< "    -c, --chromosome <string> Specify the entries and order of sequences to analyze [ALL]," << endl
		<< "                              the <string> should be a comma-delimited list (NO spaces)" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " -w 100000 -c 1,3,5,X hg18.fasta > hg18.gc.wig" << endl
		<< endl
		<< "Author: Daniel Lai <jujubix@cs.ubc.ca>" << endl;
}
