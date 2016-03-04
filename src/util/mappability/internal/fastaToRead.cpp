// ***************************************************************************
// gccounter.cpp (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
// Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 26 July 2011 (DL)
// ---------------------------------------------------------------------------

#include "Fasta.h"
#include <stdlib.h>
#include <getopt.h>

void analyze(string sequence, double &gc_content, double &n_content, int &length) ;
void help();

int main (int argc, char** argv) {

// NOTE: 0-based positions internally

	string fasta_file;
	int window_size = 35;
	string sequence_list = "";
    bool memory_map = false;
    int list = 0;
	char sep = '\t';

	int c;
	while (true) {
		static struct option long_options[] = {
			{"help",		no_argument, 0, 'h'},
			{"list",		no_argument, 0, 'l'},
			{"window",		required_argument, 0, 'w'},
			{"chromosome",	required_argument, 0, 'c'},
			{0, 0, 0, 0}
		};

		int option_index = 0;
		c = getopt_long (argc, argv, "hlw:c:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c) {
			case 'w':
				window_size = atoi(optarg);
				break;

			case 'c':
				sequence_list = optarg;
				break;

			case 'h':
				help();
				exit(0);
				break;

			case 'l':
				list = 1;
				break;

			default:
				help();
				exit(1);
		}
	}

    if (optind < argc) {
        fasta_file = argv[optind];
    } else {
        cerr << "Please specify a FASTA reference file." << endl;
		help();
        exit(1);
    }

    FastaReference fr;
    fr.open(fasta_file, memory_map);

	if (list) {
		cerr << "Listing chromosome names and lengths in file: " << fasta_file << endl;
		for (size_t i = 0; i < fr.index->sequenceNames.size(); ++i) {
			string name = fr.index->sequenceNames[i];
			FastaIndexEntry seq = fr.index->entry(name);
			cout << name << sep << seq.length << endl;
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

	for (size_t i = 0; i < chromosomes.size(); ++i) {
		string name = chromosomes[i];
		FastaIndexEntry seq = fr.index->entry(name);
		long end = seq.length - window_size + 1;
		if (end <= 0) {
			cerr << "Window size " << window_size << " greater than sequence length of " << name << endl;
			return -1;
		}

		for (int pos = 0; pos < seq.length - window_size + 1; ++pos) {
			string sequence = fr.getSubSequence(name, pos, window_size);
			cout << ">" << name << ":" << pos << endl;
			cout << sequence << endl;
		}
	}

	return 0;
}

void help() {
	cerr << "Usage: fastaToRead [options] <FASTA reference>" << endl
		<< endl
		<< "Options:" << endl
		<< "    -w, --window <int>       Specify the size of the overlapping reads [1000]" << endl
		<< "    -l, --list               List all chromosomes in FASTA reference file" << endl
		<< "    -s, --sequence <string>  Specify the entries and order of sequences to analyze [ALL]," << endl
		<< "                             the <string> should be a comma-delimited list (NO spaces)" << endl
		<< "Example:" << endl
		<< "    ./fastaToRead -w 10 -s 1,3,5,X hg18.fasta > bowtie hg18 -f -" << endl
		<< endl
		<< "Author: Daniel Lai <jujubix@cs.ubc.ca>" << endl;
}
