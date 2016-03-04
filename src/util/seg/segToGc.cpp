// ***************************************************************************
// segToGC.cpp (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
// Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 02 August 2011 (DL)
// ---------------------------------------------------------------------------


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <getopt.h>

#include "Fasta.h"

using namespace std;

void analyze(string sequence, double &gc_content, double &n_content, int &length) ;
void usage(const char* app);

int main (int argc, char** argv) {

// NOTE: 0-based positions internally

	string fasta_file;
	string seg_file;
	char sep = '\t';
    bool memory_map = false;

    if (argc == 3) {
        fasta_file = argv[1];
		seg_file = argv[2];
    } else {
        cerr << "Please specify a FASTA and SEG file." << endl;
		usage(argv[0]);
        exit(1);
    }

    FastaReference fr;
    fr.open(fasta_file, memory_map);

	ifstream seg_data(seg_file.c_str());
	if (seg_data.is_open()) {
		string line;
		cout << "data" << sep << "chr" << sep << "start" << sep << "end" << sep << "n" << sep << "gc" << endl;

		// Skip header if present
		getline(seg_data, line);
		vector<string> tokens = split(line, '\t');
		if (atoi(tokens[2].c_str())) { // Position should be int > 0
			seg_data.seekg(0, ios::beg);
		}

		while (seg_data.good()) {
			getline(seg_data, line);
			if (line.empty()) break;

			vector<string> tokens = split(line, '\t');
			string chr = tokens[1];
			int start = atoi(tokens[2].c_str());
			int size = atoi(tokens[3].c_str()) - start;

			string sequence = fr.getSubSequence(chr, start, size);
			double gc = 0;
			double n = 0;
			int length = 0;
			analyze(sequence, gc, n, length);

			// NOTE: SEG has 0-based positions
			cout << tokens[0] << sep << chr << sep << start << sep << start + length << sep << n << sep << gc << endl;
		}
		seg_data.close();
	} else {
		cerr << "Cannot open SEG file: " << seg_file << endl;
		exit(-1);
	}

	return 0;
}

void analyze(string sequence, double &gc_content, double &n_content, int &length) {
	int gc_count = 0;
	int n_count = 0;
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

	gc_content = (double)gc_count / (length - n_count);
	n_content = (double)n_count / length;

	return;
}

void usage(const char* app) {
	cerr << "Usage: " << app << " <FASTA reference> <SEG file>" << endl
		<< endl
		<< "    Calculates N and GC content for segments in the SEG file." << endl
		<< "    Ensure that the chromosome and positions in the SEG file match" << endl
		<< "    the entry and ranges of the reference in the FASTA file." << endl
		<< "    *NOTE*: Positions in SEG file should be 0-based positions" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " hg18.fa example.seg" << endl
		<< endl
		<< "Author: Daniel Lai <jujubix@cs.ubc.ca>" << endl
		<< endl
		<< "References:" << endl
		<< "    SEG: http://www.broadinstitute.org/software/igv/SEG" << endl;
}
