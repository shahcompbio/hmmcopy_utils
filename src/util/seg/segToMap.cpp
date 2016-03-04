// ***************************************************************************
// segToMap.cpp (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
// Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 29 July 2011 (DL)
// ---------------------------------------------------------------------------

// NOTE: Some loss in precision is observed...

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <getopt.h>

#include "split.h"

extern "C" {
	#include "common.h"
	#include "bbiFile.h"
	#include "bigWig.h"
}

using namespace std;

void usage();

int main (int argc, char** argv) {

// NOTE: 0-based positions internally

	string bigwig_file;
	string seg_file;
	char sep = '\t';

    if (argc == 3) {
        bigwig_file = argv[1];
		seg_file = argv[2];
    } else {
        cerr << "Please specify a BigWig and SEG file." << endl;
		usage();
        exit(1);
    }

    struct bbiFile *bwf = bigWigFileOpen((char *)bigwig_file.c_str());

	ifstream seg_data(seg_file.c_str());
	if (seg_data.is_open()) {
		string line;
		cout << "data" << sep << "chr" << sep << "start" << sep << "end" << sep << "mappability" << endl;

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
			double value = bigWigSingleSummary(bwf, (char *)tokens[1].c_str(),
				atoi(tokens[2].c_str()), atoi(tokens[3].c_str()), bbiSummaryTypeFromString("mean"), 0);
			cout << tokens[0] << sep << tokens[1] << sep << tokens[2] << sep << tokens[3] << sep << value << endl;

		}
		seg_data.close();
	} else {
		cerr << "Cannot open SEG file: " << seg_file << endl;
		exit(-1);
	}

	bigWigFileClose(&bwf);

	return 0;
}

void usage() {
	cerr << "Usage: segToMap <BigWig file> <SEG file>" << endl
		<< endl
		<< "    Calculates average mappability for segments in the SEG file." << endl
		<< "    Ensure that the chromosome and positions in the SEG file match" << endl
		<< "    the name and ranges of the data in the BigWig mappability file" << endl
		<< "    *NOTE*: Positions in SEG file must be 0-based positions" << endl
		<< endl
		<< "Example:" << endl
		<< "    ./segToMap hg18.bw example.seg" << endl
		<< endl
		<< "Author: Daniel Lai <jujubix@cs.ubc.ca>" << endl
		<< endl
		<< "References:" << endl
		<< "    BigWig: http://genome.ucsc.edu/goldenPath/help/bigWig.html" << endl
		<< "    SEG: http://www.broadinstitute.org/software/igv/SEG" << endl;
}

