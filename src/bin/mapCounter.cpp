// ***************************************************************************
// mapCounter.cpp (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>
// Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 05 August 2011 (DL)
// ---------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
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

void usage(const char* app);

int main (int argc, char** argv) {

	// NOTE: 0-based positions internally

	string bigwig_file;
	int window_size = 1000;
	string sequence_list = "";
    int list = 0;
	char sep = '\t';
	bool seg = 0;

	int c;
	while (true) {
		static struct option long_options[] = {
			{"seg",			no_argument, 0, 's'},
			{"help",		no_argument, 0, 'h'},
			{"list",		no_argument, 0, 'l'},
			{"window",		required_argument, 0, 'w'},
			{"chromosome",	required_argument, 0, 'c'},
			{0, 0, 0, 0}
		};

		int option_index = 0;
		c = getopt_long (argc, argv, "hlsw:c:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c) {
			case 's':
				seg = 1;
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
        bigwig_file = argv[optind];
    } else {
        cerr << "Please specify a BigWig file." << endl;
		usage(argv[0]);
        exit(1);
    }

    struct bbiFile *bwf = bigWigFileOpen((char *)bigwig_file.c_str());
	struct bbiChromInfo *info = bbiChromList(bwf);
	vector<string> chromosomes;
	while (info != NULL) {
		chromosomes.push_back(info->name);
		info = info->next;
	}
	bbiChromInfoFreeList(&info);

	if (list) {
		cerr << "Listing chromosome names and lengths in file: " << bigwig_file << endl;
		for (int i = 0; i < chromosomes.size(); ++i) {
			cout << chromosomes[i] << sep << bbiChromSize(bwf, (char *)chromosomes[i].c_str()) << endl;
		}
		exit(0);
	}

	if (!sequence_list.empty()) {
		bool quit = 0;
		chromosomes = split(sequence_list, ",");

		for (size_t i = 0; i < chromosomes.size(); ++i) {
			if (bbiChromSize(bwf, (char *)chromosomes[i].c_str()) == 0) {
				cerr << "Unrecognized chromosome name: " << chromosomes[i] << endl;
				quit = 1;
			}
		}
		if (quit) {
			cerr << "Use --list for list of valid chromosome names." << endl;
			exit(-1);
		}
	}

	if (seg) cout << "data" << sep << "chr" << sep << "start" << sep << "end" << sep << "mappability" << endl;

	for (size_t i = 0; i < chromosomes.size(); ++i) {
		char * name = (char *)chromosomes[i].c_str();
		struct bigWigValsOnChrom *chromVals = bigWigValsOnChromNew();

		if (bigWigValsOnChromFetchData(chromVals, name, bwf)) {
			long pos = 0;
			double * ptr = chromVals->valBuf;
			long count = 0;
			double sum = 0;
			long window = window_size;

			if (window_size > chromVals->chromSize) {
				cerr << "Window size greater than chromosome length of " << name << ", adjusting to chromosome length: " << chromVals->chromSize << endl;
				window = chromVals->chromSize;
			}

			// NOTE: WIG has 1-based positiosn
			if (!seg) cout << "fixedStep chrom=" << name << " start=1 step=" << window << " span=" << window << endl;

			while (pos < chromVals->chromSize) {
				sum += (*ptr);
				++ptr;
				++count;
				++pos;
				if (count == window || pos == chromVals->chromSize) {
					if (seg) {
						// NOTE: SEG has 0-based positions
						cout << "map" << sep << name << sep << pos - count << sep << pos << sep << sum/count << endl;
					} else {
						cout << sum/count << endl;
					}
					count = 0;
					sum = 0;
				}
			}
		} else {
			cerr << "No data on chromosome: " << name << endl;
			exit(-1);
		}

		bigWigValsOnChromFree(&chromVals);
	}

	bigWigFileClose(&bwf);

	return 0;
}

void usage(const char* app) {
	cerr << "Usage: " << app << " [options] <BigWig file>" << endl
		<< endl
		<< "Options:" << endl
		<< "    -s, --seg                 Outputs in SEG format" << endl
		<< "    -w, --window <int>        Specify the size of non-overlapping windows [1000]" << endl
		<< "    -l, --list                List all chromosomes in BigWig file" << endl
		<< "    -c, --chromosome <string> Specify the entries and order of sequences to analyze [ALL]," << endl
		<< "                              the <string> should be a comma-delimited list (NO spaces)" << endl
		<< "    -h, --help                This help message" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " -w 100000 -c 1,3,5,X hg18.bw > hg18.map.wig" << endl
		<< endl
		<< "Author: Daniel Lai <jujubix@cs.ubc.ca>" << endl;
}

