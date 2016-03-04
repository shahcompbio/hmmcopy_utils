// ***************************************************************************
// readCounter.cpp 
// (c) 2011 Daniel Lai <jujubix@cs.ubc.ca>, Gavin Ha <gha@bccrc.ca>
// Shah Lab, Department of Molecular Oncology, BC Cancer Research Agency
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 14 December, 2012 
// ---------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <getopt.h>

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "split.h"

using namespace BamTools;
using namespace std;

void usage(const char* app);

int main(int argc, char *argv[]) {

	string bamfile;
	long window_size = 1000;
	string sequence_list = "";
    int list = 0;
	int quality = 0;
	char sep = '\t';
	int build = 0;
	string sample = "reads";
	bool seg = 0;

	int c;
	while (true) {
		static struct option long_options[] = {
			{"seg",			no_argument, 0, 's'},
			{"help",		no_argument, 0, 'h'},
			{"list",		no_argument, 0, 'l'},
			{"build",		no_argument, 0, 'b'},
			{"quality",		required_argument, 0, 'q'},
			{"window",		required_argument, 0, 'w'},
			{"chromosome",	required_argument, 0, 'c'},
			{0, 0, 0, 0}
		};

		int option_index = 0;
		c = getopt_long (argc, argv, "shlbw:c:q:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c) {
			case 's':
				seg = 1;
				break;

			case 'b':
				build = 1;
				break;

			case 'q':
				quality = atoi(optarg);
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
        bamfile = argv[optind];
    } else {
        cerr << "Please specify a BAM file." << endl;
		usage(argv[0]);
        return -1;
    }

	BamReader reader;
	if (!reader.Open(bamfile) ) {
	    cerr << "Could not open input BAM files." << endl;
	    return -1;
	}

	if (build) {
		cerr << "Building index for " << bamfile << ", please hold" << endl;

		if (reader.CreateIndex()) {
			cerr << "Index creation successful" << endl;
			return 0;
		} else {
			cerr << "Index creation failed" << endl;
			return -1;
		}
	}

	if(!reader.LocateIndex()) {
		cerr << "Could not locate valid BAM index for input file. Build index with -b" << endl;
		return -1;
	}

	vector<RefData> all = reader.GetReferenceData();

	if (list) {
		cerr << "Listing chromosome names and lengths in file: " << bamfile << endl;
		for (size_t i = 0; i < all.size(); ++i) {
			cout << all[i].RefName << sep << all[i].RefLength << endl;
		}
		return 0;
	}

	vector<RefData> chromosomes;
	if (!sequence_list.empty()) {
		bool quit = 0;
		vector<string> targets = split(sequence_list, ",");
		for (size_t i = 0; i < targets.size(); ++i) {
			int refID = reader.GetReferenceID(targets[i]);
			if (refID < 0) {
				cerr << "Unrecognized chromosome name: " << targets[i] << endl;
				quit = 1;
			} else {
				assert(targets[i].compare(all[refID].RefName) == 0);
				chromosomes.push_back(all[refID]);
			}
		}
		if (quit) {
			cerr << "Use --list for list of valid chromosome names." << endl;
			exit(-1);
		}
	} else {
		chromosomes = all;
	}

	if (seg) cout << "data" << sep << "chr" << sep << "start" << sep << "end" << sep << "count" << endl;

	for(size_t i = 0; i < chromosomes.size(); ++i) {
		string refName = chromosomes[i].RefName;
		long refLength = chromosomes[i].RefLength;
		int refID = reader.GetReferenceID(refName);
		long window = window_size;
		// cout << refID << ": " << refName << " (" << refLength << ")" << endl;

		if (window_size > refLength) {
			cerr << "Window size greater than chromosome length of " << refName << ", adjusting to chromosome length: " << refLength << endl;
			window = refLength;
		}

		// NOTE: WIG has 1-based positiosn
		if (!seg) cout << "fixedStep chrom=" << refName << " start=1 step=" << window << " span=" << window << endl;

		long start = 0;
		long end = start + window;
		if (reader.SetRegion(refID, 0, refID, refLength)) {
				long count = 0;
				BamAlignment read;
				while(reader.GetNextAlignmentCore(read)) {
					while(read.Position > end) {

						if (seg) {
							// NOTE: SEG has 0-based positions
							cout << sample << sep << refName << sep << start << sep << end << sep << count << endl;
						} else {
							cout << count << endl;
						}

						count = 0;
						start += window;
						end = (start + window > refLength) ? refLength : start + window;
					}
					if (read.MapQuality >= quality && !read.IsDuplicate()) { 
						++count;
						// cout << count << " " << read.Position << " (" << read.MapQuality << ")" << endl;
					}
				}

				while(end < refLength) {
					if (seg) {
						// NOTE: SEG has 0-based positions
						cout << sample << sep << refName << sep << start << sep << end << sep << count << endl;
					} else {
						cout << count << endl;
					}
					count = 0;
					start += window;
					end = (start + window > refLength) ? refLength : start + window;
				}

				if (seg) {
					// NOTE: SEG has 0-based positions
					cout << sample << sep << refName << sep << start << sep << end << sep << count << endl;
				} else {
					cout << count << endl;
				}
		} else {
			cerr << "Cannot retrieve chromosome " << refName << " from BAM file, file may be corrupt." << endl;
			exit(-1);
		}
	}

	reader.Close();
	return 0;
}

void usage(const char* app) {
	cerr << "Usage: " << app << " [options] <BAM file>" << endl
		<< endl
		<< "Options:" << endl
		<< "    -s, --seg                 Outputs in SEG format" << endl
		<< "    -w, --window <int>        Specify the size of non-overlapping windows [1000]" << endl
		<< "    -q, --quality <int>       Specify the mapping quality value below which reads are ignored" << endl
		<< endl
		<< "    -l, --list                List all chromosomes in BAM reference file" << endl
		<< "    -c, --chromosome <string> Specify the entries and order of sequences to analyze [ALL]," << endl
		<< "                              the <string> should be a comma-delimited list (NO spaces)" << endl
		<< endl
		<< "    -b, --build               Build BAM index for file (same index format as SAMtools)" << endl
		<< endl
		<< "Example:" << endl
		<< "    " << app << " -w 100 -c 1,3,5,X aligned_reads.bam > readcounts.wig" << endl
		<< endl
		<< "Author: Daniel Lai <jujubix@cs.ubc.ca>" << endl;
}
