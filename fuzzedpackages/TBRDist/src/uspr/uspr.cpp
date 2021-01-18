/*******************************************************************************
uspr.cpp

Usage: uspr [OPTIONS]
Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of unrooted binary trees from STDIN in newick format.
Supports arbitrary leaf labels. See the README for more information.

Copyright 2018 Chris Whidden
cwhidden@fredhutch.org
https://github.com/cwhidden/uspr
May 1, 2018
Version 1.0.1

This file is part of uspr.

uspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

uspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uspr.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/

// includes
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <memory>
#include <ctime>
#include <cstdlib>
#include "utree.h"
#include "unode.h"
#include "uforest.h"
#include "tbr.h"
#include "uspr.h"

using namespace std;

// constants
//

bool DEFAULT_OPTIMIZATIONS = true;
bool PRINT_mAFS = false;
bool COUNT_mAFS = false;

bool ALL_DISTANCES = true;
bool COMPUTE_TBR_APPROX = false;
bool COMPUTE_TBR = false;
bool COMPUTE_REPLUG = false;
bool COMPUTE_USPR = false;
string USAGE =
"uspr, version 1.0.1\n"
"\n"
"usage: uspr [OPTIONS]\n"
"Calculate the subtree prune and regraft (uSPR) distance and related distances\n"
"between pairs of unrooted binary tres from STDIN. Input one tree in newick\n"
"format per line. Supports arbitrary labels. See the README for more\n"
"information.\n"
"\n"
"Copyright 2018 Chris Whidden\n"
"cwhidden@fredhutch.org\n"
"https://github.com/cwhidden/uspr\n"
"May 1, 2018\n"
"Version 1.0.1\n"
"\n"
"This program comes with ABSOLUTELY NO WARRANTY.\n"
"This is free software, and you are welcome to redistribute it\n"
"under certain conditions; See the README for details.\n"
"\n"
"Basic options\n"
"\n"
"-h --help              Print program information and exit.\n"
"\n"
"--tbr-approx\n"
"--tbr\n"
"--replug\n"
"--uspr                 By default, uspr will compute all 4 distances. If any of these\n"
"                       options are specified then uspr will compute only the specified\n"
"                       distances.\n"
"                       \n"
"--print-mAFs           Print all maximal agreement forests found during execution of\n"
"                       the program.\n"
"--count-mAFs           Count all maximal agreement forests found during execution of\n"
"                       the program.\n"
"\n"
"Algorithm Options\n"
"\n"
"--no-approx-estimate\n"
"--no-tbr-estimate\n"
"--no-replug-estimate   Disable the TBR approximation, TBR distance, or replug distance\n"
"                       heuristics when computing the SPR distance. In most cases these\n"
"                       options will greatly increase the time required by uspr.\n"
"                       \n"
"--no-opt\n"
"--no-protect-b         Disable all MAF optimizations or just the edge protection\n"
"                       optimization for enumerating agreement forests. In most\n"
"                       cases these options will greatly increase the time required\n"
"                       by uspr.\n"
"\n";


// function prototypes


int main(int argc, char *argv[]) {
	int max_args = argc-1;
	while (argc > 1) {
		char *arg = argv[--argc];
		if (strcmp(arg, "--print-mAFs") == 0) {
			PRINT_mAFS = true;
		}
		if (strcmp(arg, "--count-mAFs") == 0) {
			COUNT_mAFS = true;
		}
		else if (strcmp(arg, "--keep-labels") == 0) {
			KEEP_LABELS = true;
				Rcout << "KEEP_LABELS=true" << endl;
		}
		else if (strcmp(arg, "--no-opt") == 0) {
			DEFAULT_OPTIMIZATIONS = false;
		}
		else if (strcmp(arg, "--no-protect-b") == 0) {
			OPTIMIZE_PROTECT_B = false;
		}
		else if (strcmp(arg, "--tbr-approx") == 0) {
			COMPUTE_TBR_APPROX = true;
			ALL_DISTANCES = false;
		}
		else if (strcmp(arg, "--tbr") == 0) {
			COMPUTE_TBR = true;
			ALL_DISTANCES = false;
		}
		else if (strcmp(arg, "--replug") == 0) {
			COMPUTE_REPLUG = true;
			ALL_DISTANCES = false;
		}
		else if (strcmp(arg, "--uspr") == 0) {
			COMPUTE_USPR = true;
			ALL_DISTANCES = false;
		}
		else if (strcmp(arg, "--no-approx-estimate") == 0) {
			USE_TBR_APPROX_ESTIMATE = false;
		}
		else if (strcmp(arg, "--no-tbr-estimate") == 0) {
			USE_TBR_ESTIMATE = false;
		}
		else if (strcmp(arg, "--no-replug-estimate") == 0) {
			USE_REPLUG_ESTIMATE = false;
		}
		else if (strcmp(arg, "--help") == 0 ||
				strcmp(arg, "-h") == 0 ||
				strcmp(arg, "-help") == 0) {
			Rcout << USAGE;
			return 0;
		}
	}

	if (DEFAULT_OPTIMIZATIONS == false) {
		OPTIMIZE_2B = false;
		OPTIMIZE_PROTECT_A = false;
		OPTIMIZE_PROTECT_B = false;
		OPTIMIZE_BRANCH_AND_BOUND = false;
		Rcout << "NO OPTIMIZATIONS" << endl;
	}

	if (ALL_DISTANCES) {
		COMPUTE_TBR_APPROX = true;
		COMPUTE_TBR = true;
		COMPUTE_REPLUG = true;
		COMPUTE_USPR = true;
	}

	// label maps to allow string labels
	map<string, int> label_map= map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();

	// set random seed
	srand(unsigned(time(0)));

	// read input trees
	string T1_line = "";
	string T2_line = "";
	while (getline(cin, T1_line) && getline(cin, T2_line)) {
		// load into data structures
		uforest F1 = uforest(T1_line, &label_map, &reverse_label_map);
		F1.normalize_order();
		uforest F2 = uforest(T2_line, &label_map, &reverse_label_map);
		F2.normalize_order();
		Rcout << "T1: " << F1.str(false, &reverse_label_map) << endl;
		Rcout << "T2: " << F2.str(false, &reverse_label_map) << endl;
		// compute TBR distance
		if (COMPUTE_TBR_APPROX) {
		Rcout << "a_TBR: " << tbr_high_lower_bound(F1, F2) << " <= d_TBR <= " << tbr_low_upper_bound(F1, F2) << endl;
		}
		if (COMPUTE_TBR) {
			uforest *MAF1 = NULL;
			uforest *MAF2 = NULL;
			int distance = tbr_distance(F1, F2, false, &MAF1, &MAF2);
			Rcout << "d_TBR = " << distance << endl;
			if (MAF1 != NULL) {
				Rcout << "F1: " << MAF1->str(false, &reverse_label_map) << endl;
				delete MAF1;
			}
			if (MAF2 != NULL) {
				Rcout << "F2: " << MAF2->str(false, &reverse_label_map) << endl;
				delete MAF2;
			}
		}
		int count;
		if (PRINT_mAFS) {
			count = tbr_print_mAFs(F1, F2);
			Rcout << count << " mAFs" << endl;
		}
		else if (COUNT_mAFS) {
			count = tbr_count_mAFs(F1, F2);
			Rcout << count << " mAFs" << endl;
		}

		if (COMPUTE_REPLUG) {
			uforest *MAF1 = NULL;
			uforest *MAF2 = NULL;
			int d_replug = replug_distance(F1, F2, false, &MAF1, &MAF2);
			Rcout << "d_R = " << d_replug << endl;
			if (MAF1 != NULL) {
				Rcout << "F1: " << MAF1->str(false, &reverse_label_map) << endl;
				delete MAF1;
			}
			if (MAF2 != NULL) {
				Rcout << "F2: " << MAF2->str(false, &reverse_label_map) << endl;
				delete MAF2;
			}
		}

		if (COMPUTE_USPR) {
			int d_uspr = uspr_distance(F1, F2);
			Rcout << "d_USPR = " << d_uspr << endl;
		}

	}
}
