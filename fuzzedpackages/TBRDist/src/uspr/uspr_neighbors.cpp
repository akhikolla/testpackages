/*******************************************************************************
uspr_neighbors.cpp

Usage: uspr_neighbors [OPTIONS]
Calculate the set of trees that can be obtained from a given unrooted binary tree by a sin
The tree is read from standard input and must be in newick format.
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
bool IGNORE_ORIGINAL = false;
bool SIZE_ONLY = false;

string USAGE =
"uspr_neighbors, version 1.0.1\n"
"\n"
"usage: uspr_neighbors [OPTIONS]\n"
"Calculate all trees that can be obtained from an unrooted binary tree using a\n"
"single subtree prune and regraft operation.  Input one tree in newick format on\n"
"STDIN\n"
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
"--ignore_original      Do not output or count the input tree.\n"
"\n"
"--size_only            Count the number of SPR neighbors instead of printing them.\n";

// function prototypes


int main(int argc, char *argv[]) {
	int max_args = argc-1;
	while (argc > 1) {
		char *arg = argv[--argc];
		if (strcmp(arg, "--size_only") == 0) {
			SIZE_ONLY = true;
		}
		else if (strcmp(arg, "--ignore_original") == 0) {
			IGNORE_ORIGINAL = true;
		}
		else if (strcmp(arg, "--help") == 0 ||
				strcmp(arg, "-h") == 0 ||
				strcmp(arg, "-help") == 0) {
			Rcout << USAGE;
			return 0;
		}
	}

	// label maps to allow string labels
	map<string, int> label_map= map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();

	// set random seed
	srand(unsigned(time(0)));

	// read input trees
	string T_line = "";
	while (getline(cin, T_line)) {
		// load into data structures
		uforest F1 = uforest(T_line, &label_map, &reverse_label_map);
		F1.normalize_order();
		if (!IGNORE_ORIGINAL) {
			Rcout << F1.str(false, &reverse_label_map) << endl;
		}

		// compute neighbors
		list<utree> neighbors = get_neighbors(&F1);
		if (SIZE_ONLY) {
			int size = neighbors.size();
			Rcout << size << " neighbors" << endl;
		}
		else {
			for (utree T: neighbors) {
				Rcout << T.str(false, &reverse_label_map) << endl;
			}
		}
	}
}
