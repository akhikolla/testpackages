/*
 * readcols.cpp
 *
 *  Created on: Aug 11, 2014
 *      Author: kaiyin
 */

#include "readcols.h"



//' Read columns of a whitespace delimited text file
//'
//' @param fn input filepath
//' @param colsel a vector of target column numbers
//' @param nFirstSkipLines Integer. Number of lines to skip in the beginning
//' @param nSkipUnit Integer M. Let the function read one line out of every M
//' @return A matrix of strings from selected columns
//' @export
// [[Rcpp::export]]
Rcpp::CharacterMatrix readcols(std::string fn,
		std::vector<long> colsel, long nFirstSkipLines,
		long nSkipUnit) {

	if (colsel.empty()) {
		throw std::string("You didn't select any column!");
	}


	long nc_file = ncols(fn);
	long nr_file = countlines(fn);
	// skipping too many lines is an error
	if(nFirstSkipLines >= nr_file) {
		throw std::string("More lines skip than total number of lines.");
	}
	long nr_left = nr_file - nFirstSkipLines;
	long nr = (long) (nr_left / nSkipUnit);
	long nc = colsel.size();



	{
		long remainder = (long) ((nr_file - nFirstSkipLines) % nSkipUnit);
		if (remainder != 0) {
			nr++;
		}
		else {
		}
	}

	// c++ is 0-based, adjust for it
	for (unsigned long i = 0; i < colsel.size(); i++) {
		// indices must be in the right range. i.e. positve integers
		if(colsel[i] < 1) {
			throw std::string("Column index must be positive integer.");
		}
		colsel[i]--;
	}

	long colsel_max = *std::max_element(colsel.begin(), colsel.end());
	if (((long)colsel_max) >= nc_file) {
		throw std::string("Some col number(s) are out of range!");
	}


	Rcpp::CharacterMatrix res(nr, nc);
	std::ifstream infile(fn.c_str());
	std::string tmpline;

	// skip lines in the beginning
	for (long lineIter = 0; lineIter < nFirstSkipLines; lineIter++) {
		getline(infile, tmpline);
	}


	long rowIter = 0;
	for (long lineIter = 0; lineIter < nr_left; lineIter++) {
		std::string line;
		getline(infile, line);
		if (lineIter % nSkipUnit == 0) {
			std::istringstream lineStream(line);
			long colIter = 0;
			for (long wordIter = 0; wordIter <= (long)colsel_max; wordIter++) {
				std::string tmpword;
				lineStream >> tmpword;
				if (std::find(colsel.begin(), colsel.end(), wordIter)
						!= colsel.end()) {
					res(rowIter, colIter) = tmpword;
					colIter++;
				}
			}
			rowIter++;
		}
	}

    return res;

}

