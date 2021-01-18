/*
 * ncols.cpp
 *
 *  Created on: Jul 26, 2014
 *      Author: kaiyin
 */



#include "ncols.h"

//' Counts the number of columns of whitespace delimited file.
//'
//' @param fn Input filepath
//' @return A integer for the number of columns
//' @export
// [[Rcpp::export]]
long ncols(std::string fn) {
    fileExists(fn);
    std::ifstream in_file(fn.c_str());
    std::string tmpline;
    std::getline(in_file, tmpline);
    std::vector<std::string> strs;
    long numCols = 0;
    std::istringstream lineStream(tmpline);
    std::string tmpword = "";
    while(lineStream >> tmpword) {
    	if(not tmpword.empty()) {
    		numCols++;
    	}
    }
    return numCols;
}
