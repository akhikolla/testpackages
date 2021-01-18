/*
 * readcols.h
 *
 *  Created on: Aug 11, 2014
 *      Author: kaiyin
 */

#ifndef READCOLS_H_
#define READCOLS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include "ncols.h"
#include "countlines.h"


Rcpp::CharacterMatrix readcols(std::string fn,
		std::vector<long> colsel, long nFirstSkipLines,
		long nSkipUnit);

#endif /* READCOLS_H_ */
