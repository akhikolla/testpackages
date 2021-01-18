// Copyright (C) 2014 Mohammad H. Ferdosi
//
// HSPhase is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// HSPhase program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


/*
 * Ohd.h
 *
 *  Created on: 02/08/2013
 *      Author: mhf
 */

#ifndef OHD_H_
#define OHD_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
RcppExport SEXP ohd(SEXP genotype, SEXP unique);

#endif /* OHD_H_ */
