// structs.h: Imprecise Classification Trees
//
// Copyright (C) 2018  Paul Fink
//
// This file is part of imptree.
//
// imptree is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// imptree is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with imptree.  If not, see <https://www.gnu.org/licenses/>.

#ifndef RCPP_IMPTREE_STRUCTS_H
#define RCPP_IMPTREE_STRUCTS_H

#include <Rcpp.h>
#include "translation.h"
#include "enums.h"

struct ProbInterval {
  int obs;
  std::vector<int> freq;
  std::vector<double> lower;
  std::vector<double> upper;
  
  std::string to_string(const int nsmall = 6, const std::string &sep = "\t") const;
  Rcpp::NumericMatrix toMatrix() const;
};


struct Config {
  double s; // strictly positive
  double gamma; // in (0;1)
  double tbase; // in [-1,2]
  int minbucket; // non-negative int
  int maxdepth; // non-negative int
  EntropyCorrection ec;
  SplitMetric sm;
  IpType ip;
};

struct Data {
  Data() {};
  Data(const Rcpp::IntegerMatrix & mat);
  
  Rcpp::IntegerMatrix data;
  int classidx;
  Rcpp::IntegerVector nlevels;
  Rcpp::List labels;
  Rcpp::CharacterVector varnames;
  
};
#endif /*RCPP_IMPTREE_STRUCTS_H*/
