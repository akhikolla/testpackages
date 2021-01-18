// structs.cpp: Imprecise Classification Trees
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

#include <Rcpp.h>
#include "structs.h"

// Constructor for Data object
Data::Data(const Rcpp::IntegerMatrix & mat) : data(mat)
{
  classidx = data.attr("classidx");
  nlevels = data.attr("nlevels");
  labels = data.attr("labels");
  varnames = Rcpp::colnames(data);
}

std::string ProbInterval::to_string(const int nsmall, const std::string &sep) const {
  
  int ncat = freq.size();
  std::ostringstream out;
  out << std::fixed << std::setprecision(nsmall);
  for(int i = 0; i < ncat - 1; ++i) {
    out << "[" << lower[i] << ";" << upper[i] << "]" << sep;
  }
  out << "[" << lower[ncat-1] << ";" << upper[ncat-1] << "]";
  return out.str();
}

Rcpp::NumericMatrix ProbInterval::toMatrix() const {
  
  Rcpp::NumericMatrix m(freq.size(), 3);

  Rcpp::NumericVector freqr = Rcpp::wrap(freq);
  Rcpp::NumericVector lowerr = Rcpp::wrap(lower);
  Rcpp::NumericVector upperr = Rcpp::wrap(upper);
  m(Rcpp::_,0) = freqr;
  m(Rcpp::_,1) = lowerr;
  m(Rcpp::_,2) = upperr;
  
  Rcpp::NumericMatrix m2 = Rcpp::transpose(m);
  Rcpp::rownames(m2) = Rcpp::CharacterVector::create("Frequency", "Lower", "Upper");
  
  return m2;
}
