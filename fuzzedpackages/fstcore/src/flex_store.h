/*
 fstcore - R bindings to the fstlib library

 Copyright (C) 2017-present, Mark AJ Klik

 This file is part of the fstcore R package.

 The fstcore R package is free software: you can redistribute it and/or modify it
 under the terms of the GNU Affero General Public License version 3 as
 published by the Free Software Foundation.

 The fstcore R package is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
 for more details.

 You should have received a copy of the GNU Affero General Public License along
 with the fstcore R package. If not, see <http://www.gnu.org/licenses/>.

 You can contact the author at:
 - fstcore R package source repository : https://github.com/fstpackage/fstcore
*/


#ifndef FASTSTORE_H
#define FASTSTORE_H

// [[Rcpp::interfaces(cpp)]]

#include <Rcpp.h>


// [[Rcpp::export]]
int fstlib_version();


// [[Rcpp::export]]
SEXP fststore(Rcpp::String fileName, SEXP table, SEXP compression, SEXP uniformEncoding);


// [[Rcpp::export]]
SEXP fstmetadata(Rcpp::String fileName);


// [[Rcpp::export]]
SEXP fstretrieve(Rcpp::String fileName, SEXP columnSelection, SEXP startRow, SEXP endRow);


#endif  // FASTSTORE_H
