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


#ifndef OPEN_MP_H
#define OPEN_MP_H

// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>


// [[Rcpp::export]]
SEXP getnrofthreads();


// [[Rcpp::export]]
int setnrofthreads(SEXP nrOfThreads);


// [[Rcpp::export]]
SEXP hasopenmp();


// [[Rcpp::export]]
void restore_after_fork(bool restore);


// [[Rcpp::init]]
int avoid_openmp_hang_within_fork(DllInfo *dll);


#endif  // OPEN_MP_H
