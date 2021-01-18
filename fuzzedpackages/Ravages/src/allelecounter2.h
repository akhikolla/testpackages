#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>

#ifndef COUNTALLELES2
#define COUNTALLELES2

using namespace Rcpp;
using namespace RcppParallel;

struct allelecounter2 : public Worker {
  // input 
  const uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol; 
  const size_t nrow;
  const size_t nlevels;
  std::vector<int> group; // facteur à nlevels niveaux
  std::vector<bool> inverse;
  //output
  int * R;

  //constructeur
  allelecounter2(const uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, const size_t nlevels, std::vector<int> group, std::vector<bool> inverse);
  //constructeur pour le split
  allelecounter2(allelecounter2 & Q, Split);
  // destructeur
  ~allelecounter2();
  //worker
  void operator()(size_t beg, size_t end);
  // join
  void join(const allelecounter2 & Q);

};

#endif
