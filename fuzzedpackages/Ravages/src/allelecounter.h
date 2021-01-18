#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>

#ifndef COUNTALLELES
#define COUNTALLELES

using namespace Rcpp;
using namespace RcppParallel;

struct allelecounter : public Worker {
  // input 
  const uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol; 
  const size_t nrow;
  const size_t nlevels;
  std::vector<int> group; // facteur à nlevels niveaux
  //output
  int * R;

  //constructeur
  allelecounter(const uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, const size_t nlevels, std::vector<int> group);
  //constructeur pour le split
  allelecounter(allelecounter & Q, Split);
  // destructeur
  ~allelecounter();
  //worker
  void operator()(size_t beg, size_t end);
  // join
  void join(const allelecounter & Q);
};

#endif
