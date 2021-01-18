#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "allelecounter.h"

using namespace Rcpp;
using namespace RcppParallel;

//constructeur
allelecounter::allelecounter(const uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, const size_t nlevels, std::vector<int> group) 
          : data(data), ncol(ncol), true_ncol(true_ncol), nrow(nrow), nlevels(nlevels), group(group) {
  R = new int[2*nlevels*nrow];
  std::fill(R, R+2*nlevels*nrow, 0);
}

//constructeur pour le split
allelecounter::allelecounter(allelecounter & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), nrow(Q.nrow), nlevels(Q.nlevels), group(Q.group) {
  R = new int[2*nlevels*nrow];
  std::fill(R, R+2*nlevels*nrow, 0); 
}

// destructeur
allelecounter::~allelecounter() {
  delete [] R;
}

//worker
void allelecounter::operator()(size_t beg, size_t end) {
  for(size_t i = beg; i < end; i++) {
    for(size_t j = 0; j < true_ncol; j++) {
      uint8_t x = data[i][j];
      for(int ss = 0; ss < 4 && (4*j + ss) < ncol; ss++) {
        if( (x&3) != 3 ) {
          R[ 2*(i*nlevels + group[4*j+ss]-1) ] += (x&3);
          R[ 2*(i*nlevels + group[4*j+ss]-1) + 1] += 2-(x&3);
        }
        x >>= 2;
      }
    }
  }
}

// join
void allelecounter::join(const allelecounter & Q) {
  std::transform(R, R + 2*nlevels*nrow, Q.R, R, std::plus<int>());
}



