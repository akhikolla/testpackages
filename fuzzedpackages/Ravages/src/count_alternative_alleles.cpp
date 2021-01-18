#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"

using namespace Rcpp;
using namespace RcppParallel;

// caa = count alternative alleles
struct caa_p : public Worker {
  // input 
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol; 
  const size_t nrow;
  std::vector<bool> inverse;
  //output
  int * R;

 //constructeur
 caa_p(uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, std::vector<bool> inverse) 
            : data(data), ncol(ncol), true_ncol(true_ncol), nrow(nrow), inverse(inverse) {
    R = new int[ncol];
    std::fill(R, R+ncol, 0); 
  }

  //constructeur pour le split
  caa_p(caa_p & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), nrow(Q.nrow), inverse(Q.inverse) {
    R = new int[ncol];
    std::fill(R, R+ncol, 0); 
  }

  // destructeur
  ~caa_p() {
    delete [] R;
  }

  //worker
  void operator()(size_t beg, size_t end) {
    int gg[4];
    gg[3] = 0;
    for(size_t i = beg; i < end; i++) {
      if(inverse[i]) {
        gg[0] = 2; gg[1] = 1; gg[2] = 0;
      } else {
        gg[0] = 0; gg[1] = 1; gg[2] = 2;
      }
      size_t k = 0;
      for(size_t j = 0; j < true_ncol; j++) {
        uint8_t x = data[i][j];
        for(int ss = 0; ss < 4 && (4*j + ss) < ncol; ss++) {
          R[k++] += gg[x&3];
          x >>= 2;
        }
      }
    }
  }

  // join
  void join(const caa_p & Q) {
    std::transform(R, R + ncol, Q.R, R, std::plus<int>());
  }

};


//[[Rcpp::export]]
IntegerVector count_alternative_alleles(XPtr<matrix4> p_A, LogicalVector which_snps, LogicalVector inverse) {
  int n = p_A->nrow; // nb snps
  int m = p_A->ncol; // nb inds

  if(which_snps.length() != p_A->nrow)
    stop("Dimensions mismatch");

  int nb_snps = sum(which_snps);
  uint8_t ** data = new uint8_t * [nb_snps];
  std::vector<bool> inv(nb_snps);
  size_t k = 0;
  for(size_t i = 0; i < n; i++) {
    if(which_snps[i]) {
      inv[k] = inverse[i];
      data[k++] = p_A->data[i];
    }
  }

  caa_p X(data, p_A->ncol, p_A->true_ncol, p_A->nrow, inv);
  parallelReduce(0, nb_snps, X);
  delete [] data;

  IntegerVector R(m);
  std::copy(X.R, X.R+m, R.begin());

  return R;
}

RcppExport SEXP oz_count_alternative_alleles(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(count_alternative_alleles(p_A, which_snps, inverse));
    return rcpp_result_gen;
END_RCPP
}
