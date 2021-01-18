#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"

using namespace Rcpp;
using namespace RcppParallel;

// 
struct ploc : public Worker {
  // input 
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol; 
  const size_t nrow;
  const size_t nlevels;
  std::vector<int> group; // facteur à nlevels niveaux
  std::vector<bool> inverse;
  //output
  int * R;

 //constructeur
 ploc(uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, const size_t nlevels, std::vector<int> group, std::vector<bool> inverse) 
            : data(data), ncol(ncol), true_ncol(true_ncol), nrow(nrow), nlevels(nlevels), group(group), inverse(inverse) {
    R = new int[nlevels*nrow];
    std::fill(R, R+nlevels*nrow, 0); 
  }

  //constructeur pour le split
  ploc(ploc & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), nrow(Q.nrow), nlevels(Q.nlevels), group(Q.group), inverse(Q.inverse) {
    R = new int[nlevels*nrow];
    std::fill(R, R+nlevels*nrow, 0); 
  }

  // destructeur
  ~ploc() {
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
      for(size_t j = 0; j < true_ncol; j++) {
        uint8_t x = data[i][j];
        for(int ss = 0; ss < 4 && (4*j + ss) < ncol; ss++) {
          R[ i*nlevels + group[4*j+ss]-1 ] += gg[x&3];
          x >>= 2;
        }
      }
    }
  }

  // join
  void join(const ploc & Q) {
    std::transform(R, R + nlevels*nrow, Q.R, R, std::plus<int>());
  }

};


//[[Rcpp::export]]
IntegerMatrix alt_alleles_by_factor(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector group, LogicalVector inverse) {
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

  std::vector<int> gr(m);
  for(size_t i = 0; i < m; i++)
    gr[i] = group[i];

  int nlevels = as<CharacterVector>(group.attr("levels")).size();

  ploc X(data, p_A->ncol, p_A->true_ncol, p_A->nrow, nlevels, gr, inv);
  parallelReduce(0, nb_snps, X);
  delete [] data;

  IntegerMatrix R(nlevels,nb_snps);
  std::copy(X.R, X.R + nlevels*nb_snps, R.begin());

  return R;
}

RcppExport SEXP oz_alt_alleles_by_factor(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP grSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(grSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(alt_alleles_by_factor(p_A, which_snps, group, inverse));
    return rcpp_result_gen;
END_RCPP
}
