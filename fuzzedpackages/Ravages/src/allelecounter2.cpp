#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "allelecounter2.h"

using namespace Rcpp;
using namespace RcppParallel;

 //constructeur
 allelecounter2::allelecounter2(const uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, const size_t nlevels, std::vector<int> group, std::vector<bool> inverse) 
            : data(data), ncol(ncol), true_ncol(true_ncol), nrow(nrow), nlevels(nlevels), group(group), inverse(inverse) {
    R = new int[2*nlevels*nrow];
    std::fill(R, R+2*nlevels*nrow, 0); 
  }

  //constructeur pour le split
  allelecounter2::allelecounter2(allelecounter2 & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), nrow(Q.nrow), nlevels(Q.nlevels), group(Q.group), inverse(Q.inverse) {
    R = new int[2*nlevels*nrow];
    std::fill(R, R+2*nlevels*nrow, 0); 
  }

  // destructeur
  allelecounter2::~allelecounter2() {
    delete [] R;
  }

  //worker
  void allelecounter2::operator()(size_t beg, size_t end) {
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
          R[ 2*(i*nlevels + group[4*j+ss]-1) ] += gg[x&3];
          R[ 2*(i*nlevels + group[4*j+ss]-1) + 1] += gg[2-(x&3)];
          x >>= 2;
        }
      }
    }
  }

  // join
  void allelecounter2::join(const allelecounter2 & Q) {
    std::transform(R, R + 2*nlevels*nrow, Q.R, R, std::plus<int>());
  }




List alleles_by_factor(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector group, LogicalVector inverse) {
  int n = p_A->nrow; // nb snps
  int m = p_A->ncol; // nb inds

  if(which_snps.length() != p_A->nrow)
    stop("Dimensions mismatch");

  int nb_snps = sum(which_snps);
  int nlevels = as<CharacterVector>(group.attr("levels")).size();

  // une sécurité...  
  if(nb_snps == 0 || nlevels == 0) {
    IntegerMatrix R(nlevels, nb_snps);
    IntegerMatrix S(nlevels, nb_snps);
    List L;
    L["minor"] = R;
    L["major"] = S;
    return L;
  }
 

  // extraction des données pertinentes...
  const uint8_t ** data = new const uint8_t * [nb_snps];
  std::vector<bool> inv(nb_snps);
  size_t k = 0;
  for(size_t i = 0; i < n; i++) {
    if(which_snps[i]) {
      inv[k] = inverse[i];
      data[k++] = p_A->data[i];
    }
  }
  //******************************************

  std::vector<int> gr(m);
  for(size_t i = 0; i < m; i++)
    gr[i] = group[i];

  allelecounter2 X(data, p_A->ncol, p_A->true_ncol, nb_snps, nlevels, gr, inv);
  parallelReduce(0, nb_snps, X);
  delete [] data;

  IntegerMatrix R(nlevels,nb_snps);
  IntegerMatrix S(nlevels,nb_snps);
  for(size_t i = 0; i < nlevels*nb_snps; i++) {
    R[i] = X.R[2*i];
    S[i] = X.R[2*i+1];
  }
  
  List L;
  L["minor"] = R;
  L["major"] = S;
  return L;
}

RcppExport SEXP oz_alleles_by_factor(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP grSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(grSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(alleles_by_factor(p_A, which_snps, group, inverse));
    return rcpp_result_gen;
END_RCPP
}
