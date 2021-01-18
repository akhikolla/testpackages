#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"

using namespace Rcpp;

XPtr<matrix4> random_bed_matrix(NumericMatrix maf, NumericVector size) {
  int nb_pop = maf.nrow();
  if(size.length() != nb_pop)
    stop("Dimensions mismatch");
  int nrow = maf.ncol(); // nb snp
  int ncol = sum(size);  // nb inds
  XPtr<matrix4> p_A(new matrix4(nrow, ncol));
  for(int i = 0; i < nrow; i++) {
    int k = 0;
    for(int pop = 0; pop < nb_pop; pop++) {
      double q = maf(pop, i);
      double pr1 = (1-q)*(1-q), pr2 = pr1 + 2*(1-q)*q;
      int s = size[pop];
      for(int j = 0; j < s; j++) {
        double r = Rf_runif(0.0, 1.0);
        if(r < pr1) { 
          p_A->set(i,k++,0);
        } else if(r < pr2) {
          p_A->set(i,k++,1);
        } else {
          p_A->set(i,k++,2);
        }
      }
    }
  }
  return p_A;
}

XPtr<matrix4> new_bed_matrix(int nrow, int ncol) {
  XPtr<matrix4> p_A(new matrix4(nrow, ncol));
  return p_A;
}

void random_filling_bed_matrix(XPtr<matrix4> p_A, NumericMatrix maf, NumericVector size, int firstrow) {  
  int nb_pop = maf.nrow();
  if(size.length() != nb_pop)
    stop("Dimensions mismatch");
  int nrow = maf.ncol(); // nb snp
  int ncol = sum(size);  // nb inds
  if(ncol != p_A->ncol || nrow + firstrow > p_A-> nrow) { 
    stop("Dimensions mismatch");
  }
  for(int i = 0; i < nrow; i++) {
    int k = 0;
    int row = firstrow + i;
    for(int pop = 0; pop < nb_pop; pop++) {
      double q = maf(pop, i);
      double pr1 = (1-q)*(1-q), pr2 = pr1 + 2*(1-q)*q;
      int s = size[pop];
      for(int j = 0; j < s; j++) {
        double r = Rf_runif(0.0, 1.0);
        if(r < pr1) { 
          p_A->set(row,k++,0);
        } else if(r < pr2) {
          p_A->set(row,k++,1);
        } else {
          p_A->set(row,k++,2);
        }
      }
    }
  }
}


//[[Rcpp::export]]
void random_filling_bed_matrix_noHW(XPtr<matrix4> p_A, NumericMatrix proba_g0, NumericMatrix proba_g1, NumericVector size, int firstrow) {  
  int nb_pop = size.length();
  int nrow = proba_g0.ncol(); // nb snp
  int ncol = sum(size);  // nb inds
  if(proba_g0.nrow() != nb_pop || proba_g1.nrow() != nb_pop || nrow != proba_g1.ncol())
    stop("Dimensions mismatch");
  if(ncol != p_A->ncol || nrow + firstrow > p_A-> nrow) { 
    stop("Dimensions mismatch");
  }

  // construction matrice proba_g0 + proba_g1
  NumericMatrix pp2( clone(proba_g1) );
  for(int i = 0; i < nb_pop; i++) {
    for(int j = 0; j < nrow; j++) {
      pp2(i,j) += proba_g0(i,j);
      if( pp2(i,j) > 1.0 + 1e-14 ) 
        stop("Sum of probas > 1!");
    }
  }

  for(int i = 0; i < nrow; i++) {
    int k = 0;
    int row = firstrow + i;
    for(int pop = 0; pop < nb_pop; pop++) {
      double pr1 = proba_g0(pop,i); 
      double pr2 = pp2(pop,i);
      int s = size[pop];
      for(int j = 0; j < s; j++) {
        double r = Rf_runif(0.0, 1.0);
        if(r < pr1) { 
          p_A->set(row,k++,0);
        } else if(r < pr2) {
          p_A->set(row,k++,1);
        } else {
          p_A->set(row,k++,2);
        }
      }
    }
  }
}






RcppExport SEXP oz_random_bed_matrix(SEXP mafSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(random_bed_matrix(maf, size));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP oz_new_bed_matrix(SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(new_bed_matrix(nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP oz_random_filling_bed_matrix(SEXP p_ASEXP, SEXP mafSEXP, SEXP sizeSEXP, SEXP firstrowSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type firstrow(firstrowSEXP);
    random_filling_bed_matrix(p_A, maf, size, firstrow);
    return R_NilValue;
END_RCPP
}

RcppExport SEXP oz_random_filling_bed_matrix_noHW(SEXP p_ASEXP, SEXP proba_g0SEXP, SEXP proba_g1SEXP, SEXP sizeSEXP, SEXP firstrowSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type proba_g0(proba_g0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type proba_g1(proba_g1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type firstrow(firstrowSEXP);
    random_filling_bed_matrix_noHW(p_A, proba_g0, proba_g1, size, firstrow);
    return R_NilValue;
END_RCPP
}


