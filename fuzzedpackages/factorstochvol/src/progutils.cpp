/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2016-2020
 *     Darjus Hosszejni Copyright (C) 2019-2020
 *  
 *  This file is part of the R package factorstochvol: Bayesian Estimation
 *  of (Sparse) Latent Factor Stochastic Volatility Models
 *  
 *  The R package factorstochvol is free software: you can redistribute
 *  it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or any
 *  later version of the License.
 *  
 *  The R package factorstochvol is distributed in the hope that it will
 *  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with the R package factorstochvol. If that is not the case,
 *  please refer to <http://www.gnu.org/licenses/>.
 */

#include "progutils.h"

double logdnormquot(double x, double y, double mu, double sigma) {
 return ((y-mu)*(y-mu) - (x-mu)*(x-mu)) / (2*sigma*sigma);
}

double logspecialquot(double x, double y, double alpha, double beta, double c) {
 return (alpha/c) * (x - y) - beta * (exp(x/c) - exp(y/c));
}


void store(const Rcpp::NumericMatrix &curfacload, Rcpp::NumericVector &facload,
           const Rcpp::NumericMatrix &curf,       Rcpp::NumericVector &f,
           const Rcpp::NumericMatrix &curh,       Rcpp::NumericVector &h,
           const Rcpp::NumericVector &curh0,      Rcpp::NumericMatrix &h0,
           const Rcpp::NumericMatrix &curpara,    Rcpp::NumericVector &para,
           const Rcpp::NumericVector &curlambda2, Rcpp::NumericMatrix &lambda2,
           const Rcpp::NumericMatrix &curtau2,    Rcpp::NumericVector &tau2,
           const Rcpp::NumericVector &curbeta, Rcpp::NumericVector &beta,
           const arma::umat &curmixind,  Rcpp::IntegerVector &mixind,
           const bool auxstore, const int thintime, const int where) {
 
 std::copy(curfacload.begin(), curfacload.end(), facload.begin() + where * curfacload.length());
 std::copy(curpara.begin(), curpara.end(), para.begin() + where * curpara.length());
  
 if (thintime == 1) { // store everything

  std::copy(curf.begin(), curf.end(), f.begin() + where * curf.length());
  std::copy(curh.begin(), curh.end(), h.begin() + where * curh.length());
 
 } else if (thintime == -1) { // store only t = T

  for (int i = 0; i < curf.nrow(); i++) {
   f(where*curf.nrow() + i) = curf(i, curf.ncol()-1);
  }
  
  for (int i = 0; i < curh.ncol(); i++) {
   h(where*curh.ncol() + i) = curh(curh.nrow()-1, i);
  }

 } else if (thintime > 1) { // store every thintimeth point in time
  
  int tmp = curf.ncol()/thintime;
  int tmpp = where * curf.nrow() * tmp;
  
  for (int j = 0; j < tmp; ++j) {
   int tmppp = j*thintime;
   int tmpppp = tmpp + j*curf.nrow();
   
   for (int i = 0; i < curf.nrow(); ++i) {
    f(tmpppp + i) = curf(i, tmppp);
   }
  }

  tmpp = where * curh.ncol() * tmp;
  
  for (int i = 0; i < curh.ncol(); ++i) {
   int tmpppp = tmpp + i*tmp; 
   
   for (int j = 0; j < tmp; ++j) {
    h(tmpppp + j) = curh(j*thintime, i);
   }
  }
 }
 
 std::copy(curh0.begin(), curh0.end(), h0.begin() + where * curh0.length());
 if (beta.size() > 0) {
   std::copy(curbeta.begin(), curbeta.end(), beta.begin() + where * curbeta.length());
 }
 
 if (auxstore) { // store mixture probabilities, mixture indicators, shrinkage hyperparas, h0
  std::copy(curmixind.begin(), curmixind.end(), mixind.begin() + where * curmixind.n_elem);
  std::copy(curlambda2.begin(), curlambda2.end(), lambda2.begin() + where * curlambda2.length());
  std::copy(curtau2.begin(), curtau2.end(), tau2.begin() + where * curtau2.length());
 }
}
