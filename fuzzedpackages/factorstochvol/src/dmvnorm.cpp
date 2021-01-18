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

#include <RcppArmadillo.h>

using namespace Rcpp;

RcppExport SEXP dmvnorm(const SEXP x_in, const SEXP means_in, const SEXP vars_in, const SEXP log_in) {

 // note: SEXP to Rcpp conversion REUSES memory unless "clone"d
 // Rcpp to Armadillo conversion allocates NEW memory unless deact'd
 
 const bool loga = as<bool>(log_in);

 NumericMatrix x_(x_in);
 arma::mat x(x_.begin(), x_.nrow(), x_.ncol(), false);

 NumericMatrix means_(means_in);
 arma::mat means(means_.begin(), means_.nrow(), means_.ncol(), false);
 
 NumericVector vars_(vars_in);
 IntegerVector vars_dims = vars_.attr("dim");
 int dim = vars_dims(0);
 int reps = vars_dims(2);
 arma::cube vars(vars_.begin(), dim, dim, reps, false);

 NumericVector out_(reps);
 arma::vec out(out_.begin(), reps, false);

 arma::mat root(dim, dim);
 arma::vec tmp(dim);
 double tracelog;
 
 const double normalizer = -(dim/2.) * log(2. * M_PI);
 const arma::mat I = arma::eye<arma::mat>(dim, dim);
 
 for (int i = 0; i < reps; i++) {
  try {
   root = arma::trans(arma::solve(arma::trimatu(arma::chol(vars.slice(i))), I));
   tracelog = arma::sum(arma::log(root.diag()));
   tmp = root * (x.col(i) - means.col(i));
   out(i) = normalizer - .5 * arma::sum(tmp%tmp) + tracelog;
  } catch (...) {
   out(i) = R_NegInf;
  }
 }

 if (!loga) out_ = exp(out_);
 return wrap(out_);
}
