// ----------------------------------------------------------------------
// This file is part of the R package bayesImageS. It contains utility
// functions common to both Metropolis-Hastings and sequential Monte
// Carlo algorithms for the hidden Potts model.
// Copyright (C) 2013-2015  Matthew Moores
//
// bayesImageS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bayesImageS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------
#ifndef _bayesImageS_POTTS_UTIL_H
#define _bayesImageS_POTTS_UTIL_H

#include <RcppArmadillo.h>

arma::uvec unsign(const Rcpp::IntegerVector & x);

arma::umat unsignMx(const Rcpp::IntegerMatrix & m);

arma::rowvec rgamma(const arma::rowvec & shape, const arma::rowvec & rate);

arma::rowvec rnorm(const arma::rowvec & mean, const arma::rowvec & stddev);

arma::mat dnorm(const Rcpp::NumericVector & yunique, const arma::uvec & ymatch, const arma::rowvec & mean, const arma::rowvec & stddev);

arma::umat randomIndices(const unsigned n, int k);

unsigned sum_ident(const arma::umat & z, const arma::umat & neigh, const std::vector<arma::uvec> & blocks);

double sum_logs(arma::vec log_prob);

double interp(double val, unsigned idx, const arma::mat & path);

arma::rowvec gibbsMeans(const arma::rowvec & nZ, const arma::rowvec & sumY,
                        const arma::rowvec & pr_mu, const arma::rowvec & pr_mu_tau,
                        const arma::rowvec & sigma);

arma::rowvec gibbsStdDev(const arma::rowvec & nZ, const arma::rowvec & sumY,
                         const arma::rowvec & sqDiff, const arma::rowvec & pr_sd_nu,
                         const arma::rowvec & pr_sd_SS, const arma::rowvec & mean);

arma::rowvec gibbsDirichlet(const arma::rowvec & nZ, const arma::rowvec & pr_lambda);

void gibbsLabels(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                 arma::umat & z, arma::umat & alloc, const double beta,
                 const arma::mat & log_xfield);

void gibbsLabelsNoData(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                 arma::umat & z, arma::umat & alloc, const double beta);

void swLabelsNoData(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                    const double beta, const double k, arma::umat & z, arma::umat & alloc);

// Computes the number of neighbouring pixels allocated to component j, for pixel i. 
void neighbj(arma::mat & ne, arma::uvec & e, const arma::umat & z, const arma::umat & neigh);

// pseudo log-likelihood of the Potts model
double pseudolike(const arma::mat & ne, const arma::uvec & e, const double b, const unsigned n, const unsigned k);

void classify(arma::umat & z, arma::umat & alloc, const arma::rowvec & lambda, const arma::mat & log_xfield);

void updateStats(const arma::colvec & y, const arma::umat & z,
                 arma::rowvec & nZ, arma::rowvec & sumY, arma::rowvec & sqDiff);
                 
double rwmh(const double mean, const double stddev, const double prior[2]);

#endif
