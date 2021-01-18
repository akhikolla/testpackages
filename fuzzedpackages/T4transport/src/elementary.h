#ifndef _T4transport_ELEMENTARY_H
#define _T4transport_ELEMENTARY_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

/* ELEMENTARY OPERAITONS
 * (01) compute_pdist2     : compute pairwise distance matrix
 * (02) cpp_subgrad_weight : Cuturi & Doucet (2014)'s notation of (1/lbd)*h(P)
 *      cpp_subgrad_plan
 *      cpp_subgrad_both     return a field of {weight} and {plan} in arma::mat
 */

arma::mat compute_pdist2(arma::mat& X, arma::mat& Y);
arma::vec cpp_subgrad_weight(arma::vec a, arma::vec b, arma::mat M, double lambda);
arma::mat cpp_subgrad_plan(arma::vec a, arma::vec b, arma::mat M, double lambda);
arma::field<arma::mat> cpp_subgrad_both(arma::vec a, arma::vec b, arma::mat M, double lambda);

#endif
