#include <Rcpp.h>
#include <cmath>

#ifndef INCLUDED_SAMPLING_H
#define INCLUDED_SAMPLING_H

using namespace Rcpp;
using namespace std;

const double MINUS_LOG_SQRT_2_PI = -0.9189385;

inline double sq(double x) { return x * x; }

inline double log_likelihood_normal(
	double mean,
	double sd,
	double value
){ return MINUS_LOG_SQRT_2_PI - log(sd) - 0.5 * sq((value - mean)/sd); }

inline int random_integer(int exc_max)
{
	return (int)(unif_rand() * (double)exc_max) % exc_max;
}

inline double logit(double value) {
	return log(value) - log(1.0 - value);
}

inline double expit(double value) {
	return 1.0 - 1.0/(1.0+exp(value));
}

inline double logit_beta(
	double shape1,
	double shape2,
	double value
){ return value - 2.0 * log(1.0 + exp(value)) + log(R::dbeta(expit(value), shape1, shape2, true)); }

// [[Rcpp::export]]
List bevimed_mc(
	int its,
	LogicalVector y,
	IntegerVector var_block_start_index,
	IntegerVector var_block_stop_index,
	IntegerVector cases,
	IntegerVector counts,
	IntegerVector min_ac,
	double tau_shape1,
	double tau_shape2,
	double pi_shape1,
	double pi_shape2,
	double z_shape1,
	double z_shape2,
	LogicalMatrix z0,
	bool estimate_logit_z_rate,
	NumericVector logit_z_rates,
	NumericVector logit_z_rate_proposal_sds,
	NumericVector z_weights,
	bool estimate_phi,
	NumericVector log_phis,
	double log_phi_mean,
	double log_phi_sd,
	NumericVector log_phi_proposal_sds,
	NumericVector t,
	int swaps,
	bool annealing,
	int tandem_variant_updates,
	IntegerVector y1_case_block_start_index,
	IntegerVector y1_case_block_stop_index,
	IntegerVector y1_variants,
	bool return_z_trace,
	bool return_x_trace
);

#endif
