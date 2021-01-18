#ifndef EMVS_H
#define EMVS_H

#include <math.h> 
#include <stdio.h>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

double beta(double x, double y);
double delogit(double x);
vec density_norm(vec &x, double mu, double sigma);
vec density_norm_log(vec &x, double mu, double sigma);
vec pnorm(vec x);
vec pnorm_log(vec x);
double erfc_log(double x);
double M_p(vec &post, double a, double b);
mat conj_E_beta_binom(vec &beta_k, double sigma_k, double v0, double v1, double theta, double t);
vec conj_M_beta(mat &XtY, mat &X, vec &Y, mat &XtX, vec &inv_var);
double conj_M_sigma(mat &Y, mat &X, vec &beta_k, vec &inv_var, double eta, double lambda);
double log_g(uvec &gamma, mat &X, mat &Y, double eta, double lambda, double v0, double v1, const string &type, double a, double b);
double log_prior(uvec& gamma, const string &type, double a, double b, int n);
mat ind_E_beta_binom(vec &beta_k, double v0, double v1, double theta, double t);
vec ind_M_beta(mat &XtY, mat &X, vec &Y, mat &XtX, vec &inv_var, double sigma_k);
double ind_M_sigma(mat &Y, mat &X, vec &beta_k, double eta, double lambda);

#endif
