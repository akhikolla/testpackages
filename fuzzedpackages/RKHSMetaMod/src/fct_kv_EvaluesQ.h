#ifndef GramMatrix_H
#define GramMatrix_H
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
NumericVector k0linear(NumericVector u, double v);
NumericVector k0quad(NumericVector u, double v);
NumericVector int_1v(NumericVector x);
double int_1d(double x);
NumericVector k0matern(NumericVector u, double v);
NumericVector k0brownian(NumericVector u, double v);
NumericVector int_1gv(NumericVector x);
double int_1gd(double x);
NumericVector k0gaussian(NumericVector u, double v);
NumericVector k0(NumericVector u, double v, String kernel);
StringVector concatenate(StringVector a);
StringVector namesGrp(int d, int Dmax, List index);
SEXP calc_Kv(NumericMatrix X, String kernel, int Dmax, bool correction, bool verbose, double tol);
#endif
