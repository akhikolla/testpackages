#include <Rcpp.h>
using namespace Rcpp;

double CritEval(NumericMatrix X0, int q, int crit);
double MD2(NumericMatrix X0, int q);
double CD2(NumericMatrix X0, int q);
double WD2(NumericMatrix X0, int q);
double maximin(NumericMatrix X0, int q);
