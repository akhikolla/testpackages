#include <RcppArmadillo.h>
using namespace arma;

RcppExport
SEXP kmeans(SEXP coords, SEXP centroids);
