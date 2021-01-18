#include "RcppArmadillo.h"
#include "../inst/include/Countr_types.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

arma::Col <unsigned> seq0(unsigned to, unsigned by) {
  // unsigned from = 0;
  unsigned length = to + 1;

  Col <unsigned> res (length, fill::zeros);
  for (unsigned i = 0; i < length ; i ++)
    res(i) =  i*by;

  return(res);
}
