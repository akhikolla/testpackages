#ifndef BIGSTATSR_COLSTATS_HPP_INCLUDED
#define BIGSTATSR_COLSTATS_HPP_INCLUDED

/******************************************************************************/

#include <Rcpp.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

namespace bigstatsr {

template <class C>
ListOf<NumericVector> bigcolvars(C macc) {

  size_t n = macc.nrow();
  size_t m = macc.ncol();

  NumericVector res(m), res2(m);
  double x, xSum, xxSum;
  size_t i, j;

  for (j = 0; j < m; j++) {
    xSum = xxSum = 0;
    for (i = 0; i < n; i++) {
      x = macc(i, j);
      xSum += x;
      xxSum += x*x;
    }
    res[j] = xxSum - xSum * xSum / n;
    res2[j] = xSum;
  }

  return List::create(_["sum"] = res2,
                      _["var"] = res/(n-1));
}

}

/******************************************************************************/

#endif // #ifndef BIGSTATSR_COLSTATS_HPP_INCLUDED
