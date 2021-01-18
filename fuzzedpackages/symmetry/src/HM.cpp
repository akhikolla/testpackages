#include <Rcpp.h>
#include <cstdint>
#include <cmath>
#include <helpers.h>
using namespace Rcpp;


// [[Rcpp::export]]
double HM_Cpp(const NumericVector& X, double a) {
  int n = X.size();
  NumericVector Xs = clone(X).sort();
  double TS_sum = 0;
  double i,j;
  double sqdiff;
  double wjk, logwjk;

  for (i = 0; i < n - 1; i++) {
    wjk = (i + 1) * (n - (i + 1)) / n / n * (i + 1) * (n - (i + 1)) / n / n;
    logwjk = log(wjk);
    TS_sum += - (2*a + logwjk) / (a*a * logwjk * logwjk);
    TS_sum -= - (2*a + logwjk) / ((Xs[i] + Xs[i])*(Xs[i] + Xs[i]) +
                                  a*a * logwjk * logwjk);

  }

  for(i = 0; i < n - 1; i++) {
    for(j = 0; j < i; j++) {
      // non diagonal are counted twice
      wjk = (i + 1) * (n - (i + 1)) / n / n * (j + 1) * (n - (j + 1)) / n / n;
      logwjk = log(wjk);
      TS_sum += - 2 * (2*a * logwjk) / ((Xs[i] - Xs[j])*(Xs[i] - Xs[j]) +
                                        a*a * logwjk * logwjk);
      TS_sum -= - 2 * (2*a * logwjk) / ((Xs[i] + Xs[j])*(Xs[i] + Xs[j]) +
                                        a*a * logwjk * logwjk);
    }
  }

  double TS_value = TS_sum / (2 * n * n);

  return TS_value;
}

