#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MI_Cpp(const NumericVector& X) {
  NumericVector Xs = clone(X).sort();
  double N = X.size();
  double med_Xs = median(Xs);
  double D = std::pow(N, 0.2) * (Xs[N/2 + 0.5 * std::pow(N,0.8) - 1] -
                 Xs[N/2 - 0.5 * std::pow(N, 0.8)]);
  NumericVector Xsub = Xs[Xs <= med_Xs];
  double g = mean(Xs) - 2/N*sum(Xsub);
  double S = 4 * var(Xs) + D*D - 4 * D * g;
  return sqrt(N) * 2 * (mean(Xs)-med_Xs) / sqrt(S);
}
