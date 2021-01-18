#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
double COX(NumericVector tps, NumericVector cureE, NumericVector gene, NumericVector delta)
 {
  int n = tps.size();
  NumericVector S0(n), S1(n), xx = cureE*gene;

  double ss = 0;
  int i, j = n-1;
  while(j >= 0) {
    ss += cureE(j);
    for(i = j-1; i >= 0 && tps(i) == tps(j); i--) ss+= cureE(i);
    for(int k = i+1; k <= j; k++) S0(k) = ss;
    j = i;
  }

  ss = 0;
  j = n-1;
  while(j >= 0) {
    ss += xx(j);
    for(i = j-1; i >= 0 && tps(i) == tps(j); i--) ss+= xx(i);
    for(int k = i+1; k <= j; k++) S1(k) = ss;
    j = i;
  }

  NumericVector Kt0 = cumsum(delta/S0), Kt0minus(n);
  double der = Kt0(n-1); // = Theta0
  for(i = n-2; i >= 0 && Kt0(i) == der; i--);
  double pen = Kt0(i);
  for(i = 1; i < n; i++) Kt0minus(i) = (Kt0(i-1) == der)?pen:Kt0(i-1);

  NumericVector poids(n);
  for(i = 0; i < n; i++) poids(i) = (Kt0minus(i) == der)?0:(1+log1p(-Kt0minus(i)/der)); // le test doit être rendu inutile par la manip précédente

  double csm = 0, csn = 0;
  double W = 0, Sig = 0;
  for(i = 0; i < n; i++) {
    csm += delta(i)*poids(i)/S0(i);
    csn += delta(i)*poids(i)/S0(i)*S1(i)/S0(i);
    double a = delta(i)*poids(i)*(gene(i) - S1(i)/S0(i)) - csm*gene(i) + csn;
    W += a;
    Sig += a*a;
  }

  return W*W/Sig;
}

RcppExport SEXP iBST_COX(SEXP tpsSEXP, SEXP cureESEXP, SEXP geneSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type tps(tpsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cureE(cureESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gene(geneSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    __result = Rcpp::wrap(COX(tps, cureE, gene, delta));
    return __result;
END_RCPP
}

