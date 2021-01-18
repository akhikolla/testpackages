#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using std::log;
using std::exp;
using std::max;
using std::abs;
using std::sqrt;
using std::pow;

using namespace Rcpp; 




// **********************************************************//
//     	           Likelihood function            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull (arma::vec Y,
                       arma::vec Y0,
                       arma::vec eXB, 
                       arma::vec delta,
                       arma::vec C,
                       arma::vec LY,
                       double rho) {
  
  double r = exp(-rho);
  
  arma::vec lambda = eXB;
  arma::uvec ids1 = find(lambda == 0);
  lambda.elem(ids1).fill(exp(-740));
  arma::uvec ids2 = find(lambda == arma::datum::inf);
  lambda.elem(ids2).fill(exp(700));
  
  arma::vec eXBY = eXB % Y;
  arma::uvec ids3 = find(eXBY == arma::datum::inf);
  eXBY.elem(ids3).fill(exp(700));
  arma::uvec ids4 = find(eXBY == 0);
  eXBY.elem(ids4).fill(exp(-740));
  
  arma::vec eXBYr = pow(eXBY, r);
  arma::uvec ids5 = find(eXBYr == arma::datum::inf);
  eXBYr.elem(ids5).fill(exp(700));  
  
  arma::vec eXBY0 = eXB % Y0;
  arma::uvec ids6 = find(eXBY0 == arma::datum::inf);
  eXBY0.elem(ids6).fill(exp(700));
  arma::uvec ids7 = find(eXBY0 == 0);
  eXBY0.elem(ids7).fill(exp(-740));
  
  arma::vec eXBY0r = pow(eXBY0, r);
  arma::uvec ids8 = find(eXBY0r == arma::datum::inf);
  eXBY0r.elem(ids8).fill(exp(700)); 
  arma::uvec ids9 = find(eXBY0r == 0);
  eXBY0r.elem(ids9).fill(exp(-740));
  
  arma::vec st = exp(-eXBYr);
  
  arma::vec st0 = exp(-eXBY0r);
  arma::uvec ids10= find(st0 == 0);
  st0.elem(ids10).fill(exp(-740));
  arma::uvec ids11 = find(st0 == arma::datum::inf);
  st0.elem(ids11).fill(exp(700));
  
  arma::vec S = st / st0;
  
  arma::vec lst0 = log(st0);
  arma::uvec ids18 = find(lst0 == -arma::datum::inf);
  lst0.elem(ids18).fill(-740);
  
  arma::vec f = log(r) + r * log(lambda) + (r - 1) * log(Y) + log(S);
  arma::uvec ids12 = find(f == arma::datum::inf);
  f.elem(ids12).fill(exp(700)); 
  arma::uvec ids13 = find(f == -arma::datum::inf);
  f.elem(ids13).fill(-740);
  
  arma::vec ldelta = log(1 - delta);
  arma::uvec ids14 = find(ldelta == -arma::datum::inf);
  ldelta.elem(ids14).fill(-740);
  arma::uvec ids15 = find(ldelta == arma::datum::inf);
  ldelta.elem(ids15).fill(exp(700)); ; 
  
  arma::vec sdelta = (1 - delta) % S;
  arma::uvec ids17 = find(sdelta == 0);
  sdelta.elem(ids17).fill(exp(-740));
  
  arma::vec cens(lambda.size());
  for(unsigned int i=0; i<cens.size(); ++i) {
    if ( LY[i]==0 || C[i]==0 ) {
      cens[i] = log(delta[i] + sdelta[i]);
    } else {
      cens[i] = 0;
    }
  }
  
  arma::vec nocens(lambda.size());
  for(unsigned int i=0; i<nocens.size(); ++i) {
    if ( C[i]==1 && LY[i]==1 ) {
      nocens[i] = ldelta[i] + f[i];
      arma::uvec ids16 = find(nocens == arma::datum::inf);
      nocens.elem(ids16).fill(exp(700)); 
    } else {
      nocens[i] = 0;
    }
  }
  
  arma::vec llik = nocens + cens;
  
  return sum(llik);
}
