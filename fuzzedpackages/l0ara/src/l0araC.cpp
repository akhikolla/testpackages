#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List l0araC(arma::mat x, arma::vec y, String family, double lam, int maxit, double eps) {
  // initialization
  int n = x.n_rows;
  int m = x.n_cols;
  arma::mat w = zeros(m, 1);
  arma::mat old_w(m, 1);
  arma::mat Xt(n, m);
  arma::mat xw(n, 1);
  arma::mat s1(n, 1);
  arma::mat A(n, n);
  arma::mat a = zeros(m, 1);
  arma::mat z(n, 1);

  int iter = 1;
  double ybar = mean(y);
  if(family=="gaussian"){
    w(0) = ybar;
  }
  if(family=="inv.gaussian"){
    w(0) = 1/pow(ybar, 2);
  }
  if(family=="gamma"){
    w(0) = 1/ybar;
  }
  if(family=="logit"){
    w(0) = log(ybar/(1-ybar));
  }
  if(family=="poisson"){
    w(0) = log(ybar);
  }
  Xt = x;
  while(iter < maxit) {
    old_w = w;
    xw = x*w;
    if(family=="gaussian"){
      s1 = xw;
      A = eye(n, n);
    }
    if(family=="poisson"){
      s1 = exp(xw);
      A = diagmat(vec(s1));
    }
    if(family=="logit"){
      s1 = 1/(1+exp(-xw));
      a = s1 % (1-s1);
      A = diagmat(vec(a));
    }
    if(family=="gamma"){
      s1 = 1/xw;
      a = -s1 % s1;
      A = diagmat(vec(a));
    }
    if(family=="inv.gaussian"){
      s1 = 1/sqrt(xw);
      a = -pow(s1, 3)/2;
      A = diagmat(vec(a));
    }
    z = A*xw + y-s1;
    if(n > m) {
      w = solve(trans(Xt)*A*x+lam*eye(m,m), trans(Xt)*z);
    } else {
      w = trans(Xt)*solve(A*x*trans(Xt)+lam*eye(n,n), z);
    }
    Xt = repmat(trans(w % w), n, 1) % x;
    iter += 1;

    // check for convergence
    if(iter >= maxit){
      warning("Did not converge. Increase maxit.");
    }
    if(norm(w-old_w, 2) < eps) {
      break;
    }
  }
  for(int i=0; i<m; i++){
    if(std::abs(w(i, 0)) < 1e-3){
      w(i,0) = 0;
    }
  }
  return List::create(Named("beta") = w, Named("iter") = iter);
}
