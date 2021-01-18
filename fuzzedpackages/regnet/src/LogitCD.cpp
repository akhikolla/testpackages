#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"LogitCD.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;

arma::vec Network(const arma::mat& x, const arma::vec& y, double lam1, double lam2, arma::vec b, double r, const arma::mat& a, int n, int p)
{
  arma::vec y0 = x * b;
  arma::vec pi = 1/(1+exp(-y0)), t = (y -pi)*4;
  for(int k = 0; k < (p+1); k++){
    double bold = b(k);
    double l = arma::accu(x.col(k) % t)/n + b(k);
    int _p = p - 1;
    if(k == 0) b(k) = l;
    else{
      int m = k;
      if(k > _p) m = _p;
      double z = l*0.25 + lam2 * arma::as_scalar(a.row(k-1).subvec(m, _p) * b.subvec(m+1, p));
      double u = 0.25 + lam2 * arma::accu(abs(a.row(k-1).subvec(m, _p)));
      if(std::abs(z) > (r*lam1*u)){
        b(k) = z/u;
      }
      else{
        b(k) = Soft(z, lam1)/(u - 1/r);
      }
    }
    t -= x.col(k) * (b(k) - bold);
  }
  return(b);
}

arma::vec MCP(const arma::mat& x, const arma::vec& y, double lambda, arma::vec b, double r, int n, int p)
{
  arma::vec y0 = x * b;
  arma::vec pi = 1/(1+exp(-y0)), t = (y -pi)*4;
  for(int k = 0; k < (p+1); k++){
    double bold = b(k);
    double l = arma::accu(x.col(k) % t)/n + b(k);
    if(k == 0) b(k) = l;
    else if(std::abs(l) > (r*lambda)) b(k) = l;
    else{
      double z = l*0.25;
      b(k) = Soft(z, lambda)/(0.25 - 1/r);
    }
    t -= x.col(k) * (b(k) - bold);
  }
  return(b);
}

arma::vec Elastic(const arma::mat& x, const arma::vec& y, double lambda, arma::vec b, double alpha, int n, int p)
{
  arma::vec y0 = x * b;
  arma::vec pi = 1/(1+exp(-y0)), t = (y -pi)*4;

  for(int k = 0; k < (p+1); k++){
    double bold = b(k);
    double l = arma::accu(x.col(k) % t)/n + b(k);
    if(k == 0) b(k) = l;
    else{
      double z = l*0.25;
      b(k) = Soft(z, lambda*alpha)/(lambda*(1-alpha)+0.25);
    }
    t -= x.col(k) * (b(k) - bold);
  }
  return(b);
}
