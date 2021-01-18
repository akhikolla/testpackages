#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"NonRobust.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;

arma::vec LSNet(const arma::mat& x, const arma::vec& y, double lam1, double lam2, arma::vec b, double r, const arma::mat& a, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    //std::cout << "m: " << m <<std::endl;
    t += x.col(m) * b(m);
    double inp = arma::accu(square(x.col(m))), net1 = 0, net2 = 0;
    if(m < p-1){
      net1 = lam2 * as_scalar(a.row(m).subvec(m+1, p-1) * b.subvec(m+1, p-1));
      net2 = lam2 * arma::accu(arma::abs(a.row(m).subvec(m+1, p-1)));
    }
    double z = arma::accu(x.col(m) % t) + net1;
    double u = inp + net2;

    if(std::abs(z) > (r*lam1*u)){
      b(m) = z/u;
    }
    else{
      b(m) = Soft(z, lam1)/(u - 1/r);
    }
    t -= x.col(m) * b(m);
  }
  return(b);
}

arma::vec LSMCP(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, double r, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    //std::cout << "m: " << m <<std::endl;
    t += x.col(m) * b(m);
    double inp = arma::accu(square(x.col(m)));
    double z = arma::accu(x.col(m) % t);
    
    if(std::abs(z) > (r*lam1*inp)){
      b(m) = z/inp;
    }else{
      b(m) = Soft(z, lam1)/(inp - 1/r);
    }
    t -= x.col(m) * b(m);
  }
  return(b);
}

arma::vec LSLasso(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    //std::cout << "m: " << m <<std::endl;
    t += x.col(m) * b(m);
    double inp = arma::accu(square(x.col(m)));
    double z = arma::accu(x.col(m) % t);
    
    b(m) = Soft(z, lam1)/inp;
    t -= x.col(m) * b(m);
  }
  return(b);
}
