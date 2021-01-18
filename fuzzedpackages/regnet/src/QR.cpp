#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"QR.h"

using namespace Rcpp;
using namespace arma;


arma::vec QRWMR(const arma::mat& x, const arma::vec& y, arma::vec b)
{
  int p = x.n_cols, n = x.n_rows;
  arma::vec t = y - x * b;
  arma::vec u(n, fill::none), w(n, fill::none);
  for(int m = 0; m < p; m++){
    //std::cout << "m: " << m <<std::endl;
    t += x.col(m) * b(m);
    //arma::vec t = y - x * b + x.col(m) * b(m);
    u = t/x.col(m);
    w = arma::abs(x.col(m))/n;
    
    uvec index = sort_index(u);
    arma::vec _w = w(index), _u = u(index);

    double TotalWeight = accu(w), SUM = 0;
    int j = 0;
    do{
      SUM += _w(j)/TotalWeight;
      j++;
    }while(SUM <= 0.5);
    b(m) = _u(j-1);
    t -= x.col(m) * b(m);
  }
  return(b);
}
