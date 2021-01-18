#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"Robust.h"

using namespace Rcpp;
using namespace arma;

//Robust Network
arma::vec LadNet(const arma::mat& x, const arma::vec& y, double lam1, double lam2, arma::vec b, double r, const arma::mat& a, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    //std::cout << "m: " << m <<std::endl;
    t += x.col(m) * b(m);
    //arma::vec t = y - x * b + x.col(m) * b(m);
    arma::vec u(n+p-m, fill::zeros), w(n+p-m, fill::zeros);
    
    u.subvec(0,n-1) = t/x.col(m);
    u(n) = 0;
    if(m < p-1){
      u.subvec(n+1, n+p-m-1) = arma::sign(a.row(m).subvec(m+1, p-1)).t() % b.subvec(m+1, p-1);
    }
    //std::cout << "u: " << u.subvec(n+1, n+p-m-1).t() <<std::endl;
    
    w.subvec(0, n-1) = arma::abs(x.col(m))/n;
    if(std::abs(b(m)) < lam1 * r){
      w(n) = lam1 - std::abs(b(m))/r;
      //std::cout << "shrinked: "<<  lam1 << " - " << std::abs(b(m)) << "/" << r << " = " << w(n) << std::endl;
    }else{
      w(n) = 0;
    }
    if(m < p-1){
      w.subvec(n+1, n+p-m-1) = arma::abs(a.row(m).subvec(m+1, p-1)).t() * lam2;
    }
    //std::cout << "w: " << w.t() <<std::endl;
    uvec index = sort_index(u);
    arma::vec _w = w(index), _u = u(index);
    double TotalWeight = arma::accu(w), SUM = 0;
    int j = 0;
    do{
      SUM += _w(j)/TotalWeight;
      j++;
    }while(SUM <= 0.5);
    b(m) = _u(j-1);
    //std::cout << "b(m): " << b(m) <<std::endl;
    t -= x.col(m) * b(m);
  }
  return(b);
}

//Robust MCP
arma::vec LadMCP(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, double r, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    
    arma::vec u(n+1, fill::zeros), w(n+1, fill::zeros);
    
    u.subvec(0,n-1) = t/x.col(m);
    u(n) = 0;
    
    w.subvec(0, n-1) = arma::abs(x.col(m))/n;
    
    if(std::abs(b(m)) < lam1 * r){
      w(n) = lam1 - std::abs(b(m))/r;
    }else{
      w(n) = 0;
    }
    
    uvec index = sort_index(u);
    arma::vec _w = w(index), _u = u(index);
    double TotalWeight = arma::accu(w), SUM = 0;
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

//Robust Lasso
arma::vec LadLasso(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, int n, int p)
{
  arma::vec t = y - x * b;
  for(int m = 0; m < p; m++){
    t += x.col(m) * b(m);
    
    arma::vec u(n+1, fill::zeros), w(n+1, fill::zeros);
    
    u.subvec(0,n-1) = t/x.col(m);
    u(n) = 0;
    
    w.subvec(0, n-1) = arma::abs(x.col(m))/n;
    w(n) = lam1;
    
    uvec index = sort_index(u);
    arma::vec _w = w(index), _u = u(index);
    double TotalWeight = arma::accu(w), SUM = 0;
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
