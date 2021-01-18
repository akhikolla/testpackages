#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"LogitCD.h"
#include"RunLogistic.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec RunLogit(arma::mat& x, arma::vec& y, double lamb1, double lamb2, arma::vec b, double r, arma::mat& a, int p, 
                    double alpha, char method)
{
  int count = 0, n = x.n_rows;
  arma::vec bnew;
  while(count < 20){
    if(method == 'n'){
        bnew = Network(x, y, lamb1, lamb2, b, r, a, n, p);
    }else if(method == 'm'){
        bnew = MCP(x, y, lamb1, b, r, n, p);
    }else{
        bnew = Elastic(x, y, lamb1, b, alpha, n, p);
    }
    double dif = arma::accu(arma::abs(b - bnew))/n;
    if(dif < 0.001) break;
    else{
      b = bnew;
      count++;
    }
  }
  return(bnew);
}

// [[Rcpp::export]]
arma::vec RunNet(arma::mat& x, arma::vec& y, double lamb1, double lamb2, arma::vec b, double r, arma::mat& a, int p)
{
  int count = 0, n = x.n_rows;
  arma::vec bnew;
  while(count < 20){
    bnew = Network(x, y, lamb1, lamb2, b, r, a, n, p);
    double dif = arma::accu(arma::abs(b - bnew))/n;
    if(dif < 0.001) break;
    else{
      b = bnew;
      count++;
    }
  }
  return(bnew);
}

// [[Rcpp::export]]
arma::vec RunMCP(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double r, int p)
{
  int count = 0, n = x.n_rows;
  arma::vec bnew;
  while(count < 20){
    bnew = MCP(x, y, lambda, b, r, n, p);
    double dif = arma::accu(arma::abs(b - bnew))/n;
    if(dif < 0.001) break;
    else{
      b = bnew;
      count++;
    }
  }
  return(bnew);
}

// [[Rcpp::export]]
arma::vec RunElastic(arma::mat& x, arma::vec& y, double lambda, arma::vec b, double alpha, int p)
{
  int count = 0, n = x.n_rows;
  arma::vec bnew;
  while(count < 20){
    bnew = Elastic(x, y, lambda, b, alpha, n, p);
    double dif = arma::accu(arma::abs(b - bnew))/n;
    if(dif < 0.001) break;
    else{
      b = bnew;
      count++;
    }
  }
  return(bnew);
}
