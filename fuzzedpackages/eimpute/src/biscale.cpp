#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
using namespace std;
// [[Rcpp::export]]
Eigen::MatrixXd vec2mat(Eigen::VectorXd x, int type, int num){
  // type = 1 denote column, type = 2 denote row
  int veclen = x.size();
  Eigen::MatrixXd Y;
  if(type == 1){
    Y.setZero(veclen, num);
    for(int i = 0; i < num; i++){
      Y.col(i) = x;
    }
  }else{
    Y.setZero(num, veclen);
    for(int i = 0; i < num; i++){
      Y.row(i) = x;
    }
  }
  return Y;
}

// [[Rcpp::export]]
List biscale_alt(Eigen::MatrixXd x, Eigen::MatrixXd ind, Eigen::VectorXd obsrow, Eigen::VectorXd obscol, int max_it, double tol, Eigen::VectorXd alpha, Eigen::VectorXd beta, Eigen::VectorXd tau, Eigen::VectorXd gamma,
            bool row_mean, bool col_mean, bool row_std, bool col_std) {
  int m = x.rows();
  int n = x.cols();
  int obnum = ind.rows();
  double err = 10;
  double tempconst =  0;
  int iter = 0;
  Eigen::MatrixXd temp;
  temp.setZero(m,n);
  Eigen::VectorXd alphaold;
  Eigen::VectorXd betaold;
  Eigen::VectorXd tauold;
  Eigen::VectorXd gammaold;
  while (err > tol && iter < max_it){
    alphaold = alpha;
    betaold = beta;
    tauold = tau;
    gammaold = gamma;
    // beta
    if(col_mean == TRUE){
      tempconst = tau.array().inverse().sum();
      for(int i = 0; i < obnum; i++){
        temp(ind(i,0), ind(i, 1)) = (x(ind(i,0), ind(i, 1)) - alphaold(ind(i,0)))/tau(ind(i,0));
      }
      beta = temp.colwise().sum()/tempconst;
    }

    // alpha
    if(row_mean == TRUE){
      tempconst = gamma.array().inverse().sum();
      for(int i = 0; i < obnum; i++){
        temp(ind(i,0), ind(i, 1)) = (x(ind(i,0), ind(i, 1)) - betaold(ind(i,1)))/gamma(ind(i,1));
      }
      alpha = temp.rowwise().sum()/tempconst;
    }

    // gamma
    if(col_std == TRUE){
      for(int i = 0; i < obnum; i++){
        temp(ind(i,0), ind(i, 1)) = (x(ind(i,0), ind(i, 1)) - alphaold(ind(i,0)) - betaold(ind(i,1))) * (x(ind(i,0), ind(i, 1)) - alphaold(ind(i,0)) - betaold(ind(i,1)));
      }
      temp = temp.array() / vec2mat(tau.array().square(), 1, n).array();
      gamma = temp.colwise().sum();
      gamma = gamma.array()  / obscol.array();
      gamma = gamma.array().sqrt();
    }
    // tau
      if(row_std == TRUE){
        for(int i = 0; i < obnum; i++){
          temp(ind(i,0), ind(i, 1)) = (x(ind(i,0), ind(i, 1)) - alphaold(ind(i,0)) - betaold(ind(i,1))) * (x(ind(i,0), ind(i, 1)) - alphaold(ind(i,0)) - betaold(ind(i,1)));
        }
        temp = temp.array() / vec2mat(gamma.array().square(), 2, m).array();
        tau = temp.rowwise().sum();
        tau = tau.array() / obsrow.array();
        tau = tau.array().sqrt();
      }
    err = (alpha - alphaold).array().square().sum() > (beta - betaold).array().square().sum() ? (alpha - alphaold).array().square().sum() : (beta - betaold).array().square().sum();
    err = (tau - tauold).array().square().sum() > err ? (tau - tauold).array().square().sum() : err;
    err = (gamma - gammaold).array().square().sum() > err ? (gamma - gammaold).array().square().sum() : err;
    iter = iter + 1;

 }
  return Rcpp::List::create(alpha, beta, tau, gamma);
}




