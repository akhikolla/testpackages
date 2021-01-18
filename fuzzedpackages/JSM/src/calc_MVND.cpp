#include <RcppEigen.h> 

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

double calc_MVND(const Eigen::Map<Eigen::VectorXd> & x, const Eigen::Map<Eigen::VectorXd> & mu, const Eigen::Map<Eigen::MatrixXd> & K){ 
 
  //  This function implements:
  //  dmvnorm(x, mu, K)

  double p = x.size();
  Eigen::LLT<Eigen::MatrixXd> LLT_of_K(K); // compute the Cholesky decomposition of K
  if ( !LLT_of_K.info() ) { 
    Eigen::MatrixXd Rooti = LLT_of_K.matrixLLT().triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(p,p)); 
    double quads = (Rooti * (x-mu)).array().square().sum();
    double res   =  exp( -((p/2.0)*log(2.0*M_PI))  + Rooti.diagonal().array().log().sum() - 0.5*quads); 
    return ( res );
  } else {
    return ( -1.0 );
  }
}

