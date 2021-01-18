#include <cmath>
#include <iostream>
#include "binomial.h"

namespace lps {
  double Binomial::eval(const arma::colvec& beta)
  {
    checkDim(beta, p);
    eta = X * beta;
    arma::colvec tmp = -Y % eta + log(1 + exp(eta));
    return sum(tmp) / static_cast<double> (n);
  }

  void Binomial::gradient(arma::colvec& ret,
			  const arma::colvec& beta,
			  const arma::uvec& index) 
  {
    checkDim(beta, p);

    // evaluate gradient of the squre according to index
    arma::mat subX = subMatrix(X, index);
    eta = subX * beta.elem(index);
    ret = trans(subX) * (-Y + exp(eta) / (1 + exp(eta))) / static_cast<double> (n);
  }

  void Binomial::hessian(arma::mat& ret,
			 const arma::colvec& beta,
			 const arma::uvec& index) 
  {
    checkDim(beta, p);

    // evaluate Hessian matrix
    arma::mat subX = subMatrix(X, index);
    arma::colvec tmp = exp(eta) / pow(1 + exp(eta), 2);
    ret = trans(subX) * (repmat(tmp, 1, index.n_rows) % subX) 
      / static_cast<double> (n);
  }
}
