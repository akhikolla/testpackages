#include <cmath>
#include "gaussian.h"

namespace lps {
  double Gaussian::eval(const arma::colvec& beta)
  {
    checkDim(beta, p);
    return pow(norm(Y - X * beta, 2), 2) / static_cast<double> (n);
  }

  void Gaussian::gradient(arma::colvec& ret,
			  const arma::colvec& beta,
			  const arma::uvec& index)
  {
    checkDim(beta, p);

    // evaluate gradient of the square according to index
    arma::mat subX = subMatrix(X, index);
    ret = -trans(subX) * (Y - subX * beta.elem(index)) / static_cast<double> (n);
  }

  void Gaussian::hessian(arma::mat& ret,
			 const arma::colvec& beta,
			 const arma::uvec& index)
  {
    checkDim(beta, p);

    // evaluate hessian matrix
    arma::mat subX = subMatrix(X, index);
    ret = trans(subX) * subX / static_cast<double> (n);
  }
}
