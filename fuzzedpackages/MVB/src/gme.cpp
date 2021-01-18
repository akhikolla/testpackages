#include "gme.h"

namespace lps {
  // container to implement recent m iters for L-BFGS
  void gme::setOrder(int input) {
    order = input;
    ptr -> toSetOrder(order);
    p = ptr -> getDim();
    numCol = ptr -> getNumCol();
    beta = arma::zeros <arma::colvec> (p, 1);
    bi = arma::zeros<arma::mat> (numLevel, numCol);
    bi.randu();
    ptr -> setbi(bi);
    sigma = arma::ones <arma::colvec> (numCol, 1);
    index = arma::zeros<arma::uvec> (p, 1);
    for (unsigned i = 0; i < p; i++)
      index(i) = i;
    augY = ptr -> getAug();
  }

  int gme::fit () {
    unsigned iter;
    // fit with BFGS
    arma::mat invHess = arma::eye<arma::mat> (p, p);
    double alpha = 1.0;
    const double alphaIncr = 2.0;
    const double alphaDecr = .8;
    const double alphaMax = 1e8;
    //const double alphaMin = 1e-3;
    ptr -> eval(beta);

    for (iter = 0; iter < 500; iter++) {
      // alpha *= alphaDecr;
      alpha = 1;
      arma::colvec grad;
      evalGrad(grad, beta);
      if (norm(grad, 2) < 1e-3) break;
      ptr -> mean(mean);
      ptr -> variance(var);
      double current = eval(beta);
      // reevaluate bi and sigma
      for (unsigned i = 0; i < n; i++)
	for (unsigned j = 0; j < numCol; j++) {
	  unsigned pos = i * numCol + j;
	  bi(Z(i), j) += sigma(j) * sigma(j) * (augY(i, j) - mean(pos))
	    / static_cast<double> (n);
	}
      ptr -> setbi(bi);
      std(bi, sigma);

      // arma::colvec step = -invHess * grad;
      evalHess(invHess, beta);
      arma::colvec step = -solve(invHess, grad);
      arma::colvec coef0 = beta;
      // implement line search
      beta = coef0 + step / alpha;
      double val = eval(beta);
      /*
      while (val > current) {
	std::cout << "alpha = " << alpha << std::endl;
	if (alpha > alphaMax) return -1;
	alpha *= alphaIncr;
	beta = coef0 + step / alpha;
	val = eval(beta);
	std::cout << val << std::endl;
      }
      */
      arma::colvec grad0 = grad;
      evalGrad(grad, beta);
      // compute new invHess
      /*
      arma::colvec s = beta - coef0;
      arma::colvec y = grad - grad0;
      arma::colvec tmp = y.t() * s;
      double rho = 1 / tmp(0);
      std::cout << "rho = " << rho << std::endl;
      std::cout << s << y << std::endl;
      invHess = (arma::eye<arma::mat> (p, p) - rho * s * y.t()) * invHess *
       (arma::eye<arma::mat> (p, p) - rho * y * s.t()) + rho * s * s.t();
      */
      // reevaluate sigma
    }
    iterSpent = iter;
    return iterSpent;
  }

  double gme::eval(const arma::colvec& vec) {
    double ret = 0;
    arma::mat gg, gg2;
    gfun(gg);
    g2(gg2);
    for (unsigned i = 0; i < numLevel; i++)
      for (unsigned j = 0; j < numCol; j++)
	ret -= (2.0 * log(sigma(j)) + log(std::abs(gg2(i, j))) - gg(i, j)) /
	  static_cast<double> (n); 

    return ret;
  }

  void gme::evalGrad (arma::colvec& ret,
		      const arma::colvec& vec) {
    ptr -> gradient(ret, vec, index);
  }

  void gme::evalHess (arma::mat& ret,
		      const arma::colvec& vec) {
    ptr -> hessian(ret, vec, index);
  }
  void gme::gfun(arma::mat& ret) {
    double loss = ptr -> eval(beta);
    ret = -.5 * bi % bi;
    for (unsigned j = 0; j < numCol; j++)
      ret.col(j) = ret.col(j) / sigma(j) / sigma(j);
    ret = ret - loss * static_cast<double> (n);
  }

  void gme::g1(arma::mat& ret) {
    ret = arma::zeros <arma::mat> (numLevel, numCol);
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = 0; j < numCol; j++) {
	unsigned pos = i * numCol + j;
	ret(Z(i), j) += augY(i, j) - mean(pos);
      }
    for (unsigned i = 0; i < numLevel; i++)
      for (unsigned j = 0; j < numCol; j++)
	ret(i, j) -= bi(i, j) / sigma(j) / sigma(j);
  }

  void gme::g2(arma::mat& ret) {
    ret = arma::zeros <arma::mat> (numLevel, numCol);
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = 0; j < numCol; j++) {
	unsigned pos = i * numCol + j;
	ret(Z(i), j) -= var(pos, pos);
      }
    for (unsigned i = 0; i < numLevel; i++)
      for (unsigned j = 0; j < numCol; j ++)
	ret(i, j) -= 1.0 / sigma(j) / sigma(j);
    ret /= static_cast<double> (n);
  }

  void gme::std(const arma::mat& matrix,
		arma::colvec& ret) {
    for (unsigned i = 0; i < matrix.n_cols; i++) {
      double mm = 0;
      for (unsigned j = 0; j < matrix.n_rows; j++)
	mm += matrix(j, i);
      mm /= static_cast<double> (matrix.n_rows);
      mm = 0;
      double ss = 0;
      for (unsigned j = 0; j < matrix.n_rows; j++)
	ss += (matrix(j, i) - mm) * (matrix(j, i) - mm);
      ss /= static_cast<double> (matrix.n_rows);
      ret(i) = std::sqrt(ss);
    }
  }
	      
}
