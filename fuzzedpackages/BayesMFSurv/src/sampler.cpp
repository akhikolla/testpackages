#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using std::log;
using std::exp;
using std::max;
using std::abs;
using std::sqrt;
using std::pow;

using namespace Rcpp;



// **********************************************************//
//     	           Likelihood function            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull (arma::vec Y,
                    arma::vec eXB,
                    arma::vec alpha,
                    arma::vec C,
                    double lambda) {
  arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
  arma::vec dexp2 = pow(eXB % Y, lambda - 1);
  arma::vec dexp3 = pow(eXB % Y, lambda);
  arma::vec llik1 = (1 - alpha) % dexp1 + lambda * alpha % eXB % dexp2 % dexp1;
  arma::vec lalpha = log(1 - alpha);

  arma::uvec ids0 = find(llik1 == 0);
  llik1.elem(ids0).fill(exp(-740));
  arma::uvec ids1 = find(llik1 == arma::datum::inf);
  llik1.elem(ids1).fill(exp(700));
  arma::uvec ids2 = find(dexp3 == arma::datum::inf);
  dexp3.elem(ids2).fill(exp(700));
  arma::uvec ids3 = find(lalpha == -arma::datum::inf);
  lalpha.elem(ids3).fill(-740);
  arma::vec llik = C % log(llik1) + (1 - C) % (lalpha - dexp3);
  return sum(llik);
}


// **********************************************************//
//     	             Likelihood function                     //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull2(arma::vec Y,
                    arma::vec Y0,
					          arma::vec eXB,
				          	arma::vec alpha,
				          	arma::vec C,
					          double lambda) {
  arma::vec eXBY = eXB % Y;
  arma::uvec ids = find(eXBY == 0);
  eXBY.elem(ids).fill(exp(-740));
  arma::uvec ids3 = find(alpha == 0);
  alpha.elem(ids3).fill(exp(-740));
  arma::vec dexp0 = pow(eXB % Y0, lambda);
  arma::vec dexp4 = exp(dexp0);
  arma::uvec ids5 = find(dexp0 == arma::datum::inf);
  dexp0.elem(ids5).fill(exp(700));
  arma::uvec ids4 = find(dexp4 == arma::datum::inf);
  dexp4.elem(ids4).fill(exp(700));
  arma::vec dexp1 = exp(-pow(eXBY, lambda));
	arma::vec dexp2 = pow(eXBY, lambda - 1);
	arma::uvec ids6 = find(dexp2 == arma::datum::inf);
	dexp2.elem(ids6).fill(exp(700));
	arma::vec dexp3 = pow(eXBY, lambda);
	arma::vec llik1 = (1 - alpha) + lambda * alpha % eXB % dexp2 % dexp1 % dexp4;
	arma::uvec ids1 = find(llik1 == arma::datum::inf);
	llik1.elem(ids1).fill(exp(700));
	arma::uvec ids0 = find(llik1 == 0);
	llik1.elem(ids0).fill(exp(-740));
	arma::uvec ids2 = find(dexp3 == arma::datum::inf);
	dexp3.elem(ids2).fill(exp(700));

	arma::vec llik = C % log(llik1) + (1 - C) % (log(alpha) - dexp3 + dexp0);
	return sum(llik);
}

