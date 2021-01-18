#ifndef GME_H_
#define GME_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <utility>
//#include <Rcpp.h>
#include "mvbernoulli.h"

namespace lps {
  class gme {
  private:
    // data for response, fixed and random effect design matrices
    arma::mat Y, X;
    arma::mat augY;
    arma::uvec Z;
    // dimension of fixed and random
    unsigned n, p, q, K;
    // fixed and random coefficients
    arma::colvec beta, sigma;
    arma::mat bi;
    arma::colvec mean;
    arma::mat var;
    // iterations 
    unsigned iterSpent;
    // number of levels
    unsigned numLevel, numCol;
    int order;
    MVBernoulli* ptr;
    arma::uvec index;
    arma::uvec counts;

    // evaluation of gradient
    double eval(const arma::colvec&);
    void evalGrad(arma::colvec&, const arma::colvec&);
    void evalHess(arma::mat&, const arma::colvec&);
    
    // utility functions
    void gfun(arma::mat&);
    void g1(arma::mat&);
    void g2(arma::mat&);

    void std(const arma::mat&, arma::colvec&);
  public:
    // constructor
    explicit gme(const arma::mat& inputY,
		 const arma::mat& inputX,
		 const arma::uvec& inputZ) :
    Y(inputY), X(inputX), Z(inputZ), 
      n(inputX.n_rows), K(inputY.n_cols) {
      numLevel = static_cast<unsigned> (Z.max()) + 1;
      counts = arma::zeros <arma::umat> (numLevel, 1);
      for (unsigned i = 0; i < numLevel; i++) {
	arma::uvec tmp = find(Z == i);
	counts(i) = tmp.n_rows;
      }
      ptr = new MVBernoulli(Y, X);
      ptr -> setLevel(numLevel, Z);
    }

    void setOrder(int);
    int fit ();

    arma::colvec& getCoef() {return beta;};
    unsigned getIter() {return iterSpent;};
    arma::colvec& getSigma() {return sigma;};
    arma::mat& getB() { return bi;}; 
  };
}

#endif // GME_H_
