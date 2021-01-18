#ifndef BINOMIAL_H_
#define BINOMIAL_H_

#include "uniDistri.h"

namespace lps {
  // the class to implement univariate logistic regression
  class Binomial : public UniDistri {
  private:
    // eta function to be used in eval, gradient and hessian
    arma::colvec eta;
  public:
    explicit Binomial (const arma::mat& inputY,
		       const arma::mat& inputX):
    UniDistri(inputY, inputX) {};
    virtual double eval(const arma::colvec&);
    virtual void gradient(arma::colvec&,
			 const arma::colvec&,
			 const arma::uvec&);
    virtual void hessian(arma::mat&,
			 const arma::colvec&,
			 const arma::uvec&);
    virtual void mean(arma::colvec&) {};
    virtual void variance(arma::mat&) {};
    virtual arma::colvec getlinear() const {
      return eta;
    }
    virtual ~Binomial() {};
  };
}

#endif // BINOMIAL_H_

