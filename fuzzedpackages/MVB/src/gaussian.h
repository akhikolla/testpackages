#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include "uniDistri.h"

namespace lps {
  // class to implement least square loss
  class Gaussian : public UniDistri {
  public:
    explicit Gaussian (const arma::mat& inputY,
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
      return arma::zeros<arma::mat> (2, 1);
    }
    virtual ~Gaussian() {};
  };
}

#endif // GAUSSIAN_H_
