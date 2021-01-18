#ifndef PENALTY_H_
#define PENALTY_H_

#include <RcppArmadillo.h>

namespace lps {
  class penalty {
  public:
    virtual double eval(const arma::colvec&,
			const arma::uvec&) const = 0;
    virtual ~penalty() {};
  };

  class l1 : public penalty {
  public:
    virtual double eval(const arma::colvec& beta,
			const arma::uvec& index) const
    {
      return arma::norm(beta.elem(index), 1);
    }
    virtual ~l1() {};
  };

}

#endif // PENALTY_H_
