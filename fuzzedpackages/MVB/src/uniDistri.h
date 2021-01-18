#ifndef UNIDISTRI_H_
#define UNIDISTRI_H_

#include "loss.h"

namespace lps {
  class UniDistri : public Loss {
  protected:
    // response matrix (need a new copy because input is a matrix)
    arma::colvec Y;
  private:
  public:
    // always assume the data is store as response first
    explicit UniDistri (const arma::mat& inputY,
			const arma::mat& inputX) :
    Loss(inputX), Y(inputY.col(0)) {};
    virtual void addRand(const arma::colvec& err) {
      Y += err;
    }
    virtual unsigned getDim() const {return X.n_cols;};
    virtual unsigned getNumCol() const {return 1;};
    virtual unsigned getK() const {return 1; };
    virtual arma::colvec getY() const {return Y;};
  };
}

#endif // UNIDISTRI_H_

