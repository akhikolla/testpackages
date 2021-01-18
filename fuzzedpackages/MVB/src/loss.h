#ifndef LOSS_H_
#define LOSS_H_

#include <RcppArmadillo.h>

namespace lps {
  // loss class as an abstract class
  class Loss {
  protected:
    // design matrix
    const arma::mat& X;
    const unsigned p, n;
    // function to be used to raise exception when dimension does not match
    static void checkDim(const arma::colvec vec, unsigned dim)
    {
      if (vec.n_rows != dim)
	throw std::logic_error ("Dimension does not agree!");
    }

    // function to generate submatrix for given index
    static arma::mat subMatrix(const arma::mat& input,
			       const arma::uvec& index) {
      arma::mat ret(input.n_rows, index.n_rows);
      for (unsigned i = 0; i < ret.n_cols; i++)
	ret.col(i) = input.col(index(i));
      return ret;
    }

  public:
    // constructor for read in data
    explicit Loss(const arma::mat& input) : X(input), p(input.n_cols),
      n(input.n_rows) {};
    // evaluation of loss
    virtual double eval(const arma::colvec&) = 0;
    virtual void gradient(arma::colvec&,
			  const arma::colvec&,
			  const arma::uvec&) = 0;
    virtual void hessian(arma::mat&,
			 const arma::colvec&,
			 const arma::uvec&) = 0;
    virtual void mean(arma::colvec&) = 0;
    virtual void variance(arma::mat&) = 0;
    // return linear predictor
    virtual arma::colvec getlinear() const = 0;
    // add random perturbation to Y
    virtual void addRand(const arma::colvec&) = 0;
    // CORRECT ME: Any other way for me to specify order for MVB?
    virtual void toSetOrder(const unsigned) {};
    virtual ~Loss() {};
    // get the dimension of the unknow coefficients
    virtual unsigned getDim() const = 0;
    // get the number of augmented columns in response
    virtual unsigned getNumCol() const = 0;
    // get the number of columns in response
    virtual unsigned getK() const = 0;
    // get vectorized Y
    virtual arma::colvec getY() const = 0;
  };
}

#endif // LOSS_H_
