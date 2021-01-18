#ifndef LPS_H_
#define LPS_H_

// #include <memory>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <utility>
//#include <Rcpp.h> // for Rcout
#include "gaussian.h"
#include "binomial.h"
#include "mvbernoulli.h"
#include "distriFactory.h"
#include "penalty.h"

// tuning criteron
enum TUNING {AIC, BIC, GACV, BGACV};

namespace lps {
  // register possible distributions
  namespace {
    // Gaussian
    DistriHelper<Gaussian> registerGauss("gaussian");
    // Binomial
    DistriHelper<Binomial> registerBinom("binomial");
    // Multivariate Bernoulli
    DistriHelper<MVBernoulli> registerMVB("mvbernoulli");
  }

  class lps {
  private:
    // name of the distribution
    const std::string distri;
    // smart pointer to Loss and penalty object
    // NOTES: to used unique_ptr in C++11
    Loss* ptrLoss;
    penalty* ptrPenalty;
    // data
    arma::mat Y, X;
    // numrow for the number of parameters in mvb beta
    const unsigned n;
    unsigned numRow;
    // scores for the corresponding lambda for grid search
    arma::colvec scores;
    unsigned p; // not constant because this depends on distribution
    // results for fitted cofficients
    arma::mat beta;
    // parameters
    struct _param_ {
      bool output;        // wheter to print fitting process
      unsigned printIter; // print frequency
      double gptol;       // convergence tolerance
      unsigned maxIter;   // maximum number of iterations
      int tuneMethod;     // tuning method
      double sigma;       // fraction of inactive set to evaluate
      bool useNewton;     // whether to use newton step
      unsigned numRep;    // number of replicates
      double epSigma;     // std for randomized GACV
    } param;
    // indices to be penalized and not but estimated
    // penalized is set to be matrix to match dimension of lambda
    std::vector<arma::uvec> penalized;
    arma::uvec constants;
    // iteration results, number of iterations used
    arma::uvec iterSpent;

    // general printing
    void printIter(int iter, double gpnorm, bool print = 0) const {
      if (!param.output) return;
      if (iter == 0) {
	Rcpp::Rcout << "fit started" << std::endl;
      }
      if (iter % param.printIter == 0 || print) {
	Rcpp::Rcout << "iteration " << iter << " gpnorm = ";
	Rcpp::Rcout << gpnorm << std::endl;
      }
    }
    arma::mat epsilon;

    // print for lasso pattern search
    void print(const arma::colvec& lambda, unsigned iter, double gpnorm,
	       double alpha, unsigned numNonZero, bool converge = 0) const
    {
      if (iter == 1) {
	Rcpp::Rcout << std::endl;
	Rcpp::Rcout << "lambda = " << trans(lambda) << std::endl;
      } 
      if (iter % param.printIter == 0 || iter == 1 || converge) {
	Rcpp::Rcout << "iter " << iter << " gpnorm = ";
	Rcpp::Rcout << gpnorm << " nonzeros = " << numNonZero;
	Rcpp::Rcout << "(" << 100 * numNonZero / static_cast<double>(p);
	Rcpp::Rcout << "%)   alpha = " << alpha << std::endl;
      }
    }

    // initialization in constructors
    void construct() {
      ptrLoss = DistriFactory::instance().createLoss(distri, Y, X);
      p = ptrLoss -> getDim();
      numRow = p / ptrLoss -> getNumCol();
      ptrPenalty = new l1;
      setParam();
      epsilon.resize(n * ptrLoss -> getK(), param.numRep);
      epsilon.randu();
      epsilon *= param.epSigma;
    }
    // interface for fitting in lasso pattern search
    unsigned solveLPS(const arma::colvec&, const arma::colvec&, arma::colvec&);
    double objectiveFunc(const arma::colvec& beta_vec,
			 const arma::colvec& lambda) const {
      double ret = ptrLoss -> eval(beta_vec);
      for (unsigned i = 0; i < lambda.n_rows; i++) 
	ret += lambda(i) * ptrPenalty -> eval(beta_vec, penalized[i]);
      return ret;
    }
    double firstOrderMove(arma::colvec&, arma::colvec&,
			  bool&, unsigned, const arma::colvec&, double);
    double newtonStep(arma::colvec&, arma::colvec&, const arma::colvec&, double,
		      const arma::uvec&, const arma::uvec& );
    // merge two pre-sorted uvec
    arma::uvec merge(const arma::uvec vec1, const arma::uvec vec2) {
      arma::uvec ret = arma::zeros<arma::uvec> (vec1.n_rows + vec2.n_rows, 1);
      unsigned pos, pos1, pos2;
      pos = pos1 = pos2 = 0;
      while (pos1 != vec1.n_rows && pos2 != vec2.n_rows) {
	if (vec1(pos1) < vec2(pos2)) ret(pos++) = vec1(pos1++);
	else if (vec1(pos1) > vec2(pos2)) ret(pos++) = vec2(pos2++);
	else { ret(pos++) = vec1(pos1++); pos2++;}
      }
      while (pos1 != vec1.n_rows) ret(pos++) = vec1(pos1++);
      while (pos2 != vec2.n_rows) ret(pos++) = vec2(pos2++);
      ret.reshape(pos, 1);
      return ret;
    }

    // Newton Raphson algorithm to fit generalized linear model
    std::pair<arma::colvec, unsigned> newtonRaphson(bool, arma::uvec) const;

    // utility function for Nelder Mead
    double evalLambda(const arma::colvec&, arma::colvec&);
    // search nearest lambda for the given one
    int inline bestLambda(const arma::mat&, const arma::colvec&) const;

  public:
    // constructor for lps, accept input to Y and X
    explicit lps(const std::string& distribution,
		 const arma::mat& inputY,
		 const arma::mat& inputX) :
    distri(distribution), Y(inputY), X(inputX),
      n(inputX.n_rows) {
      construct();
      constants = arma::zeros<arma::umat>(p, 1);
      for (unsigned i = 0; i < p; i++)
	constants(i) = i;
    }

    // constructor with initial constants vector
    explicit lps (const std::string& distribution,
		  const arma::mat& inputY,
		  const arma::mat& inputX,
		  arma::uvec index) :
    distri(distribution), Y(inputY), X(inputX),
      n(inputX.n_rows){
      construct();
      constants = index;
    }

    // return iterations spent and scores
    const arma::uvec& getIter() const { return iterSpent; };
    const arma::colvec& getScore() const { return scores; };

    // set parameter values with default
    void setParam(bool output_ = 0,
		  unsigned printIter_ = 100,
		  double gptol_ = 1e-6,
		  unsigned maxIter_ = 500,
		  int tuneMethod_ = AIC,
		  double sigma_ = 0.1,
		  bool useNewton_ = 1,
		  unsigned numRep_ = 20,
		  double epSigma_ = 0.01) {
      param.output = output_;
      param.printIter = printIter_;
      param.gptol = gptol_;
      param.maxIter = maxIter_;
      param.tuneMethod = tuneMethod_;
      param.sigma = sigma_;
      param.useNewton = useNewton_;
      param.numRep = numRep_;
      param.epSigma = epSigma_;
    }

    arma::mat& getCoef() {return beta;};
    // set constant vector
    void setConst()
    {
      constants.reshape(Y.n_cols, 1);
      for (unsigned i = 0; i < constants.n_rows; i++)
	constants(i) = i * numRow;
      for (unsigned i = 0; i < ptrLoss -> getNumCol(); i++) {
	if (i < constants.n_rows) {
	  arma::uvec tmp = arma::zeros<arma::uvec> (ptrLoss -> getNumCol() - 1, 1);
	  for (unsigned j = 1; j < ptrLoss -> getNumCol(); j++)
	    tmp(j - 1) = j;
	  penalized.push_back(tmp);
	} else { 
	  arma::uvec tmp = arma::zeros<arma::uvec> (ptrLoss -> getNumCol(), 1);
	  for (unsigned j = 0; j < ptrLoss -> getNumCol(); j++)
	    tmp(j) = j;
	  penalized.push_back(tmp);
	}
      }
    }

    void setConst(const arma::uvec& index) 
    {
      constants = index;
      for (unsigned col = 0; col < ptrLoss -> getNumCol(); col++) {
	arma::uvec tmp = arma::zeros<arma::umat> (numRow, 1);
	unsigned pos = 0;
	for (unsigned i = col * numRow; i < col * numRow + numRow; i++) {
	  // search for appearance in constants (MAY need improve when constants is large
	  bool valid = 1;
	  for (unsigned j = 0; j < constants.n_rows; j++)
	    if (constants(j) == i) {
	      valid = 0;
	      break;
	    }
	  if (valid) tmp(pos++) = i;
	}
	tmp.reshape(pos, 1);
	penalized.push_back(tmp);
      } 

    }

    // to run newton as interface
    void runNewton() {
      beta = arma::zeros<arma::mat>(p, 1); 
      arma::uvec tmp = arma::zeros<arma::umat> (1, 1);
      std::pair<arma::colvec, unsigned> ret;
      ret = newtonRaphson(0, tmp);
      beta.col(0) = ret.first;
      iterSpent.reshape(1, 1);
      iterSpent(0) = ret.second;
    }
    // return dimension of the design matrix
    unsigned getDim() const {return ptrLoss -> getDim();};
    // set tuning method
    void setTune(int method) {param.tuneMethod = method;};
    // tuning with various criterion
    double tune(const arma::colvec&, int);
    // element-wise step fit
    void stepfit();
    // lasso pattern search for given lambda series
    void gridSearch(const arma::mat&);
    // Nelder Mead search of the result
    void nelderMead(const arma::colvec& lambda);
    // set order for multivariate Bernoulli
    void setOrder(const unsigned order = 2) {
      ptrLoss -> toSetOrder(order);
      p = ptrLoss -> getDim();
      numRow = p / ptrLoss -> getNumCol();
      constants = arma::zeros<arma::umat>(p, 1);
      for (unsigned i = 0; i < p; i++)
	constants(i) = i;
    }
  };
}

#endif // LPS_H_
