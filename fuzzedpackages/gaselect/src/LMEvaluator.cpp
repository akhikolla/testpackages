//
//  LMEvaluator.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 28.08.2013.
//
//

#include "config.h"

#include <limits>
#include <RcppArmadillo.h>

#include "LMEvaluator.h"
#include "Logger.h"

LMEvaluator::LMEvaluator(const arma::mat &X, const arma::colvec &y, const LMEvaluator::Statistic statistic, const VerbosityLevel &verbosity, const bool addIntercept) : Evaluator(verbosity), y(y), statistic(statistic), Xdesign(X) {
	if(addIntercept) {
		arma::colvec intercept(X.n_rows);
		intercept.ones();
		this->Xdesign.insert_cols(0, intercept);
	}

	this->r2denom = arma::accu(arma::square(this->y - arma::mean(this->y)));
}

double LMEvaluator::evaluate(arma::uvec &columnSubset) {
	double ret = 0.0;
	columnSubset += 1;
	columnSubset.insert_rows(0, 1);
	
	arma::mat Xsub(this->Xdesign.cols(columnSubset));
	try {
		arma::colvec coef = arma::solve(Xsub, this->y);
		arma::colvec residuals = this->y - Xsub * coef;
		
		double RSS = arma::accu(arma::square(residuals));
		
		switch(this->statistic) {
			case BIC: {
				ret = -(Xsub.n_rows * log(RSS / Xsub.n_rows) + Xsub.n_cols * log((double) Xsub.n_rows));
				break;
			}
			case AIC: {
				ret = -(Xsub.n_rows * log(RSS / Xsub.n_rows) + Xsub.n_cols * 2.0);
				break;
			}
			case ADJ_R2: {
				double r2 = 1 - (RSS / this->r2denom);
				ret = 1 - (((Xsub.n_rows - 1) * (1 - r2)) / (Xsub.n_rows - Xsub.n_cols - 1));
				break;
			}
			case R2: {
				ret = 1 - (RSS / this->r2denom);
				break;
			}
			default:
				break;
		}
	} catch(const std::runtime_error& re) {
		throw Evaluator::EvaluatorException("The subset could not be evaluated because the system could not be solved (probably the subset is singular)");
	} catch(...) {
		throw Evaluator::EvaluatorException("The subset could not be evaluated due to an unknown error.");
	}
	
	return ret;
}
