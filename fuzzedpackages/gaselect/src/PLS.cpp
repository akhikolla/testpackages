//
//  PLS.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 15.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "config.h"
#include "Logger.h"
#include "PLS.h"

void PLS::viewSelectColumns(const arma::uvec &columns) {
	this->viewXCol = this->X.cols(columns);

	this->currentViewState = COLUMNS;
}

void PLS::viewSelectRows(const arma::uvec &rows) {
	this->viewX = this->viewXCol.rows(rows);
	this->viewY = this->Y.rows(rows);
	this->currentViewState = ROWS;
}

void PLS::viewSelectAllRows() {
	this->viewX = this->viewXCol;
	this->viewY = this->Y;
	this->currentViewState = ROWS;
}

/**
 * uint16_t ncomp ... The 0-based
 */
arma::vec PLS::predict(const arma::mat &newX, uint16_t ncomp) const {
	const arma::mat& coefs = this->getCoefficients();
	const arma::vec& intercepts = this->getIntercepts();
	if(ncomp > coefs.n_cols) {
		GAerr << "Trying to predict with " << ncomp << " components when only " << coefs.n_cols << " components are available" << std::endl;
		throw Rcpp::exception("Can not predict values for a model with more components than fit components", __FILE__, __LINE__);
	}

	--ncomp;

	arma::vec pred = newX * coefs.col(ncomp);
	pred += intercepts[ncomp];
	return pred;
}

arma::mat PLS::predict(const arma::mat &newX) const {
	const arma::mat& coefs = this->getCoefficients();
	const arma::vec& intercepts = this->getIntercepts();
	arma::mat pred(newX.n_rows, coefs.n_cols);
	for(uint16_t i = 0; i < coefs.n_cols; ++i) {
		pred.col(i) = newX * coefs.col(i);
		pred.col(i) += intercepts[i];
	}

	return pred;
}

std::unique_ptr<PLS> PLS::getInstance(PLSMethod method, const arma::mat &X, const arma::vec &Y) {
  return std::unique_ptr<PLS>(new PLSSimpls(X, Y));
	// switch(method) {
	// 	case SIMPLS:
			// return std::unique_ptr<PLS>(new PLSSimpls(X, Y));
	// }
}


