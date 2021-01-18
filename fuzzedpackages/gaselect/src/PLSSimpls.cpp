//
//  PLSSimpls.cpp
//  gaselect
//
//  Created by David Kepplinger on 16.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "config.h"

#include <iostream>
#include <stdexcept>
#include <RcppArmadillo.h>
#include "PLSSimpls.h"

const double PLSSimpls::NORM_TOL = 1e-25;

PLSSimpls::PLSSimpls(const arma::mat &X, const arma::vec &Y) : PLS(X, Y) {
}

PLSSimpls::~PLSSimpls() {
}

std::unique_ptr<PLS> PLSSimpls::clone() const {
	return std::unique_ptr<PLS>(new PLSSimpls(arma::mat(this->X), arma::mat(this->Y)));
}

inline void PLSSimpls::centerView() {
	// center X and Y

	switch(this->currentViewState) {
		case PLS::UNKNOWN:
			this->viewX = this->X;
			this->viewY = this->Y;
			break;
		case PLS::COLUMNS:
			this->viewX = this->viewXCol;
			this->viewY = this->Y;
			break;
		default:
			break;
	}

	this->Xmean = arma::mean(this->viewX);
	this->Ymean = arma::mean(this->viewY);

	this->viewX.each_row() -= this->Xmean;
	this->viewY -= this->Ymean;
}

/* Small highly optimized functions */
namespace {
/**
 * Stabilized Gram-Schmidt orthogonalization WITHOUT norming the vector! This is done in
 * the call to `manualDeflate` more efficiently.
 */
inline void stabilizedGramSchmidt(arma::mat &V, arma::uword col) {
	arma::vec v = V.unsafe_col(col);

	for (arma::uword i = 0; i < col; ++i) {
		v -= V.col(i) * arma::dot(V.col(i), v);
	}
}

/**
 * Deflate the vector s according to the formula
 * v = v / norm(v)
 * s = s - v * v.t() * s
 */
inline void manualDeflate(arma::vec& s, arma::vec &v) {
	double *vm = v.memptr();
	double *sm = s.memptr();
	arma::uword i;
	double dot = 0, norm = 0;

	for (i = 0; i < s.n_elem; ++i) {
		norm += vm[i] * vm[i];
		dot += vm[i] * sm[i];
	}

	dot /= norm;
	norm = sqrt(norm);

	for (i = 0; i < s.n_elem; ++i) {
		sm[i] -= vm[i] * dot;
		vm[i] /= norm;
	}
}

/**
 * Fast version of
 * this->coef.col(i) = this->coef.col(i - 1) + S * q / tnorm;
 * without the need of any memory copying whatsoever.
 */
inline void updateCoefs(double *Cm, const arma::vec& s, const double a, const arma::uword col) {
	const double *sm = s.memptr();
	for (arma::uword i = 0, n = col * s.n_elem, o = (col - 1) * s.n_elem; i < s.n_elem; ++o, ++n, ++i) {
		Cm[n] = Cm[o] + sm[i] * a;
	}
}
}

/**
 * ncomp is 1-based!!!
 */
void PLSSimpls::fit(uint16_t ncomp) {
	uint16_t maxNComp = ((this->viewX.n_cols < this->viewX.n_rows) ? this->viewX.n_cols : this->viewX.n_rows - 1);
	if(ncomp == 0 || ncomp > maxNComp) {
		ncomp = maxNComp;
	}

	/*
	 * Center X and Y views
	 */
	this->centerView();

	/*
	 * Init neccessary matrices and vectors
	 * Variable names are according to the original paper by S. de Jong (1993)
	 */

	this->coef.zeros(this->viewX.n_cols, ncomp);
	this->intercepts.zeros(ncomp);
	this->V.set_size(this->viewX.n_cols, ncomp);

	arma::vec S = this->viewX.t() * this->viewY; // Cross product

	// Working vectors
	arma::vec t; // X block factor scores
	double tnorm = 1.0;
	double q;

	for(uint16_t i = 0; i < ncomp; ++i) {
		arma::vec v = this->V.unsafe_col(i);
		t = this->viewX * S;

		t -= arma::mean(t); // Center y block factor scores
		tnorm = arma::norm(t, 2); // Calculate norm

		// The norm of t can be zero (or close to it). This is unacceptable.
		if (tnorm < PLSSimpls::NORM_TOL) {
			throw std::underflow_error("All block-factor scores are (almost) zero.");
		}

		t /= tnorm;  // Normalize scores

		v = this->viewX.t() * t; // Calculate x loadings
		q = arma::dot(this->viewY, t); // Calculate y loadings

		if(i > 0) {
			stabilizedGramSchmidt(this->V, i); // Make v orthogonal to previous loadings
			updateCoefs(this->coef.memptr(), S, q / tnorm, i);
//			this->coef.col(i) = this->coef.col(i - 1) + S * q / tnorm;
		} else {
			this->coef.col(i) = S * q / tnorm;
		}


		/* deflate S and norm v */
		manualDeflate(S, v);

		this->intercepts[i] = this->Ymean - arma::dot(this->Xmean, this->coef.col(i));
	}

	this->resultNComp = ncomp;
}
