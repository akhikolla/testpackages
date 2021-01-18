// Copyright (C) 2012 - 2015  Philip Rinn
// Copyright (C) 2012 - 2015  Carl von Ossietzy Universität Oldenburg
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, see <https://www.gnu.org/licenses/gpl-2.0>.

#include <RcppArmadillo.h>

inline arma::vec linreg(arma::vec x, arma::vec y, arma::vec w) {
	// y = b + a*x
    // w = vector of weights for x
	// returns a vector (b, a)

	// linear regression with weights
	// X = [1; x]   <- matrix
	// W            <- matrix of weights
	// which is in matrix form:
	// t(X) * W * X * coef = t(X) * W * y
	int nsteps = x.n_elem;
	arma::mat X = arma::ones<arma::mat>(nsteps, 2);
	X.col(1) = x;
	arma::mat W = arma::zeros<arma::mat>(nsteps, nsteps);
	// set weights
	W.diag() = w;
	arma::mat tmp = trans(X)*W;
	// catch situations where no solution is found and return R_NaN in this case
    arma::vec res(2);
    try {
		res = solve(tmp*X, tmp*y);
	}
	catch (const std::runtime_error& error) {
		res.fill(R_NaN);
	}
	return res;
}

// to make it easy to handle no given weights
inline arma::vec linreg(arma::vec x, arma::vec y) {
	return linreg(x, y, arma::ones<arma::vec>(x.n_elem));
}

