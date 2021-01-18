/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include <R.h>
#include "cor.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// C++ isnan() is not portable and gives error on Windows systems
// use R macro ISNAN() instead

// -------------------
// Pearson correlation
// -------------------

// barebones arma version of the Pearson correlation
//double corPearson(const vec& x, const vec& y) {
//	const uword n = x.n_elem;
//	// compute means
//	double meanX = 0, meanY = 0;
//	for(uword i = 0; i < n; i++) {
//		meanX += x(i);
//		meanY += y(i);
//	}
//	meanX /= n;
//	meanY /= n;
//	// compute Pearson correlation
//	// covariance and variances are computed up to scaling factor (n-1)
//	double covXY = 0, varX = 0, varY = 0;
//	for(uword i = 0; i < n; i++) {
//		double tmpX = x(i) - meanX, tmpY = y(i) - meanY;  // centered values
//		covXY += tmpX * tmpY;
//		varX += tmpX * tmpX;
//		varY += tmpY * tmpY;
//	}
//	return covXY / (sqrt(varX) * sqrt(varY));
//}
double corPearson(const vec& x, const vec& y) {
	mat corXY = cor(x, y);	// arma function cor() always returns matrix
	return corXY(0, 0);
}

// R interface
SEXP R_corPearson(SEXP R_x, SEXP R_y) {
	// convert data
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	vec y(Rcpp_y.begin(), Rcpp_y.size(), false);	// reuse memory
	// call arma version and wrap result
	return wrap(corPearson(x, y));
}


// --------------------
// Spearman correlation
// --------------------

// barebones arma version of the Spearman correlation
double corSpearman(const vec& x, const vec& y) {
	const uword n = x.n_elem;
	// make sure it returns NA in the presence of NA's
	for(uword i = 0; i < n; i++) {
		if(ISNAN(x(i)) || ISNAN(y(i))) return(NA_REAL);
	}
	// compute ranks
	vec ranksX = rank_ccaPP(x), ranksY = rank_ccaPP(y);
	// return Pearson correlation for ranks
	return corPearson(ranksX, ranksY);
}

// version with parameter that determines whether consistent estimate at the
// normal model should be computed
double corSpearman(const vec& x, const vec& y, const bool& consistent) {
	double r = corSpearman(x, y);	// Spearman correlation
	if(consistent) {
		r = 2 * sin(M_PI * r / 6);	// consistent estimate at the normal model
	}
	return r;
}

// R interface
SEXP R_corSpearman(SEXP R_x, SEXP R_y, SEXP R_consistent) {
	// convert data
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	vec y(Rcpp_y.begin(), Rcpp_y.size(), false);	// reuse memory
	bool consistent = as<bool>(R_consistent);
	// call arma version and wrap result
	return wrap(corSpearman(x, y, consistent));
}


// -------------------
// Kendall correlation
// -------------------

// naive n^2 implementation
double naiveCorKendall(const vec& x, const vec& y, const uword& n) {
	// compute sum of signs for pairs of observations
	double sum = 0, xi, xj, yi, yj;
	sword signX, signY;				// signed integers
	uword nTieX = 0, nTieY = 0;		// unsigned integers
	for(uword j = 0; j < n; j++) {
		xj = x(j); yj = y(j);
		for(uword i = 0; i < j; i++) {
			// obtain sign and update number of ties for pair of x values
			xi = x(i);
			if(xi > xj) {
				signX = 1;
			} else if(xi < xj) {
				signX = -1;
			} else {
				signX = 0;
				nTieX++;
			}
			// obtain sign and update number of ties for pair of y values
			yi = y(i);
			if(yi > yj) {
				signY = 1;
			} else if(yi < yj) {
				signY = -1;
			} else {
				signY = 0;
				nTieY++;
			}
			// update sum for numerator
			sum += signX * signY;
		}
	}
	// ties need to be subtracted from number of pairs for denominator
	uword nPairs = n * (n-1) / 2;
	double dX = nPairs - nTieX, dY = nPairs - nTieY;
	// return Kendall correlation
	return sum / (sqrt(dX) * sqrt(dY));
}

// barebones arma version of the Kendall correlation
double corKendall(const vec& x, const vec& y) {
	const uword n = x.n_elem;
	// make sure it returns NA in the presence of NA's
	for(uword i = 0; i < n; i++) {
		if(ISNAN(x(i)) || ISNAN(y(i))) return(NA_REAL);
	}
	// call naive version for small n and fast version for larger n
	double r;
	if(n < 30) {
		r = naiveCorKendall(x, y, n);
	} else {
		r = fastCorKendall(x, y, n);
	}
	return r;
}

// version with parameter that determines whether consistent estimate at the
// normal model should be computed
double corKendall(const vec& x, const vec& y, const bool& consistent) {
	double r = corKendall(x, y);	// Kendall correlation
	if(consistent) {
		r = sin(M_PI * r / 2);		// consistent estimate at the normal model
	}
	return r;
}

// R interface
SEXP R_corKendall(SEXP R_x, SEXP R_y, SEXP R_consistent) {
	// convert data
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	vec y(Rcpp_y.begin(), Rcpp_y.size(), false);	// reuse memory
	bool consistent = as<bool>(R_consistent);
	// call arma version and wrap result
	return wrap(corKendall(x, y, consistent));
}


// --------------------
// quadrant correlation
// --------------------

// barebones arma version of the quadrant correlation
double corQuadrant(const vec& x, const vec& y) {
	const uword n = x.n_elem;
	// make sure it returns NA in the presence of NA's
	for(uword i = 0; i < n; i++) {
		if(ISNAN(x(i)) || ISNAN(y(i))) return(NA_REAL);
	}
	// count number of observations in 1./3. and 2./4. quadrants, respectively
	double medX = median(x), medY = median(y);
	sword onethree = 0, twofour = 0;
	for(uword i = 0; i < n; i++) {
		// this should be numerically stable
		double xi = x(i), yi = y(i);
	    if(((xi > medX) && (yi > medY)) || ((xi < medX) && (yi < medY))) {
	    	onethree++;
	    } else if(((xi > medX) && (yi < medY)) || ((xi < medX) && (yi > medY))) {
	    	twofour++;
	    }
	}
	// return quadrant correlation
	return ((double) (onethree - twofour)) / ((double) (onethree + twofour));
}

// version with parameter that determines whether consistent estimate at the
// normal model should be computed
double corQuadrant(const vec& x, const vec& y, const bool& consistent) {
	double r = corQuadrant(x, y);	// quadrant correlation
	if(consistent) {
		r = sin(M_PI * r / 2);		// consistent estimate at the normal model
	}
	return r;
}

// R interface
SEXP R_corQuadrant(SEXP R_x, SEXP R_y, SEXP R_consistent) {
	// convert data
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	vec y(Rcpp_y.begin(), Rcpp_y.size(), false);	// reuse memory
	bool consistent = as<bool>(R_consistent);
	// call arma version and wrap result
	return wrap(corQuadrant(x, y, consistent));
}


// ----------------------
// Huber-type M-estimator
// ----------------------

// barebones arma version of the Huber-type M-estimator
double corM(const vec& x, const vec& y, const double& prob,
		const string& initial, const double& tol) {
	const uword n = x.n_elem;
	// note that (1 - prob) needs to be used in the call to Rf_qchisq()
	double d = Rf_qchisq(1 - prob, 2, 0, 0);
	// compute initial center, covariance matrix and correlation
	double centerX, centerY, r;
	mat xy(n, 2), covMat(2, 2);
	// compute initial correlation
	if(initial == "pearson") {
		r = corPearson(x, y);
	} else if(initial == "quadrant") {
		r = corQuadrant(x, y, true);
	} else if(initial == "spearman") {
		r = corSpearman(x, y, true);
	} else if(initial == "kendall") {
		r = corKendall(x, y, true);
	} else if(initial == "spearman") {
		r = corSpearman(x, y, true);
	} else {
		error("method not available");	// should never happen
	}
	if((1 - abs(r)) > tol) {
		// covariance matrix would be singular otherwise
		// build initial covariance matrix
		if(initial == "pearson") {
			centerX = mean(x); centerY = mean(y);	// means of x and y
			covMat(0, 0) = var(x);
			covMat(1, 1) = var(y);

		} else {
			// check whether the MAD of x is 0 and fall back on the
			// standard deviation in that case
			double tmp = mad(x, centerX);
			if(tmp == 0) {
				centerX = mean(x);
				covMat(0, 0) = var(x);
			} else {
				covMat(0, 0) = pow(tmp, 2);
			}
			// check whether the MAD of y is 0
			tmp = mad(y, centerY);
			if(tmp == 0) {
				centerY = mean(y);
				covMat(1, 1) = var(y);
			} else {
				covMat(1, 1) = pow(tmp, 2);
			}
		}
		covMat(1, 0) = r * sqrt(covMat(0, 0) * covMat(1, 1));
		covMat(0, 1) = covMat(1, 0);
		xy = join_rows(x - centerX, y - centerY);	// centered x and y
		// iteratively compute M-estimator
		double previousFisherZ = R_NegInf, fisherZ = atanh(r);	// Fisher transform
		while((abs(fisherZ - previousFisherZ) > tol) && ((1 - abs(r)) > tol)) {
			// second condition ensures that covariance matrix is not singular
			// compute Mahalanobis distances
			mat invCovMat = inv(covMat);			// inverse covariance matrix
			vec md = sum((xy * invCovMat) % xy, 1);	// squared Mahalanobis distances
			// compute weights based on Mahalanobis distances
			vec w(n);
			for(uword i = 0; i < n; i++) {
				if(md(i) > d) {
					w(i) = sqrt(d / md(i));
				} else {
					w(i) = 1;
				}
			}
			// update location estimates
			double sumW = sum(w);
			centerX = dot(w, x) / sumW;
			centerY = dot(w, y) / sumW;
			xy = join_rows(x - centerX, y - centerY);	// centered x and y
			// update scatter matrix
			mat wxy = join_rows(w % xy.unsafe_col(0), w % xy.unsafe_col(1));
			covMat = trans(wxy) * wxy / (n * prob);
			// update correlation
			r = covMat(1, 0) / sqrt(covMat(0, 0) * covMat(1, 1));
			// compute Fisher transformation
			previousFisherZ = fisherZ;
			fisherZ = atanh(r);
		}
	}
	return r;
}

// R interface to corM()
SEXP R_corM(SEXP R_x, SEXP R_y, SEXP R_prob, SEXP R_initial, SEXP R_tol) {
	// convert data to arma types
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	vec y(Rcpp_y.begin(), Rcpp_y.size(), false);	// reuse memory
	double prob = as<double>(R_prob);
	string initial = as<string>(R_initial);			// convert character string
	double tol = as<double>(R_tol);
	// call arma version and wrap result
	return wrap(corM(x, y, prob, initial, tol));
}
