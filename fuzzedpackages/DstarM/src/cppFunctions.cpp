// for speed!
#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// include VOSS & VOSS code
extern "C" {
	#include "rtdists.h"
}

using namespace arma;

/* unused alternative for trapz */
// [[Rcpp::export]]
arma::vec simpsonC(const arma::vec& x, const arma::mat& fx) {

	// 1D integration.
	// input:
	// x: vector with grid of stuff to integrate
	// fx: matrix containing what to integrate. Every column is seen as a new thing to integrate.

	// output:
	// vector with integral of every column.

	const int n = x.n_elem;  // size of grid
	const int p = fx.n_cols; // number of integrals

	int n2 = 2*n - 1; // size of interpolated grid
	arma::vec x2 = arma::linspace<arma::vec>(x(0), x(n-1), n2); // interpolated grid
	arma::mat fx2(n2, p);
	arma::vec tmp(n2);

	for (int i = 0; i < p; ++i) {
		arma::interp1(x, fx.col(i), x2, tmp, "*linear");
		fx2.col(i) = tmp;
	}

	double h = (x2(2) - x2(1)) / 3.0;

	arma::vec out(p, arma::fill::zeros);

	for (int i = 0; i < (n - 1); ++i) {

		out += trans(fx2.row(2*i) + 4.0 * fx2.row(2*i + 1) + fx2.row(2*i + 2));

	}

	return out * h;

}

// [[Rcpp::export]]
arma::vec dunifc(const arma::vec& x, const double& a, const double& b) {

	// OUTPUT:
	// uniform density function
	// INPUT:
	// x: points to evaluate the density on.
	// a: lower bound of uniform density.
	// b: upper bound of uniform density.

	arma::vec out(x.n_elem);
	double d = 1.0 / (b-a);
	for (arma::uword i = 0; i < x.n_elem; ++i) {

		if (x(i) >= a && x(i) <= b) {
			out(i) = d;
		} else {
			out(i) = 0.0;
		}

	}

	return out;

}

// [[Rcpp::export]]
arma::vec convolveC(const arma::vec& x, const arma::vec& y) {

	// x, y: equal length.
	// output has length x.

	arma::vec z = arma::conv(x, y);
	return z.subvec(0, size(x));
	//return arma::conv(x, y, "same");

}

// [[Rcpp::export]]
arma::mat convolveC2(arma::mat& x, arma::mat& y) {

	const int nr = x.n_rows;
	const int nc = x.n_cols;

	arma::mat out(2*nr - 1, nc);

	for (int i = 0; i < nc; ++i) {

		out.col(i) = arma::conv(x.col(i), y.col(i));

	}

	return out;

}

// [[Rcpp::export]]
double chisqC(const arma::vec& tt, const arma::vec& a, const arma::vec& b) {

	arma::vec vals = pow(a - b, 2) / (a + b + 1e-10);
	return arma::as_scalar(arma::trapz(tt, vals));

}

// [[Rcpp::export]]
double rObjC3(arma::vec& r, arma::vec& tt, arma::vec& a, arma::vec& bb, arma::vec& lenPre, arma::vec& lenPost) {

	// both lenPre and lenPost
	arma::vec bb0 = join_cols(lenPre, r);
	arma::vec bb1 = join_cols(bb0, lenPost);
	arma::vec bb2 = arma::conv(bb1, bb);
	arma::vec bb3 = bb2.rows(0, a.size()-1);
	return chisqC(tt, a, bb3);

}

// [[Rcpp::export]]
double rObjC2(arma::vec& r, arma::vec& tt, arma::vec& a, arma::vec& bb, arma::vec& lenPre) {

	// no lenPost
	arma::vec bb1 = join_cols(lenPre, r);
	arma::vec bb2 = arma::conv(bb1, bb);
	arma::vec bb3 = bb2.rows(0, a.size()-1);
	return chisqC(tt, a, bb3);

}

// [[Rcpp::export]]
double rObjC1(arma::vec& r, arma::vec& tt, arma::vec& a, arma::vec& bb, arma::vec& lenPost) {

	// no lenPre
	arma::vec bb1 = join_cols(r, lenPost);
	arma::vec bb2 = arma::conv(bb1, bb);
	arma::vec bb3 = bb2.rows(0, a.size()-1);
	return chisqC(tt, a, bb3);

}

// [[Rcpp::export]]
double rObjC0(arma::vec& r, arma::vec& tt, arma::vec& a, arma::vec& bb) {

	// no lenPre or lenPost
	arma::vec bb2 = arma::conv(r, bb);
	arma::vec bb3 = bb2.rows(0, a.size()-1);
	return chisqC(tt, a, bb3);

}


// [[Rcpp::export]]
double nthMomentSC(const arma::vec& x, const arma::vec& fx, const int& nth) {

	return arma::as_scalar(arma::trapz(x, pow(x, nth) % fx));
}

// [[Rcpp::export]]
double nthCMomentSC(const arma::vec& x, const arma::vec& fx, const int& nth) {

	double ex = arma::as_scalar(arma::trapz(x, x % fx));
	arma::vec dif = x - ex;
	return arma::as_scalar(arma::trapz(x, pow(dif, nth) % fx));

}

// [[Rcpp::export]]
arma::vec getVarC(arma::mat Pdf, const arma::vec& tt, const arma::mat& mm2) {

	// sum degenerate pdfs
	Pdf = Pdf * mm2;
	Pdf = Pdf * arma::diagmat(1.0 / arma::trapz(tt, Pdf));

	const int nc = Pdf.n_cols;
	arma::vec out(nc);
	for (int i = 0; i < nc; ++i) {

		out(i) = nthCMomentSC(tt, Pdf.col(i), 2);

	}

	return out;

}

// [[Rcpp::export]]
bool oscCheckC(const arma::mat& x) {

	// x: matrix where every column is a pdf.
	// checks if any pdf is multimodal, stops at first encounter (false for no error, true for error).

	const int nr = x.n_rows;
	const int nc = x.n_cols;
	int i, j;

	for (int c = 0; c < nc; c++) {

		// (re)set i
		i = 1;

		// loop and check if still ascending the density
		while ((i < nr) && (x(i-1, c) <= x(i, c))) ++i;

		// first mode found at i
		j = i;

		// loop and check if still descending the density
		while ((j < nr) && (x(j-1, c) >= x(j, c))) j++;

		// if true implies multimodality, therefore exit
		if (j != nr) {
			return true;
		}
	}

	// if execution reaches here then only one mode was found for all the densities.
	return false;

}

// only function that requires external C code
// [[Rcpp::export]]
arma::mat getVoss(arma::vec& rt, arma::mat& pars, const double& precision) {

	// rescale parameters
	// pars.row(7) /= pars.row(0); // z

	// arma::vec a(pars.n_cols);
	// arma::vec b(pars.n_cols);
	// double by = rt(1) - rt(0);

	for (unsigned int i = 0; i < pars.n_cols; ++i) {
		// sz
		if (pars(7, i) < 0.5) {
			pars(4, i) *= 2.0 * pars(7, i);
		} else {
			pars(4, i) *= 2.0 * (1.0 - pars(7, i));
		}

		// a(i) = pars(2, i) - pars(6, i) * pars(2, i);
		// b(i) = pars(2, i) + pars(6, i) * pars(2, i);

		// pars(6, i) = 0;
		pars(2, i) += pars(6, i) / 2.0;
		// pars(2, i) = 0;
		// st0//
		// pars(6, i) = 2.0 * pars(2, i);

	}

	// initialize
	int sum = rt.n_elem;
	arma::vec dens_upp(sum, arma::fill::zeros);
	arma::vec dens_low(sum, arma::fill::zeros);
	arma::mat dens(sum, 2*pars.n_cols);

	// Rcpp::Rcout << "dens:\n" << dens.n_cols << std::endl;

	// make pointers
	int *in_numvalues = &sum;

	double *in_RTs = rt.memptr();
	double inPrec = precision;
	double *in_precision = &inPrec;
	double *out_densities_u = dens_upp.memptr();
	double *out_densities_l = dens_low.memptr();
	double *in_params;

	int j = 0;
	for (unsigned int i = 0; i < pars.n_cols; ++i) {

		// pointer to parameters.
		in_params = pars.colptr(i);

		// call VOSS code.
		dfastdm(in_numvalues, in_params, in_RTs, in_precision, out_densities_u, out_densities_l);

		// if (b(i) == 0.0 || by > (b(i) - a(i))) { // conditions for absent nondecision pdf

			//Rcpp::Rcout << "no Conv\n" << std::endl;
  	dens.col(j) = abs(dens_low);
  	dens.col(j+1) = dens_upp;

		// } //else { // manual convolution with uniform nd

			//Rcpp::Rcout << "manual Conv\n" << std::endl;

		// 	arma::vec nondecisionPdf = dunifc(rt, a(i), b(i));
		// 	dens.col(j) = convolveC(abs(dens_low), nondecisionPdf);
		// 	dens.col(j+1) = convolveC(dens_upp, nondecisionPdf);
		//
		// }

		// increment j
		j += 2;

	}

	return dens;
}

// [[Rcpp::export]]
void imposeFixationsC(arma::vec& pars, const arma::mat fixed) {

	for (unsigned int i = 0; i < fixed.n_cols; ++i) {

		pars.insert_rows(fixed(1, i), 1);

		if (fixed(0, i) == 1) { // TRUE: fixed value, otherwise function

			pars(fixed(1, i)) = fixed(2, i);

		} else if (fixed(4, i) == 0) { // addition

			pars(fixed(1, i)) = fixed(2, i) + pars(fixed(3, i));

		} else if (fixed(4, i) == 1) { // substraction

			pars(fixed(1, i)) = fixed(2, i) - pars(fixed(3, i));

		} else if (fixed(4, i) == 2) { // multiplication

			pars(fixed(1, i)) = fixed(2, i) * pars(fixed(3, i));

		} else { // division

			pars(fixed(1, i)) = fixed(2, i) / pars(fixed(3, i));

		}

		// Rcpp::Rcout << "pars = " << pars << std::endl;

	}

	// return pars;

}

// [[Rcpp::export]]
arma::mat getPdfC(arma::vec& tt, arma::mat pars, const arma::mat& mm, const bool& DstarM, const bool& oscPdf, const double& precision) {

	// Pdfs, time grid, condition matrix, yes/ no DstarM analysis, check for oscillations.
	// note that the 'A' matrix is used to report errors (instead of previous NULL).

	// add rows with zeros to parameters matrix
	arma::mat zeroes = arma::zeros<arma::mat>(1, pars.n_cols);

	if (DstarM) {

		// include d, t0, st0
		pars.insert_rows(2, zeroes);
		pars.insert_rows(2, zeroes);
		pars.insert_rows(7, zeroes);
		pars.insert_rows(7, zeroes); // additional column for zr

	} else {

		// include d
		pars.insert_rows(3, zeroes);
		pars.insert_rows(8, zeroes); // additional column for zr

	}

	// switch zr around (see rtdists.h for description)
	pars.row(8) = pars.row(4);
  pars.shed_row(4);

	// call Voss & Voss pdfs
	arma::mat pdf = getVoss(tt, pars, precision);

	if (oscPdf && oscCheckC(pdf)) {
		arma::mat A(1, 1, arma::fill::ones);
		return A;
	}

	arma::mat cor = 1.0 / arma::trapz(tt, pdf * mm);
	if (cor.has_inf()) {
		arma::mat A(1, 1, arma::fill::ones);
		return A;
	}

	arma::vec cor2(2*cor.n_elem);
	for (unsigned int i = 0; i < cor.n_elem; i ++) {
		cor2(2*i) = cor(i);
		cor2(2*i + 1) = cor(i);
	}

	// Rcpp::Rcout << "cor2: " << cor2 << std::endl;

	// return renormalized pdfs
	return pdf * arma::diagmat(cor2);

}

// [[Rcpp::export]]
double totalobjectiveC(arma::vec pars, arma::vec& tt, const arma::vec& ql, const arma::vec& ii, const arma::vec& jj, const arma::vec& varData,
                       const arma::mat& g, arma::mat restr, const arma::mat& mm, const arma::mat& mm2,
                       const bool& DstarM, const bool& oscPdf, const bool& forceRestriction, double precision,
                       const bool& anyFixed, arma::mat fixed) {

	double out = 0.0;

	if (anyFixed) {
		imposeFixationsC(pars, fixed);
	}

	// convert restr to matrix of parameters
	for (unsigned int i = 0; i < restr.n_rows; ++i) {
		for (unsigned int j = 0; j < restr.n_cols; j++) {
			restr(i, j) = pars(restr(i, j));
		}
	}

	arma::mat pdf = getPdfC(tt, restr, mm, DstarM, oscPdf, precision);

	if (pdf.n_elem == 1) {
		return 1e9;
	}

	if (!DstarM) {

		for (unsigned int i = 0; i < pdf.n_cols; ++i) {

			out += chisqC(tt, pdf.col(i), g.col(i)) * 100 * ql(i) / sum(ql);

		}

	} else {

		if (forceRestriction) {

			arma::vec varModel = getVarC(pdf, tt, mm2);
			arma::vec varNonDec = varData - varModel;
			if (arma::any(varNonDec < 0)) { // % (arma::abs(varNonDec) < pow(2, -52)))) { # redundant?

				// Rcpp::Rcout << "var restriction " << 1 << std::endl;
				return 1e9;
			}

		}

		arma::vec a;
		arma::vec b;
		for (unsigned int l = 0; l < ii.n_elem; l++) {
			a = convolveC(g.col(ii(l)), pdf.col(jj(l)));
			b = convolveC(g.col(jj(l)), pdf.col(ii(l)));

			out += chisqC(tt, a, b) * 100 * (ql(ii(l)) + ql(jj(l))) / sum(ql);

		}

	}

	return out;

}
