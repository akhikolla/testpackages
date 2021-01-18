#include "lrfunc.h"
#include "aag.h"
using namespace arma;
using namespace Rcpp;

vec lrpieabs(const vec& y, const mat& x, const vec& b, int n, double lam2, const mat& l) {
	return abs(x.t()*(y - x*b) / n - lam2*l*b);
}

vec lraen(vec b, double lam1, double lam2, const vec& w, const mat& x, const vec& y, const vec& dl, int maxiter, double cri) {
	int p = b.n_elem, n = y.n_elem, act0 = p, act1;
	double bt, loss0, loss1, softer0, softer1;
	vec rr, lw = lam1*w, xsq = trans(mean(pow(x, 2))), la1 = xsq + lam2*dl, la12 = lam2*dl;
	uvec actset(p, fill::ones);
	vec::iterator bpoint, xsqpoint, lapoint, lwpoint;
	uvec::iterator actpoint;
	rr = y - x*b;
	loss0 = dot(rr, rr) / n + dot(lw, abs(b)) + dot(b, la12 % b);
	for (int k = 0; k < maxiter; ++k) {
		act1 = 0;
		bpoint = b.begin();
		xsqpoint = xsq.begin();
		actpoint = actset.begin();
		lapoint = la1.begin();
		lwpoint = lw.begin();
		for (int q = 0; q < p; ++q) {
			if (*actpoint) {
				bt = *bpoint;
				softer0 = *xsqpoint*bt + dot(x.col(q), rr) / n;
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (*lapoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				rr -= (*bpoint - bt)*x.col(q);
			}
			++bpoint;
			++actpoint;
			++xsqpoint;
			++lapoint;
			++lwpoint;
		}
		act0 -= act1;
		loss1 = dot(rr, rr) / n + dot(lw, abs(b)) + dot(b, la12 % b);
		if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
			break;
		}
		loss0 = loss1;
		if (k == maxiter - 1) {
			Rcout << "Does not converge!" << endl;
		}
	}
	return b;
}

vec lral(vec b, double lam1, const vec& w, const mat& x, const vec& y, int maxiter, double cri) {
	int p = b.n_elem, n = y.n_elem, act0 = p, act1;
	double bt, loss0, loss1, softer0, softer1;
	vec rr, lw = lam1*w, xsq = trans(mean(pow(x, 2)));
	uvec actset(p, fill::ones);
	vec::iterator bpoint, xsqpoint, lwpoint;
	uvec::iterator actpoint;
	rr = y - x*b;
	loss0 = dot(rr, rr) / n + dot(lw, abs(b));
	for (int k = 0; k < maxiter; ++k) {
		act1 = 0;
		bpoint = b.begin();
		xsqpoint = xsq.begin();
		actpoint = actset.begin();
		lwpoint = lw.begin();
		for (int q = 0; q < p; ++q) {
			if (*actpoint) {
				bt = *bpoint;
				softer0 = *xsqpoint*bt + dot(x.col(q), rr) / n;
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / *xsqpoint;
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				rr -= (*bpoint - bt)*x.col(q);
			}
			++bpoint;
			++actpoint;
			++xsqpoint;
			++lwpoint;
		}
		act0 -= act1;
		loss1 = dot(rr, rr) / n + dot(lw, abs(b));
		if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
			break;
		}
		loss0 = loss1;
		if (k == maxiter - 1) {
			Rcout << "Does not converge!" << endl;
		}
	}
	return b;
}
