#include "logifunc.h"
#include "aag.h"
using namespace arma;
using namespace Rcpp;

double auc(int n, const vec& prob, const vec& yy) {
	if ((all(yy == 1)) | (all(yy == 0))) {
		Rcout << "Does not contain two classes!" << endl;
		return 0;
	}
	if (all(prob == prob(0))) {
		return .5;
	}
	vec rocy, y0, y1, stackx, stacky;
	rocy = yy(sort_index(prob, "descend"));
	y0 = conv_to<vec>::from(rocy == 0);
	y1 = conv_to<vec>::from(rocy == 1);
	stackx = cumsum(y0) / sum(y0);
	stacky = cumsum(y1) / sum(y1);
	return dot(diff(stackx), stacky(span(1, n - 1)));
}

double logi1(double my) {
	return log(my / (1 - my));
}

vec logipieabs(const vec& y, const mat& x, double b0, const vec& b, int n, double lam2, const mat& l) {
	return abs(x.t()*(y - (1 / (1 + exp(-b0 - x*b)))) / n - lam2*l*b);
}

vec logiaen(double b0, vec b, double lam1, double lam2, const vec& w, const mat& x, const vec& y, const vec& dl, bool intercept, int maxiter, double cri) {
	int p = b.n_elem, act0 = p, act1;
	double del, b00, bq, loss0, loss1, softer0, softer1;
	vec lw = lam1*w, ww, pp, bx = b0 + x*b, ldl = lam2*dl, tmpone;
	mat xsq2 = pow(x, 2);
	uvec actset(p, fill::ones), tmp1;
	pp = 1 / (1 + exp(-bx));
	ww = pp % (1 - pp);
	loss0 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b)) + dot(b, ldl % b);
	vec::iterator bpoint, lwpoint, ldlpoint;
	uvec::iterator actpoint;
	for (int k = 0; k < maxiter; ++k) {
		act1 = 0;
		bpoint = b.begin();
		lwpoint = lw.begin();
		ldlpoint = ldl.begin();
		actpoint = actset.begin();
		for (int q = 0; q < p; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				del = mean(ww%xsq2.col(q));
				softer0 = del*bq + mean(x.col(q) % (y - pp));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (del + *ldlpoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx += (*bpoint - bq)*x.col(q);
				pp = 1 / (1 + exp(-bx));
				ww = pp % (1 - pp);
				if (any(bx > 10)) {
					tmp1 = find(bx > 10);
					tmpone = ones<vec>(tmp1.n_elem);
					pp(tmp1) = tmpone;
					ww(tmp1) = .0001*tmpone;
				}
				if (any(bx < -10)) {
					tmp1 = find(bx < -10);
					tmpone = zeros<vec>(tmp1.n_elem);
					pp(tmp1) = tmpone;
					ww(tmp1) = .0001 + tmpone;
				}
			}
			++bpoint;
			++actpoint;
			++lwpoint;
			++ldlpoint;
		}
		if (intercept) {
			b00 = mean(y - pp) / mean(ww);
			b0 += b00;
			bx += b00;
		}
		act0 -= act1;
		loss1 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b)) + dot(b, ldl%b);
		if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
			break;
		}
		if (k == maxiter - 1) {
			Rcout << "Does not converge!" << endl;
		}
		loss0 = loss1;
	}
	b.insert_rows(0, 1);
	b(0) = b0;
	return b;
}

vec logial(double b0, vec b, double lam1, const vec& w, const mat& x, const vec& y, bool intercept, int maxiter, double cri) {
	int p = b.n_elem, act0 = p, act1;
	double del, b00, bq, loss0, loss1, softer0, softer1;
	vec lw = lam1*w, ww, pp, bx = b0 + x*b, tmpone;
	mat xsq2 = pow(x, 2);
	uvec actset = ones<uvec>(p), tmp1;
	pp = 1 / (1 + exp(-bx));
	ww = pp % (1 - pp);
	loss0 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b));
	vec::iterator bpoint, lwpoint;
	uvec::iterator actpoint;
	for (int k = 0; k < maxiter; ++k) {
		act1 = 0;
		bpoint = b.begin();
		lwpoint = lw.begin();
		actpoint = actset.begin();
		for (int q = 0; q < p; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				del = mean(ww%xsq2.col(q));
				softer0 = del*bq + mean(x.col(q) % (y - pp));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / del;
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx += (*bpoint - bq)*x.col(q);
				pp = 1 / (1 + exp(-bx));
				ww = pp % (1 - pp);
				if (any(bx > 10)) {
					tmp1 = find(bx > 10);
					tmpone = ones<vec>(tmp1.n_elem);
					pp(tmp1) = tmpone;
					ww(tmp1) = .0001*tmpone;
				}
				if (any(bx < -10)) {
					tmp1 = find(bx < -10);
					tmpone = zeros<vec>(tmp1.n_elem);
					pp(tmp1) = tmpone;
					ww(tmp1) = .0001 + tmpone;
				}
			}
			++bpoint;
			++lwpoint;
			++actpoint;
		}
		if (intercept) {
			b00 = mean(y - pp) / mean(ww);
			b0 += b00;
			bx += b00;
		}
		act0 -= act1;
		loss1 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b));
		if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
			break;
		}
		if (k == maxiter - 1) {
			Rcout << "Does not converge!" << endl;
		}
		loss0 = loss1;
	}
	b.insert_rows(0, 1);
	b(0) = b0;
	return b;
}
