#include "Coxfunc.h"
#include "aag.h"
using namespace arma;
using namespace Rcpp;

int tsign(double ti, int ci, double tj, int cj) {
	if (ti > tj) {
		return cj;
	}
	else if (ti < tj) {
		return -ci;
	}
	else {
		return cj - ci;
	}
}

vec invcumsum(int n, vec& x) {
	vec cx(n, fill::zeros);
	double init = 0;
	vec::iterator x_begin = x.begin(), x_end = x.end(), cx_end = cx.end();
	while (x_end != x_begin) {
		init += *(--x_end);
		*(--cx_end) = init;
	}
	return cx;
}

double cidx(int n, vec score, vec& time, vec& del) {
	int up = 0, down = 0, sct;
	vec::iterator t1 = time.begin() + 1, te = time.end(), d1 = del.begin() + 1, s1 = score.begin() + 1, t2, d2, s2;
	while (t1 != te) {
		t2 = time.begin();
		d2 = del.begin();
		s2 = score.begin();
		while (t2 != t1) {
			sct = tsign(*t1, *d1, *t2++, *d2++);
			up += sct*getsign(*s2++ - *s1);
			down += abs(sct);
		}
		++t1;
		++d1;
		++s1;
	}
	return .5*((double)up / down + 1);
}



void tieup(vec& hh, vec& ee, uvec tie1, uvec tie2) {
	uvec::iterator tie1b = tie1.begin(), tie2b = tie2.begin(), tie1e = tie1.end();
	double htmp = ee(*tie1b++);
	while (tie1b < tie1e) {
		if (*tie2b == *(tie2b + 1)) {
			hh(*tie1b) += htmp;
			htmp += ee(*tie1b++);
			++tie2b;
		}
		else {
			htmp = ee(*tie1b++);
			++tie2b;
		}
	}
}

void tiedown(vec& hh, vec& ee, uvec& tie1, uvec& tie2) {
	uvec::iterator tie1b = tie1.begin(), tie1e = tie1.end() - 1, tie2e = tie2.end() - 1;
	double htmp = ee(*tie1e);
	while (tie1b < tie1e) {
		if (*tie2e == *(tie2e - 1)) {
			hh(*tie1e) += htmp;
			htmp += ee(*tie1e--);
			--tie2e;
		}
		else {
			htmp = ee(*tie1e--);
			--tie2e;
		}
	}
}

vec getdd(vec del, uvec& tie1, uvec& tie2, uvec& tie3) {
	uvec::iterator tie1b = tie1.begin(), tie2b = tie2.begin(), tie3b = tie3.begin(), tie1e = tie1.end();
	int ddtmp = sum(del(span(*tie2b++, *tie3b++)));
	del(*tie1b++) = ddtmp;
	while (tie1b != tie1e) {
		if (*(tie2b - 1) == *tie2b) {
			del(*tie1b++) = ddtmp;
			++tie2b;
			++tie3b;
		}
		else {
			ddtmp = sum(del(span(*(tie2b++), *(tie3b++))));
			del(*tie1b++) = ddtmp;
		}
	}
	return del;
}


vec Coxpieabs(bool tie, int ntie, uvec& tie1, uvec& tie2, vec& del, vec& dd, mat& x, vec& b, int n, double lam2, mat& l) {
	vec eeta = exp(x*b), h00, h0, h;
	h00 = invcumsum(n, eeta);
	if (tie) {
		tieup(h00, eeta, tie1, tie2);
	}
	h0 = (1 / h00) % dd;
	h = cumsum(h0);
	if (tie) {
		tiedown(h, h0, tie1, tie2);
	}
	return abs(x.t()*(del - eeta%h) - n*lam2*l*b);
}

double devianceCox(int n, vec score, bool tie, int ntie, uvec& tie1, uvec& tie2, vec& del) {
	double etmp, mi, mi1;
	int tiei = 0;
	vec escore = exp(score), yscr = invcumsum(n, escore);
	uvec::iterator tie1b = tie1.begin(), tie2b = tie2.begin(), tie1e = tie1.end();
	if (tie) {
		etmp = escore(*tie1b++);
		while (tie1b != tie1e) {
			if (*(tie2b + 1) == *tie2b) {
				yscr(*tie1b) += escore(*tie1b);
				etmp += escore(*tie1b++);
			}
			else {
				etmp += escore(*tie1b++);
			}
			++tie2b;
		}
	}
	yscr = log(yscr);
	if (!is_finite(yscr)) {
		uvec nonfiidx = find_nonfinite(yscr);
		bool pos = any(yscr(nonfiidx) > 0), neg = any(yscr(nonfiidx) < 0), tie0;
		unsigned np;
		vec subscore, subscore2;
		if (pos) {
			np = sum(yscr(nonfiidx) > 0);
			if (tie) {
				tie0 = (tie1(0) < np);
				if (tie0) {
					tiei = (tie1(0) == 0) ? 1 : 0;
				}
			}
			else {
				tie0 = tie;
			}
			mi = max(score);
			subscore = score - mi;
			yscr(0) = mi + log(sum(exp(subscore)));
			if (np > 1) {
				for (unsigned i = 1; i < np; ++i) {
					if (tie0 && any(tie1 == i) && (tie1(tiei) != tie2(tiei))) {
						yscr(i) = yscr(i - 1);
						++tiei;
					}
					else {
						if (tie0 && any(tie1 == i) && (tie1(tiei) == tie2(tiei))) {
							++tiei;
						}
						mi1 = max(score.tail(n - i));
						subscore += mi - mi1;
						mi = mi1;
						yscr(i) = mi + log(sum(exp(subscore.tail(n - i))));
					}
				}
			}
		}
		if (neg) {
			np = sum(yscr(nonfiidx) < 0);
			if (tie) {
				tie0 = (tie1(ntie - 1) >= n - np);
				if (tie0) {
					uvec tietmp = find(tie1 >= n - np);
					tiei = tie1(tietmp(0)) == n - np ? tietmp(1) : tietmp(0);
				}
			}
			else {
				tie0 = tie;
			}
			mi = min(score.tail(np));
			subscore = score.tail(np) - mi;
			yscr(n - np) = mi + log(sum(exp(subscore)));
			if (np > 1) {
				for (unsigned i = n - np + 1; i < (unsigned)n; i++) {
					if (tie0&&any(tie1 == i) && (tie1(tiei) != tie2(tiei))) {
						yscr(i) = yscr(i - 1);
						++tiei;
					}
					else {
						if (tie0&&any(tie1 == i) && (tie1(tiei) == tie2(tiei))) {
							++tiei;
						}
						mi1 = min(score.tail(n - i));
						subscore += mi - mi1;
						mi = mi1;
						yscr(i) = mi + log(sum(exp(subscore.tail(n - i))));
					}
				}
			}
		}
	}
	return dot(del, score - yscr);
}

vec Coxnet0(int n, int p, int ntie, vec b, mat& x, bool tie, uvec& tie1, uvec& tie2, vec& del, vec& dd, mat l, int maxiter, double cri) {
	double loss0, loss1;
	vec dl = l.diag(), zeta(n, fill::zeros), eh, h00, h0, h, eta, eeta, delta;
	mat nl = n*l, deye = eye<mat>(n, n), ddd, z, xdx;
	bool ridge = all(findzerocol(l, p));
	if (p > n) {
		z = x*solve(nl, x.t());
		eta = z*zeta;
		eeta = exp(eta);
	}
	else {
		eta = x*b;
		eeta = exp(eta);
	}
	h00 = invcumsum(n, eeta);
	if (tie) {
		tieup(h00, eeta, tie1, tie2);
	}
	h0 = (1 / h00) % dd;
	h = cumsum(h0);
	if (tie) {
		tiedown(h, h0, tie1, tie2);
	}
	eh = eeta%h;
	delta = del - eh;
	ddd = diagmat(eh);
	if (p > n) {
		loss0 = dot(del, log(h00) - eta) / n + dot(zeta, z*zeta);
	}
	else {
		if (ridge) {
			loss0 = dot(del, log(h00) - eta) / n + dot(b, dl%b);
		}
		else {
			loss0 = dot(del, log(h00) - eta) / n + dot(b, l*b);
		}
	}
	if (p > n) {
		for (int i = 0; i < maxiter; ++i) {
			xdx = ddd*z;
			zeta = solve(xdx + deye, xdx*zeta + delta);
			eta = z*zeta;
			eeta = exp(eta);
			h00 = invcumsum(n, eeta);
			if (tie) {
				tieup(h00, eeta, tie1, tie2);
			}
			h0 = (1 / h00) % dd;
			h = cumsum(h0);
			if (tie) {
				tiedown(h, h0, tie1, tie2);
			}
			eh = eeta%h;
			delta = del - eh;
			ddd = diagmat(eh);
			loss1 = dot(del, log(h00) - eta) / n + dot(zeta, z*zeta);
			if (fabs((loss1 - loss0) / loss0) < cri) {
				break;
			}
			if (i == maxiter) {
				Rcout << "Does not converge" << endl;
			}
			loss0 = loss1;
		}
		b = solve(nl, x.t()*zeta);
	}
	else {
		for (int i = 0; i < maxiter; ++i) {
			xdx = x.t()*ddd*x;
			b = solve(xdx + nl, xdx*b + x.t()*delta);
			eta = x*b;
			eeta = exp(eta);
			h00 = invcumsum(n, eeta);
			if (tie) {
				tieup(h00, eeta, tie1, tie2);
			}
			h0 = (1 / h00) % dd;
			h = cumsum(h0);
			if (tie) {
				tiedown(h, h0, tie1, tie2);
			}
			eh = eeta%h;
			delta = del - eh;
			ddd = diagmat(eh);
			if (ridge) {
				loss1 = dot(del, log(h00) - eta) / n + dot(b, dl%b);
			}
			else {
				loss1 = dot(del, log(h00) - eta) / n + dot(b, l*b);
			}
			if (fabs((loss1 - loss0) / loss0) < cri) {
				break;
			}
			if (i == maxiter) {
				Rcout << "Does not converge" << endl;
			}
			loss0 = loss1;
		}
	}
	return b;
}

vec Coxaen(int n, int p, int ntie, vec b, double lam1, double lam2, vec w, mat x, bool tie, uvec& tie1, uvec& tie2, vec& del, vec&dd, vec& dl, int maxiter, double cri) {
	int act0 = p, act1;
	double bq, loss0, loss1, softer0, softer1;
	mat xsq2 = pow(x, 2);
	vec lw = lam1*w, eta = x*b, eeta = exp(eta), ldl = lam2*dl, zz, ww, bx0, xsq, eh, h00, h0, h, h02, h2;
	uvec actset(p, fill::ones);
	vec::iterator bpoint, lwpoint, ldlpoint;
	uvec::iterator actpoint;
	h00 = invcumsum(n, eeta);
	if (tie) {
		tieup(h00, eeta, tie1, tie2);
	}
	h0 = (1 / h00) % dd;
	h = cumsum(h0);
	if (tie) {
		tiedown(h, h0, tie1, tie2);
	}
	h02 = pow(h0, 2);
	h2 = cumsum(h02);
	if (tie) {
		tiedown(h2, h02, tie1, tie2);
	}
	eh = eeta%h;
	ww = eh - pow(eeta, 2) % h2;
	zz = eta + (del - eh) / ww;
	loss0 = dot(del, log(h00) - eta) / n + dot(lw, abs(b)) + dot(b, ldl%b);
	for (int k = 0; k < maxiter; ++k) {
		act1 = 0;
		actpoint = actset.begin();
		bpoint = b.begin();
		lwpoint = lw.begin();
		ldlpoint = ldl.begin();
		for (int q = 0; q < p; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				softer0 = mean(ww%x.col(q) % (zz - eta + bq*x.col(q)));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (dot(ww, xsq2.col(q)) / n + *ldlpoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx0 = (*bpoint - bq)*x.col(q);
				eta += bx0;
				eeta %= exp(bx0);
				h00 = invcumsum(n, eeta);
				if (tie) {
					tieup(h00, eeta, tie1, tie2);
				}
				h0 = (1 / h00) % dd;
				h = cumsum(h0);
				if (tie) {
					tiedown(h, h0, tie1, tie2);
				}
				h02 = pow(h0, 2);
				h2 = cumsum(h02);
				if (tie) {
					tiedown(h2, h02, tie1, tie2);
				}
				eh = eeta%h;
				ww = eh - pow(eeta, 2) % h2;
				zz = eta + (del - eh) / ww;
			}
			++bpoint;
			++lwpoint;
			++ldlpoint;
			++actpoint;
		}
		act0 -= act1;
		loss1 = dot(del, log(h00) - eta) / n + dot(lw, abs(b)) + dot(b, ldl%b);
		if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
			break;
		}
		if (k == maxiter) {
			Rcout << "Does not converge!" << endl;
		}
		loss0 = loss1;
	}
	return b;
}

vec Coxal(int n, int p, int ntie, vec b, double lam, vec w, mat x, bool tie, uvec& tie1, uvec& tie2, vec& del, vec& dd, int maxiter, double cri) {
	int act0 = p, act1;
	double bq, loss0, loss1, softer0, softer1;
	mat xsq2 = pow(x, 2);
	vec lw = lam*w, eta = x*b, eeta = exp(eta), zz, ww, bx0, xsq, eh, h00, h0, h, h02, h2;
	uvec actset(p, fill::ones);
	vec::iterator bpoint, lwpoint;
	uvec::iterator actpoint;
	h00 = invcumsum(n, eeta);
	if (tie) {
		tieup(h00, eeta, tie1, tie2);
	}
	h0 = (1 / h00) % dd;
	h = cumsum(h0);
	if (tie) {
		tiedown(h, h0, tie1, tie2);
	}
	h02 = pow(h0, 2);
	h2 = cumsum(h0);
	if (tie) {
		tiedown(h2, h02, tie1, tie2);
	}
	eh = eeta%h;
	ww = eh - pow(eeta, 2) % h2;
	zz = eta + (del - eh) / ww;
	loss0 = dot(del, log(h00) - eta) / n + dot(lw, abs(b));
	for (int k = 0; k < maxiter; k++) {
		act1 = 0;
		actpoint = actset.begin();
		bpoint = b.begin();
		lwpoint = lw.begin();
		for (int q = 0; q < p; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				softer0 = mean(ww%x.col(q) % (zz - eta + bq*x.col(q)));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = n*getsign(softer0)*softer1 / dot(ww, xsq2.col(q));
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx0 = (*bpoint - bq)*x.col(q);
				eta += bx0;
				eeta %= exp(bx0);
				h00 = invcumsum(n, eeta);
				if (tie) {
					tieup(h00, eeta, tie1, tie2);
				}
				h0 = (1 / h00) % dd;
				h = cumsum(h0);
				if (tie) {
					tiedown(h, h0, tie1, tie2);
				}
				h02 = pow(h0, 2);
				h2 = cumsum(h02);
				if (tie) {
					tiedown(h2, h02, tie1, tie2);
				}
				eh = eeta%h;
				ww = eh - pow(eeta, 2) % h2;
				zz = eta + (del - eh) / ww;
			}
			++bpoint;
			++actpoint;
			++lwpoint;
		}
		act0 -= act1;
		loss1 = dot(del, log(h00) - eta) / n + dot(lw, abs(b));
		if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
			break;
		}
		if (k == maxiter) {
			Rcout << "Does not converge!" << endl;
		}
		loss0 = loss1;
	}
	return b;
}

vec Coxaagg0(int n, int p, int ntie, vec b, double lam1, double lam2, vec w, mat x, bool tie, uvec& tie1, uvec& tie2, vec& del, vec&dd, mat l, vec dl, int maxiter, double cri) {
	if (lam1 == 0) {
		return Coxnet0(n, p, ntie, b, x, tie, tie1, tie2, del, dd, lam2*(l + diagmat(dl)), maxiter, cri);
	}
	else if (lam2 == 0) {
		return Coxal(n, p, ntie, b, lam1, w, x, tie, tie1, tie2, del, dd, maxiter, cri);
	}
	int p1, p2, act0 = p, act1;
	uvec uncorp = findzerocol(l, p), zelm, nonzelm, actset(p, fill::ones);
	if (all(uncorp)) {
		return Coxaen(n, p, ntie, b, lam1, lam2, w, x, tie, tie1, tie2, del, dd, dl, maxiter, cri);
	}
	double bq, loss0, loss1, softer0, softer1;
	mat ll, lll, xsq2, x2(n, p, fill::zeros);
	vec eta = x*b, eeta = exp(eta), ldl(p), zz, ww, bx0, xsq, eh, h00, h0, h, h02, h2, b2(p, fill::zeros), lw(p, fill::zeros);
	vec::iterator bpoint, lwpoint, ldlpoint;
	uvec::iterator actpoint;
	bool anonze = any(uncorp);
	h00 = invcumsum(n, eeta);
	if (tie) {
		tieup(h00, eeta, tie1, tie2);
	}
	h0 = (1 / h00) % dd;
	h = cumsum(h0);
	if (tie) {
		tiedown(h, h0, tie1, tie2);
	}
	h02 = pow(h0, 2);
	h2 = cumsum(h02);
	if (tie) {
		tiedown(h2, h02, tie1, tie2);
	}
	eh = eeta%h;
	ww = eh - pow(eeta, 2) % h2;
	zz = eta + (del - eh) / ww;
	if (anonze) {
		zelm = find(uncorp);
		nonzelm = find(1 - uncorp);
		p1 = zelm.n_elem;
		p2 = nonzelm.n_elem;
		x2.head_cols(p1) = x.cols(zelm);
		x2.tail_cols(p2) = x.cols(nonzelm);
		b2.head(p1) = b.elem(zelm);
		b2.tail(p2) = b.elem(nonzelm);
		ldl.head(p1) = lam2*dl.elem(zelm);
		ldl.tail(p2) = lam2*dl.elem(nonzelm);
		ll = lam2*l.submat(nonzelm, nonzelm);
		lll = lam2*(l.submat(nonzelm, nonzelm) + diagmat(dl.elem(nonzelm)));
		lw.head(p1) = lam1*w.elem(zelm);
		lw.tail(p2) = lam1*w.elem(nonzelm);
		loss0 = dot(del, log(h00) - eta) / n + dot(lw, abs(b2)) + dot(b2.head(p1), ldl.head(p1) % b2.head(p1)) + dot(b2.tail(p2), lll*b2.tail(p2));
	}
	else {
		p1 = 0;
		p2 = p;
		x2 = x;
		b2 = b;
		ldl = lam2*dl;
		ll = lam2*l;
		lll = lam2*(l + diagmat(dl));
		lw = lam1*w;
		loss0 = dot(del, log(h00) - eta) / n + dot(lw, abs(b2)) + dot(b2, lll*b2);
	}
	xsq2 = pow(x2, 2);
	for (int k = 0; k < maxiter; ++k) {
		act1 = 0;
		actpoint = actset.begin();
		bpoint = b2.begin();
		lwpoint = lw.begin();
		ldlpoint = ldl.begin();
		for (int q = 0; q < p1; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				softer0 = mean(ww%x2.col(q) % (zz - eta + bq*x2.col(q)));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (dot(ww, xsq2.col(q)) / n + *ldlpoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx0 = (*bpoint - bq)*x2.col(q);
				eta += bx0;
				eeta %= exp(bx0);
				h00 = invcumsum(n, eeta);
				if (tie) {
					tieup(h00, eeta, tie1, tie2);
				}
				h0 = (1 / h00) % dd;
				h = cumsum(h0);
				if (tie) {
					tiedown(h, h0, tie1, tie2);
				}
				h02 = pow(h0, 2);
				h2 = cumsum(h02);
				if (tie) {
					tiedown(h2, h02, tie1, tie2);
				}
				eh = eeta%h;
				ww = eh - pow(eeta, 2) % h2;
				zz = eta + (del - eh) / ww;
			}
			++bpoint;
			++lwpoint;
			++ldlpoint;
			++actpoint;
		}
		for (int q = p1; q < p; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				softer0 = mean(ww%x2.col(q) % (zz - eta + bq*x2.col(q))) - dot(ll.col(q - p1), b2.tail(p2));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (dot(ww, xsq2.col(q)) / n + *ldlpoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx0 = (*bpoint - bq)*x2.col(q);
				eta += bx0;
				eeta %= exp(bx0);
				h00 = invcumsum(n, eeta);
				if (tie) {
					tieup(h00, eeta, tie1, tie2);
				}
				h0 = (1 / h00) % dd;
				h = cumsum(h0);
				if (tie) {
					tiedown(h, h0, tie1, tie2);
				}
				h02 = pow(h0, 2);
				h2 = cumsum(h02);
				if (tie) {
					tiedown(h2, h02, tie1, tie2);
				}
				eh = eeta%h;
				ww = eh - pow(eeta, 2) % h2;
				zz = eta + (del - eh) / ww;
			}
			++bpoint;
			++lwpoint;
			++ldlpoint;
			++actpoint;
		}
		act0 -= act1;
		if (anonze) {
			loss1 = dot(del, log(h00) - eta) / n + dot(lw, abs(b2)) + dot(b2.head(p1), ldl.head(p1) % b2.head(p1)) + dot(b2.tail(p2), lll*b2.tail(p2));
			if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
				b.elem(zelm) = b2.head(p1);
				b.elem(nonzelm) = b2.tail(p2);
				break;
			}
		}
		else {
			loss1 = dot(del, log(h00) - eta) / n + dot(lw, abs(b2)) + dot(b2, lll*b2);
			if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
				b = b2;
				break;
			}
		}
		if (k == maxiter) {
			Rcout << "Does not converge!" << endl;
		}
		loss0 = loss1;
	}
	return b;
}

