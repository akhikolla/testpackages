#include "aag.h"
#include "lrfunc.h"
using namespace arma;
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::vec lrnet(const arma::mat& x, const arma::vec& y, const arma::mat& l) {
	int n = y.n_elem;
	return solve(x.t()*x + n*l, x.t()*y);
}

//[[Rcpp::export]]
arma::vec lraagg(arma::vec b, double lam1, double lam2, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, int maxiter, double cri) {
	if (lam1 == 0) {
		return solve(x.t()*x + y.n_elem*lam2*(l + diagmat(dl)), x.t()*y);
	}
	else if (lam2 == 0) {
		return lral(b, lam1, w, x, y, maxiter, cri);
	}
	
	int p = b.n_elem, n = y.n_elem, act0 = p, act1, p1, p2;
	uvec uncorp = findzerocol(l, p), zelm, nonzelm;
	if (all(uncorp)) {
		return lraen(b, lam1, lam2, w, x, y, dl, maxiter, cri);
	}
	double bq, loss0, loss1, softer0, softer1;
	mat ll, lll, x2(n, p, fill::zeros);
	vec rr = y - x*b, lw(p, fill::zeros), la(p, fill::zeros), b2(p, fill::zeros), ldl, xsq(p, fill::zeros);
	bool anonze = any(uncorp);
	uvec actset(p, fill::ones);
	vec::iterator bpoint, xsqpoint, lwpoint, lapoint;
	uvec::iterator actpoint;
	if (anonze) {
		zelm = find(uncorp);
		nonzelm = find(1 - uncorp);
		p1 = zelm.n_elem;
		p2 = nonzelm.n_elem;
		x2.head_cols(p1) = x.cols(zelm);
		x2.tail_cols(p2) = x.cols(nonzelm);
		b2.head(p1) = b.elem(zelm);
		b2.tail(p2) = b.elem(nonzelm);
		ldl = lam2*dl.elem(zelm);
		ll = lam2*l.submat(nonzelm, nonzelm);
		lll = lam2*(l.submat(nonzelm, nonzelm) + diagmat(dl.elem(nonzelm)));
		xsq.head(p1) = trans(mean(pow(x2.head_cols(p1), 2)));
		xsq.tail(p2) = trans(mean(pow(x2.tail_cols(p2), 2)));
		la.head(p1) = xsq.head(p1) + ldl;
		la.tail(p2) = xsq.tail(p2) + lam2*dl.elem(nonzelm);
		lw.head(p1) = lam1*w.elem(zelm);
		lw.tail(p2) = lam1*w.elem(nonzelm);
		loss0 = dot(rr, rr) / n + dot(lw, abs(b2)) + dot(b2.head(p1), ldl%b2.head(p1)) + dot(b2.tail(p2), lll*b2.tail(p2));
	}
	else {
		x2 = x;
		b2 = b;
		p1 = 0;
		p2 = p;
		ll = lam2*l;
		lll = lam2*(l + diagmat(dl));
		ldl = lam2*dl;
		xsq = trans(mean(pow(x, 2)));
		la = xsq + ldl;
		lw = lam1*w;
		loss0 = dot(rr, rr) / n + dot(lw, abs(b2)) + dot(b2, lll*b2);
	}
	for (int k = 0; k < maxiter; ++k) {
		act1 = 0;
		bpoint = b2.begin();
		xsqpoint = xsq.begin();
		actpoint = actset.begin();
		lapoint = la.begin();
		lwpoint = lw.begin();
		for (int q = 0; q < p1; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				softer0 = *xsqpoint*bq + dot(x2.col(q), rr) / n;
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (*lapoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				rr -= (*bpoint - bq)*x2.col(q);
			}
			++bpoint;
			++actpoint;
			++xsqpoint;
			++lapoint;
			++lwpoint;
		}
		for (int q = p1; q < p; ++q) {
			if (*actpoint) {
				bq = *bpoint;
				softer0 = *xsqpoint*bq + dot(x2.col(q), rr) / n - dot(ll.col(q - p1), b2.tail(p2));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (*lapoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				rr -= (*bpoint - bq)*x2.col(q);
			}
			++bpoint;
			++actpoint;
			++xsqpoint;
			++lapoint;
			++lwpoint;
		}
		act0 -= act1;
		if (anonze) {			
			loss1 = dot(rr, rr) / n + dot(lw, abs(b2)) + dot(b2.head(p1), ldl%b2.head(p1)) + dot(b2.tail(p2), lll*b2.tail(p2));
			if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
				b.elem(zelm) = b2.head(p1);
				b.elem(nonzelm) = b2.tail(p2);
				break;
			}
		}
		else {
			loss1 = dot(rr, rr) / n + dot(lw, abs(b2)) + dot(b2, lll*b2);
			if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
				b = b2;
				break;
			}
		}
		loss0 = loss1;
		if (k == maxiter - 1) {
			Rcout << "Does not converge!" << endl;
		}
	}
	return b;
}



//[[Rcpp::export]]
arma::cube cvlrnet1(int nf, const arma::mat& x, const arma::vec& y, const arma::mat& l1, const arma::mat& l2, arma::vec& lam2, arma::vec& bets, const arma::uvec& cvwhich, bool meas) {
	int k1 = bets.n_elem, k2 = lam2.n_elem, ntr;
	mat xtr, xva;
	cube cvlr(k1, k2, nf, fill::zeros);
	uvec trind, vaind;
	vec m0, ytr, yva;
	vec::iterator bete = bets.end(), lame = lam2.end(), lampoint, betpoint;
	cube::iterator cvlrb = cvlr.begin();
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		ntr = trind.n_elem;
		xtr = x.rows(trind);
		ytr = y(trind);
		xva = x.rows(vaind);
		yva = y(vaind);
		lampoint = lam2.begin();
		while (lampoint != lame) {
			betpoint = bets.begin();
			while (betpoint != bete) {
				m0 = solve(xtr.t()*xtr + ntr*(*lampoint)*(*betpoint*l1 + (1 - *betpoint)*l2), xtr.t()*ytr);
				*cvlrb++ = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
				++betpoint;
			}
			++lampoint;
		}
	}
	return cvlr;
}

//[[Rcpp::export]]
arma::mat cvlrnet1_pal(arma::mat xtr, arma::mat xva, arma::vec ytr, arma::vec yva, arma::mat l1, arma::mat l2, arma::vec lam2, arma::vec bets, bool meas) {
	int k1 = bets.n_elem, k2 = lam2.n_elem, ntr = ytr.n_elem;
	mat cvlr(k1, k2, fill::zeros);
	vec m0;
	vec::iterator lampoint = lam2.begin(), lame = lam2.end(), bete = bets.end(), betpoint;
	mat::iterator cvb = cvlr.begin();
	while (lampoint != lame) {
		betpoint = bets.begin();
		while (betpoint != bete) {
			m0 = solve(xtr.t()*xtr + ntr*(*lampoint)*(*betpoint*l1 + (1 - *betpoint)*l2), xtr.t()*ytr);
			*cvb++ = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
			++betpoint;
		}
		++lampoint;
	}
	return cvlr;
}

//[[Rcpp::export]]
arma::mat cvlrnet2(int nf, const arma::mat& x, const arma::vec& y, const arma::mat& l, arma::vec& lam2, const arma::uvec& cvwhich, bool meas) {
	int k = lam2.n_elem, ntr;
	mat cvlr(k, nf, fill::zeros), xtr, xva;
	uvec trind, vaind;
	vec m0, ytr, yva;
	vec::iterator lame = lam2.end(), lampoint;
	mat::iterator cvpoint = cvlr.begin();
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		ntr = trind.n_elem;
		xtr = x.rows(trind);
		ytr = y(trind);
		xva = x.rows(vaind);
		yva = y(vaind);
		lampoint = lam2.begin();
		while (lampoint != lame) {
			m0 = solve(xtr.t()*xtr + ntr*(*lampoint++)*l, xtr.t()*ytr);
			*cvpoint++ = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
		}
	}
	return cvlr;
}

//[[Rcpp::export]]
arma::vec cvlrnet2_pal(arma::mat xtr, arma::mat xva, arma::vec ytr, arma::vec yva, arma::mat l, arma::vec lam2, bool meas) {
	int k = lam2.n_elem, ntr = ytr.n_elem;
	vec cvlr(k, fill::zeros);
	vec m0;
	vec::iterator lambegin = lam2.begin(), cvpoint = cvlr.begin(), lamend = lam2.end();
	while (lambegin != lamend) {
		m0 = solve(xtr.t()*xtr + ntr*(*lambegin++)*l, xtr.t()*ytr);
		*cvpoint++ = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
	}
	return cvlr;
}

//[[Rcpp::export]]
List cvlraagg(int nf, int dfmax, const arma::vec& b, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, arma::vec& lam1, arma::vec& lam2, const arma::uvec& cvwhich, bool meas, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = b.n_elem, ntr;
	uvec trind, vaind, vind, uind, vvind;
	mat xtr, xva, ll = l + diagmat(dl);
	vec ytr, yva, lpie, apie, m0, m01, mmm;
	cube cvn = p*ones<cube>(k1, k2, nf), cvlr(k1, k2, nf, fill::zeros);
	vec::iterator lam2point, lam1point;
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		xtr = x.rows(trind);
		ytr = y(trind);
		xva = x.rows(vaind);
		yva = y(vaind);
		ntr = ytr.n_elem;
		lam2point = lam2.begin();
		for (int j = 0; j < k2; j++) {
			lam1point = lam1.begin();
			m0 = lraagg(b, *lam1point, *lam2point, w, xtr, ytr, l, dl, maxiter, cri);
			cvn(0, j, i) = sum(m0 != 0);
			cvlr(0, j, i) = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
			if (cvn(0, j, i) < dfmax) {
				for (int q = 1; q < k1; q++) {
					lpie = lrpieabs(ytr, xtr, m0, ntr, *lam2point, ll);
					apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
					++lam1point;
					vind = (lpie >= apie);
					apie = *lam1point*w;
					vvind = find(vind);
					if (any(vind)) {
						mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01(vvind) = mmm;
						lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (any(uind)) {
							while (1) {
								vind = (vind || uind);
								vvind = find(vind);
								mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
								m01 = zeros<vec>(p);
								m01(vvind) = mmm;
								lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
								uind = ((1 - vind) && (lpie > apie));
								if (!any(uind)) {
									m0 = m01;
									break;
								}
							}
						}
						else {
							m0 = m01;
						}
					}
					else {
						m01 = zeros<vec>(p);
						lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (any(uind)) {
							while (1) {
								vind = (vind || uind);
								vvind = find(vind);
								mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
								m01 = zeros<vec>(p);
								m01(vvind) = mmm;
								lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
								uind = ((1 - vind) && (lpie > apie));
								if (!any(uind)) {
									m0 = m01;
									break;
								}
							}
						}
						else {
							m0 = zeros<vec>(p);
						}
					}
					cvn(q, j, i) = sum(m0 != 0);
					cvlr(q, j, i) = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
					if (cvn(q, j, i) >= dfmax) {
						break;
					}
				}
			}
			++lam2point;
		}
	}
	List res;
	res["CV"] = cvlr;
	res["npar"] = cvn;
	return res;
}


//[[Rcpp::export]]
List cvlraagg_pal(int dfmax, arma::vec b, arma::vec w, arma::mat xtr, arma::mat xva, arma::vec ytr, arma::vec yva, arma::mat l, arma::vec dl, arma::vec lam1, arma::vec lam2, bool meas, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = b.n_elem, ntr = ytr.n_elem;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), cvn = p*ones<mat>(k1, k2), cvlr(k1, k2, fill::zeros);
	vec lpie, apie, m0, m01, mmm;
	vec::iterator lam1e = lam1.end(), lam2point = lam2.begin(), lam1point;
	mat::col_iterator cvnb, cvlrb;
	for (int i = 0; i < k2; ++i) {
		lam1point = lam1.begin();
		m0 = lraagg(b, *lam1point, *lam2point, w, xtr, ytr, l, dl, maxiter, cri);
		cvnb = cvn.begin_col(i);
		cvlrb = cvlr.begin_col(i);
		*cvnb = sum(m0 != 0);
		*cvlrb++ = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
		if (*cvnb++ < dfmax) {
			while (lam1point != lam1e) {
				lpie = lrpieabs(ytr, xtr, m0, ntr, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
					m01 = zeros<vec>(p);
					m01(vvind) = mmm;
					lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01(vvind) = mmm;
							lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = m01;
					}
				}
				else {
					m01 = zeros<vec>(p);
					lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01(vvind) = mmm;
							lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = zeros<vec>(p);
					}
				}
				*cvnb = sum(m0 != 0);
				*cvlrb++ = -(meas ? mean(pow(xva*m0 - yva, 2)) : mean(abs(xva*m0 - yva)));
				if (*cvnb++ >= dfmax) {
					break;
				}
			}
		}
		++lam2point;
	}
	List res;
	res["CV"] = cvlr;
	res["npar"] = cvn;
	return res;
}


//[[Rcpp::export]]
arma::cube sslraagg(const arma::vec& b, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, arma::vec& lam1, arma::vec& lam2, const arma::umat& sswhich, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = b.n_elem, ntr = sswhich.n_rows, nsam = sswhich.n_cols;
	uvec vind, uind, vvind;
	mat xtr, ll = l + diagmat(dl);
	vec ytr, lpie, apie, m0, m01, mmm;
	cube ssn(p, k1, k2, fill::zeros);
	vec::iterator lam1point, lam2point;
	for (int i = 0; i < nsam; ++i) {
		xtr = x.rows(sswhich.col(i));
		ytr = y(sswhich.col(i));
		lam2point = lam2.begin();
		for (int j = 0; j < k2; ++j) {
			lam1point = lam1.begin();
			m0 = lraagg(b, *lam1point, *lam2point, w, xtr, ytr, l, dl, maxiter, cri);
			ssn.slice(j).col(0) += conv_to<vec>::from(m0 != 0);
			for (int q = 1; q < k1; ++q) {
				lpie = lrpieabs(ytr, xtr, m0, ntr, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
					m01 = zeros<vec>(p);
					m01.elem(vvind) = mmm;
					lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01.elem(vvind) = mmm;
							lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = m01;
					}
				}
				else {
					m01 = zeros<vec>(p);
					lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01.elem(vvind) = mmm;
							lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = zeros<vec>(p);
					}
				}
				ssn.slice(j).col(q) += conv_to<vec>::from(m0 != 0);
			}
			++lam2point;
		}
	}
	return ssn / nsam;
}

//[[Rcpp::export]]
arma::cube sslraagg_pal(arma::vec b, arma::vec w, arma::mat xtr, arma::vec ytr, arma::mat l, arma::vec dl, arma::vec lam1, arma::vec lam2, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = b.n_elem, ntr = xtr.n_rows;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), dbhat;
	vec lpie, apie, m0, m01, mmm;
	cube ssn(p, k1, k2, fill::zeros);
	vec::iterator lam2point = lam2.begin(), lam1point;
	for (int j = 0; j < k2; ++j) {
		lam1point = lam1.begin();
		m0 = lraagg(b, *lam1point, *lam2point, w, xtr, ytr, l, dl, maxiter, cri);
		ssn.slice(j).col(0) += conv_to<vec>::from(m0 != 0);
		for (int q = 1; q < k1; ++q) {
			lpie = lrpieabs(ytr, xtr, m0, ntr, *lam2point, ll);
			apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
			++lam1point;
			vind = (lpie >= apie);
			apie = *lam1point*w;
			vvind = find(vind);
			if (any(vind)) {
				mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
				m01 = zeros<vec>(p);
				m01.elem(vvind) = mmm;
				lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01.elem(vvind) = mmm;
						lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (!any(uind)) {
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0 = m01;
				}
			}
			else {
				m01 = zeros<vec>(p);
				lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01.elem(vvind) = mmm;
						lpie = lrpieabs(ytr, xtr, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (!any(uind)) {
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0 = zeros<vec>(p);
				}
			}
			ssn.slice(j).col(q) += conv_to<vec>::from(m0 != 0);
		}
		++lam2point;
	}
	return ssn;
}

//[[Rcpp::export]]
List lraagg_search(int dfmax, arma::vec b, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, arma::vec& lam1, arma::vec& lam2, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = b.n_elem, n = x.n_rows;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), searchn = p*ones<mat>(k1, k2), searchlik(k1, k2, fill::zeros);
	vec lpie, apie, m0, m01, mmm;
	cube searchb(p, k1, k2, fill::zeros);
	List res;
	vec::iterator lam2point = lam2.begin(), lam1point;
	mat::col_iterator sb, sl;
	for (int j = 0; j < k2; ++j) {
		sb = searchn.begin_col(j);
		sl = searchlik.begin_col(j);
		lam1point = lam1.begin();
		m0 = lraagg(b, *lam1point, *lam2point, w, x, y, l, dl, maxiter, cri);
		searchb.slice(j).col(0) = m0;
		*sb = sum(m0 != 0);
		*sl++ = -.5*n*log(2 * (datum::pi) / n*sum(pow(y - x*m0, 2))) - .5*n;
		if (*sb++ < dfmax) {
			for (int q = 1; q < k1; ++q) {
				lpie = lrpieabs(y, x, m0, n, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				if ((*lam1point == 0 && *lam2point == 0)&(n <= p)) {
					break;
				}
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
					m01 = zeros<vec>(p);
					m01.elem(vvind) = mmm;
					lpie = lrpieabs(y, x, m01, n, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01.elem(vvind) = mmm;
							lpie = lrpieabs(y, x, m01, n, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = m01;
					}
				}
				else {
					m01 = zeros<vec>(p);
					lpie = lrpieabs(y, x, m01, n, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = lraagg(b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01.elem(vvind) = mmm;
							lpie = lrpieabs(y, x, m01, n, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = zeros<vec>(p);
					}
				}
				searchb.slice(j).col(q) = m0;
				*sb = sum(m0 != 0);
				*sl++ = -.5*n*log(2 * (datum::pi) / n*sum(pow(y - x*m0, 2))) - .5*n;
				if (*sb++ >= dfmax) {
					break;
				}
			}
		}
		
		++lam2point;
	}
	res["coef"] = searchb;
	res["n"] = searchn;
	res["loglik"] = searchlik;
	return res;
}

//[[Rcpp::export]]
List lraagg_search_pal(int dfmax, arma::vec b, arma::vec w, arma::mat x, arma::vec y, arma::mat l, arma::vec dl, arma::vec lam1, double lam2, int maxiter, double cri) {
	int k1 = lam1.n_elem, p = b.n_elem, n = x.n_rows;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), searchb(p, k1, fill::zeros);
	vec searchn = p*ones<vec>(k1), searchlik(k1, fill::zeros), lpie, apie, m0, m01, mmm;
	List res;
	vec::iterator lam1point = lam1.begin(), sb = searchn.begin(), sl = searchlik.begin();
	m0 = lraagg(b, *lam1point, lam2, w, x, y, l, dl, maxiter, cri);
	searchb.col(0) = m0;
	*sb = sum(m0 != 0);
	*sl++ = -.5*n*log(2 * (datum::pi) / n*sum(pow(y - x*m0, 2))) - .5*n;
	if (*sb++ < dfmax) {
		for (int q = 1; q < k1; ++q) {
			lpie = lrpieabs(y, x, m0, n, lam2, ll);
			apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
			++lam1point;
			if ((*lam1point == 0 && lam2 == 0)&(n <= p)) {
				break;
			}
			vind = (lpie >= apie);
			apie = *lam1point*w;
			vvind = find(vind);
			if (any(vind)) {
				mmm = lraagg(b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
				m01 = zeros<vec>(p);
				m01.elem(vvind) = mmm;
				lpie = lrpieabs(y, x, m01, n, lam2, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = lraagg(b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01.elem(vvind) = mmm;
						lpie = lrpieabs(y, x, m01, n, lam2, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (!any(uind)) {
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0 = m01;
				}
			}
			else {
				m01 = zeros<vec>(p);
				lpie = lrpieabs(y, x, m01, n, lam2, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = lraagg(b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01.elem(vvind) = mmm;
						lpie = lrpieabs(y, x, m01, n, lam2, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (!any(uind)) {
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0 = zeros<vec>(p);
				}
			}
			searchb.col(q) = m0;
			*sb = sum(m0 != 0);
			*sl++ = -.5*n*log(2 * (datum::pi) / n*sum(pow(y - x*m0, 2))) - .5*n;
			if (*sb++ >= dfmax) {
				break;
			}
		}
	}
	res["coef"] = searchb;
	res["n"] = searchn;
	res["loglik"] = searchlik;
	return res;
}

//[[Rcpp::export]]
List findlr1se(const arma::vec& b, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, int maxiter, double cri, const arma::mat& tune) {
	int nmod = tune.n_rows, npara, npara0;
	vec m0, parasave;
	rowvec tunesave;
	List res;
	m0 = lraagg(b, tune(0, 0), tune(0, 1), w, x, y, l, dl, maxiter, cri);
	npara0 = sum(m0 != 0);
	npara = npara0;
	parasave = m0;
	tunesave = tune.row(0);
	for (int i = 1; i < nmod; ++i) {
		m0 = lraagg(b, tune(i, 0), tune(i, 1), w, x, y, l, dl, maxiter, cri);
		npara = sum(m0 != 0);
		if (npara < npara0) {
			parasave = m0;
			tunesave = tune.row(i);
			npara0 = npara;
		}
		else {
		  npara = npara0;
		}
		if (npara == 0) {
			break;
		}
	}
	res["tune"] = tunesave;
	res["b"] = parasave;
	return res;
}
