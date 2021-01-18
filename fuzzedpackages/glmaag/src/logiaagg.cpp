#include "aag.h"
#include "logifunc.h"
using namespace arma;
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::vec loginet(double b0, arma::vec b, arma::mat x, const arma::vec& y, arma::mat l, bool intercept, int maxiter, double cri) {
	int n = y.n_elem, p = b.n_elem;
	vec zeta(n, fill::zeros), bx = b0 + x*b, ww, pp, dl, tmpone;
	double loss0, loss1;
	uvec uncorp = findzerocol(l, p), tmp1;
	bool ridge = all(uncorp);
	mat nl = n*l, deye = eye<mat>(n, n), XWX, z;
	if (p > n) {
		z = x*solve(nl, x.t());
	}
	pp = 1 / (1 + exp(-bx));
	ww = pp % (1 - pp);
	dl = l.diag();
	if (p > n) {
		loss0 = mean(ww % pow(((y - pp) / ww), 2)) + dot(zeta, z*zeta);
	}
	else {
		if (ridge) {
			loss0 = mean(ww % pow(((y - pp) / ww), 2)) + dot(b, dl%b);
		}
		else {
			loss0 = mean(ww % pow(((y - pp) / ww), 2)) + dot(b, l*b);
		}
	}
	if (p > n) {
		for (int k = 0; k < maxiter; ++k) {
			XWX = diagmat(ww)*z;
			zeta = solve(XWX + deye, XWX*zeta + y - pp);
			if (intercept) {
				b0 += mean(y - pp) / mean(ww);
			}
			bx = b0 + z*zeta;
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
			loss1 = mean(ww % pow(((y - pp) / ww), 2)) + dot(zeta, z*zeta);
			if (fabs((loss1 - loss0) / loss1) < cri) {
				break;
			}
			if (k == maxiter - 1) {
				Rcout << "Does not converge!" << endl;
			}
			loss0 = loss1;
		}
		b = solve(nl, x.t()*zeta);
	}
	else {
		for (int k = 0; k < maxiter; ++k) {
			XWX = x.t()*diagmat(ww)*x;
			b = solve(XWX + nl, XWX*b + x.t()*(y - pp));
			if (intercept) {
				b0 += mean(y - pp) / mean(ww);
			}
			bx = b0 + x*b;
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
			if (ridge) {
				loss1 = mean(ww % pow(((y - pp) / ww), 2)) + dot(b, dl%b);
			}
			else {
				loss1 = mean(ww % pow(((y - pp) / ww), 2)) + dot(b, l*b);
			}
			if (fabs((loss1 - loss0) / loss1) < cri) {
				break;
			}
			if (k == maxiter - 1) {
				Rcout << "Does not converge!" << endl;
			}
			loss0 = loss1;
		}
	}
	b.insert_rows(0, 1);
	b(0) = b0;
	return b;
}

//[[Rcpp::export]]
arma::vec logiaagg(double b0, arma::vec b, double lam1, double lam2, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, bool intercept, int maxiter, double cri) {
	if (lam1 == 0) {
		return loginet(b0, b, x, y, lam2*(l + diagmat(dl)), intercept, maxiter, cri);
	}
	else if (lam2 == 0) {
		return logial(b0, b, lam1, w, x, y, intercept, maxiter, cri);
	}
	int p = b.n_elem, n = y.n_elem, act0 = p, act1, p1, p2;
	uvec uncorp = findzerocol(l, p), zelm, nonzelm;
	if (all(uncorp)) {
		return logiaen(b0, b, lam1, lam2, w, x, y, dl, intercept, maxiter, cri);
	}
	double del, b00, bq, loss0, loss1, softer0, softer1;
	mat ll, lll, x2(n, p), xsq2;
	vec lw(p, fill::zeros), ww, pp, bx = b0 + x*b, ldl(p, fill::zeros), b2(p, fill::zeros), tmpone;
	bool anonze = any(uncorp);
	uvec actset(p, fill::ones), tmp1;
	vec::iterator bpoint, lwpoint, ldlpoint;
	uvec::iterator actpoint;
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
	if(anonze) {
		zelm = find(uncorp);
		nonzelm = find(1 - uncorp);
		p1 = zelm.n_elem;
		p2 = nonzelm.n_elem;
		x2.head_cols(p1) = x.cols(zelm);
		x2.tail_cols(p2) = x.cols(nonzelm);
		b2.head(p1) = b(zelm);
		b2.tail(p2) = b(nonzelm);
		ldl.head(p1) = lam2*dl(zelm);
		ldl.tail(p2) = lam2*dl(nonzelm);
		ll = lam2*l.submat(nonzelm, nonzelm);
		lll = lam2*l.submat(nonzelm, nonzelm) + diagmat(ldl.tail(p2));
		lw.head(p1) = lam1*w(zelm);
		lw.tail(p2) = lam1*w(nonzelm);
		loss0 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b2)) + dot(b2.head(p1), ldl.head(p1)%b2.head(p1)) + dot(b2.tail(p2), lll*b2.tail(p2));
	}
	else {
		x2 = x;
		b2 = b;
		p1 = 0;
		p2 = p;
		ll = lam2*l;
		lll = lam2*(l + diagmat(dl));
		ldl = lam2*dl;
		lw = lam1*w;
		loss0 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b2)) + dot(b2, lll*b2);
	}
	xsq2 = pow(x2, 2);
	for (int k = 0; k < maxiter; ++k) {
		if (intercept) {
			b00 = mean(y - pp) / mean(ww);
			b0 += b00;
			bx += b00;
		}
		act1 = 0;
		bpoint = b2.begin();
		lwpoint = lw.begin();
		ldlpoint = ldl.begin();
		actpoint = actset.begin();
		for (int q = 0; q < p1; ++q) {
			if (*actpoint) {
				bq = *bpoint;
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
				del = dot(ww, xsq2.col(q)) / n;
				softer0 = del*bq + dot(x2.col(q), y - pp) / n;
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (del + *ldlpoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx += (*bpoint - bq)*x2.col(q);
			}
			++actpoint;
			++bpoint;
			++lwpoint;
			++ldlpoint;
		}
		for (int q = p1; q < p; ++q) {
			if (*actpoint) {
				bq = *bpoint;
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
				del = dot(ww, xsq2.col(q)) / n;
				softer0 = del*bq + dot(x2.col(q), y - pp) / n - dot(ll.col(q - p1), b2.tail(p2));
				softer1 = fabs(softer0) - *lwpoint;
				if (softer1 > 0) {
					*bpoint = getsign(softer0)*softer1 / (del + *ldlpoint);
				}
				else {
					*bpoint = 0;
					*actpoint = 0;
					++act1;
				}
				bx += (*bpoint - bq)*x2.col(q);
			}
			++actpoint;
			++bpoint;
			++lwpoint;
			++ldlpoint;
		}
		act0 -= act1;
		if (anonze) {
			loss1 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b2)) + dot(b2.head(p1), ldl.head(p1) % b2.head(p1)) + dot(b2.tail(p2), lll*b2.tail(p2));
			if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
				b(zelm) = b2.head(p1);
				b(nonzelm) = b2.tail(p2);
				break;
			}
		}
		else {
			loss1 = mean(ww % pow(((y - pp) / ww), 2)) + dot(lw, abs(b2)) + dot(b2, lll*b2);
			if ((act0 == 0) || ((fabs((loss1 - loss0) / loss1) < cri) && (act1 == 0))) {
				b = b2;
				break;
			}
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



//[[Rcpp::export]]
List cvlogiaagg(int nf, int dfmax, arma::vec b, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, const arma::vec& w, arma::vec& lam1, arma::vec& lam2, const arma::uvec& cvwhich, bool intercept, bool meas, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = x.n_cols, ntr, nva;
	double b0, m0i, m01i;
	uvec trind, vaind, vind, uind, vvind;
	mat xtr, xva, ll = l + diagmat(dl);
	vec ytr, yva, score, lpie, apie, m0, m01, mmm;
	cube cvn = p*ones<cube>(k1, k2, nf), cvlogi(k1, k2, nf, fill::zeros);
	vec::iterator lam2point, lam1point;
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		xtr = x.rows(trind);
		ytr = y(trind);
		ntr = trind.n_elem;
		xva = x.rows(vaind);
		yva = y(vaind);
		nva = vaind.n_elem;
		b0 = intercept ? logi1(mean(ytr)) : 0;
		lam2point = lam2.begin();
		for (int j = 0; j < k2; j++) {
			lam1point = lam1.begin();
			mmm = logiaagg(b0, b, *lam1point, *lam2point, w, xtr, ytr, l, dl, intercept, maxiter, cri);
			m0i = mmm(0);
			m0 = mmm.tail(p);
			score = m0i + xva*m0;
			cvn(0, j, i) = sum(m0 != 0);
			cvlogi(0, j, i) = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
			if (cvn(0, j, i) < dfmax) {
				for (int q = 1; q < k1; q++) {
					lpie = logipieabs(ytr, xtr, m0i, m0, ntr, *lam2point, ll);
					apie = (2 * (*(lam1point + 1)) - *lam1point) * w;
					++lam1point;
					vind = (lpie >= apie);
					apie = *lam1point*w;
					vvind = find(vind);
					if (any(vind)) {
						mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
						m01 = zeros<vec>(p);
						m01i = mmm(0);
						m01(vvind) = mmm(span(1, vvind.n_elem));
						lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (any(uind)) {
							while (1) {
								vind = (vind || uind);
								vvind = find(vind);
								mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
								m01 = zeros<vec>(p);
								m01i = mmm(0);
								m01(vvind) = mmm(span(1, vvind.n_elem));
								lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
								uind = ((1 - vind) && (lpie > apie));
								if (!any(uind)) {
									m0i = m01i;
									m0 = m01;
									break;
								}
							}
						}
						else {
							m0i = m01i;
							m0 = m01;
						}
					}
					else {
						m01i = b0;
						m01 = zeros<vec>(p);
						lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (any(uind)) {
							while (1) {
								vind = (vind || uind);
								vvind = find(vind);
								mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
								m01 = zeros<vec>(p);
								m01i = mmm(0);
								m01(vvind) = mmm(span(1, vvind.n_elem));
								lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
								uind = (1 - vind) && (lpie > apie);
								if (!any(uind)) {
									m0i = m01i;
									m0 = m01;
									break;
								}
							}
						}
						else {
							m0 = zeros<vec>(p);
							m0i = b0;
						}
					}
					score = m0i + xva*m0;
					cvn(q, j, i) = sum(m0 != 0);
					cvlogi(q, j, i) = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
					if (cvn(q, j, i) >= dfmax) {
						break;
					}
				}
			}
			++lam2point;
		}
	}
	List res;
	res["CV"] = cvlogi;
	res["npar"] = cvn;
	return res;
}

//[[Rcpp::export]]
arma::cube sslogiaagg(arma::vec b, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, const arma::vec& w, arma::vec& lam1, arma::vec& lam2, const arma::umat& sswhich, bool intercept, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = x.n_cols, ntr = sswhich.n_rows, nsam = sswhich.n_cols;
	double b0, m0i, m01i;
	uvec vind, uind, vvind;
	mat xtr, ll = l + diagmat(dl);
	vec ytr, lpie, apie, m0, m01, mmm;
	cube ssn(p, k1, k2, fill::zeros);
	vec::iterator lam1point, lam2point;
	for (int i = 0; i < nsam; ++i) {
		xtr = x.rows(sswhich.col(i));
		ytr = y(sswhich.col(i));
		b0 = intercept ? logi1(mean(ytr)) : 0;
		lam2point = lam2.begin();
		for (int j = 0; j < k2; ++j) {
			lam1point = lam1.begin();
			mmm = logiaagg(b0, b, *lam1point, *lam2point, w, xtr, ytr, l, dl, intercept, maxiter, cri);
			m0i = mmm(0);
			m0 = mmm.tail(p);
			ssn.slice(j).col(0) += conv_to<vec>::from(m0 != 0);
			for (int q = 1; q < k1; q++) {
				lpie = logipieabs(ytr, xtr, m0i, m0, ntr, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
					m01 = zeros<vec>(p);
					m01i = mmm(0);
					m01(vvind) = mmm(span(1, vvind.n_elem));
					lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
							m01 = zeros<vec>(p);
							m01i = mmm(0);
							m01(vvind) = mmm(span(1, vvind.n_elem));
							lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0i = m01i;
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0i = m01i;
						m0 = m01;
					}
				}
				else {
					m01i = b0;
					m01 = zeros<vec>(p);
					lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
							m01 = zeros<vec>(p);
							m01i = mmm(0);
							m01(vvind) = mmm(span(1, vvind.n_elem));
							lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
							uind = (1 - vind) && (lpie > apie);
							if (!any(uind)) {
								m0i = m01i;
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = zeros<vec>(p);
						m0i = b0;
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
arma::cube sslogiaagg_pal(arma::vec b, arma::mat xtr, arma::vec ytr, arma::mat l, arma::vec dl, arma::vec w, arma::vec lam1, arma::vec lam2, bool intercept, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = b.n_elem, ntr = xtr.n_rows;
	double b0 = intercept ? logi1(mean(ytr)) : 0, m0i, m01i;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), dbhat;
	vec lpie, apie, m0, m01, mmm;
	cube ssn(p, k1, k2, fill::zeros);
	vec::iterator lam2point = lam2.begin(), lam1point;
	for (int j = 0; j < k2; ++j) {
		lam1point = lam1.begin();
		mmm = logiaagg(b0, b, *lam1point, *lam2point, w, xtr, ytr, l, dl, intercept, maxiter, cri);
		m0i = mmm(0);
		m0 = mmm.tail(p);
		ssn.slice(j).col(0) += conv_to<vec>::from(m0 != 0);
		for (int q = 1; q < k1; q++) {
			lpie = logipieabs(ytr, xtr, m0i, m0, ntr, *lam2point, ll);
			apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
			++lam1point;
			vind = (lpie >= apie);
			apie = *lam1point*w;
			vvind = find(vind);
			if (any(vind)) {
				mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
				m01 = zeros<vec>(p);
				m01i = mmm(0);
				m01(vvind) = mmm(span(1, vvind.n_elem));
				lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
						m01 = zeros<vec>(p);
						m01i = mmm(0);
						m01(vvind) = mmm(span(1, vvind.n_elem));
						lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (!any(uind)) {
							m0i = m01i;
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0i = m01i;
					m0 = m01;
				}
			}
			else {
				m01i = b0;
				m01 = zeros<vec>(p);
				lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
						m01 = zeros<vec>(p);
						m01i = mmm(0);
						m01(vvind) = mmm(span(1, vvind.n_elem));
						lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
						uind = (1 - vind) && (lpie > apie);
						if (!any(uind)) {
							m0i = m01i;
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0 = zeros<vec>(p);
					m0i = b0;
				}
			}
			ssn.slice(j).col(q) += conv_to<vec>::from(m0 != 0);
		}
		++lam2point;
	}
	return ssn;
}

//[[Rcpp::export]]
List logiaagg_search(int dfmax, arma::vec b, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, const arma::vec& w, arma::vec& lam1, arma::vec& lam2, bool intercept, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = b.n_elem, n = x.n_rows;
	double b0 = intercept ? logi1(mean(y)) : 0, m0i, m01i;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), searchn = p*ones<mat>(k1, k2), searchlik(k1, k2, fill::zeros), searchint(k1, k2, fill::zeros);
	vec lpie, apie, m0, m01, mmm, score;
	cube searchb(p, k1, k2, fill::zeros);
	List res;
	mat::col_iterator sb, sl, si;
	vec::iterator lam2point = lam2.begin(), lam1point;
	for (int j = 0; j < k2; ++j) {
		sb = searchn.begin_col(j);
		sl = searchlik.begin_col(j);
		si = searchint.begin_col(j);
		lam1point = lam1.begin();
		mmm = logiaagg(b0, b, *lam1point, *lam2point, w, x, y, l, dl, intercept, maxiter, cri);
		m0i = mmm(0);
		m0 = mmm.tail(p);
		searchb.slice(j).col(0) = m0;
		*sb = sum(m0 != 0);
		score = m0i + x*m0;
		*sl++ = dot(y, score) - sum(log(1 + exp(score)));
		*si++ = m0i;
		if (*sb++ < dfmax) {
			for (int q = 1; q < k1; q++) {
				lpie = logipieabs(y, x, m0i, m0, n, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
					m01 = zeros<vec>(p);
					m01i = mmm(0);
					m01(vvind) = mmm(span(1, vvind.n_elem));
					lpie = logipieabs(y, x, m01i, m01, n, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
							m01 = zeros<vec>(p);
							m01i = mmm(0);
							m01(vvind) = mmm(span(1, vvind.n_elem));
							lpie = logipieabs(y, x, m01i, m01, n, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0i = m01i;
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0i = m01i;
						m0 = m01;
					}
				}
				else {
					m01i = b0;
					m01 = zeros<vec>(p);
					lpie = logipieabs(y, x, m01i, m01, n, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
							m01 = zeros<vec>(p);
							m01i = mmm(0);
							m01(vvind) = mmm(span(1, vvind.n_elem));
							lpie = logipieabs(y, x, m01i, m01, n, *lam2point, ll);
							uind = (1 - vind) && (lpie > apie);
							if (!any(uind)) {
								m0i = m01i;
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = zeros<vec>(p);
						m0i = b0;
					}
				}
				searchb.slice(j).col(q) = m0;
				*sb = sum(m0 != 0);
				score = m0i + x*m0;
				*sl++ = dot(y, score) - sum(log(1 + exp(score)));
				*si++ = m0i;
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
	res["interc"] = searchint;
	return res;
}

//[[Rcpp::export]]
List logiaagg_search_pal(int dfmax, arma::vec b, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, const arma::vec& w, arma::vec& lam1, double lam2, bool intercept, int maxiter, double cri) {
	int k1 = lam1.n_elem, p = b.n_elem, n = x.n_rows;
	double b0 = intercept ? logi1(mean(y)) : 0, m0i, m01i;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), searchb = zeros<mat>(p, k1);
	vec searchn = p*ones<vec>(k1), searchlik(k1, fill::zeros), searchint(k1, fill::zeros), lpie, apie, m0, m01, mmm, score;
	List res;
	vec::iterator lam1point = lam1.begin(), sb = searchn.begin(), sl = searchlik.begin(), si = searchint.begin();
	mmm = logiaagg(b0, b, *lam1point, lam2, w, x, y, l, dl, intercept, maxiter, cri);
	m0i = mmm(0);
	m0 = mmm.tail(p);
	searchb.col(0) = m0;
	*sb = sum(m0 != 0);
	score = m0i + x*m0;
	*sl++ = dot(y, score) - sum(log(1 + exp(score)));
	*si++ = m0i;
	if (*sb++ < dfmax) {
		for (int q = 1; q < k1; q++) {
			lpie = logipieabs(y, x, m0i, m0, n, lam2, ll);
			apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
			++lam1point;
			vind = (lpie >= apie);
			apie = *lam1point*w;
			vvind = find(vind);
			if (any(vind)) {
				mmm = logiaagg(b0, b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
				m01 = zeros<vec>(p);
				m01i = mmm(0);
				m01(vvind) = mmm(span(1, vvind.n_elem));
				lpie = logipieabs(y, x, m01i, m01, n, lam2, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = logiaagg(b0, b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
						m01 = zeros<vec>(p);
						m01i = mmm(0);
						m01(vvind) = mmm(span(1, vvind.n_elem));
						lpie = logipieabs(y, x, m01i, m01, n, lam2, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (!any(uind)) {
							m0i = m01i;
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0i = m01i;
					m0 = m01;
				}
			}
			else {
				m01i = b0;
				m01 = zeros<vec>(p);
				lpie = logipieabs(y, x, m01i, m01, n, lam2, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						mmm = logiaagg(b0, b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), y, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
						m01 = zeros<vec>(p);
						m01i = mmm(0);
						m01(vvind) = mmm(span(1, vvind.n_elem));
						lpie = logipieabs(y, x, m01i, m01, n, lam2, ll);
						uind = (1 - vind) && (lpie > apie);
						if (!any(uind)) {
							m0i = m01i;
							m0 = m01;
							break;
						}
					}
				}
				else {
					m0 = zeros<vec>(p);
					m0i = b0;
				}
			}
			searchb.col(q) = m0;
			*sb = sum(m0 != 0);
			score = m0i + x*m0;
			*sl++ = dot(y, score) - sum(log(1 + exp(score)));
			*si++ = m0i;
			if (*sb++ >= dfmax) {
				break;
			}
		}
	}
	res["coef"] = searchb;
	res["n"] = searchn;
	res["loglik"] = searchlik;
	res["interc"] = searchint;
	return res;
}


//[[Rcpp::export]]
List cvlogiaagg_pal(int dfmax, arma::vec b, arma::mat xtr, arma::mat xva, arma::vec ytr, arma::vec yva, arma::mat l, arma::vec dl, arma::vec w, arma::vec lam1, arma::vec lam2, bool intercept, bool meas, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, p = xtr.n_cols, ntr = ytr.n_elem, nva = yva.n_elem;
	double b0, m0i, m01i;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), cvn = p*ones<mat>(k1, k2), cvlogi(k1, k2, fill::zeros);
	vec score, lpie, apie, m0, m01, mmm;
	b0 = intercept ? logi1(mean(ytr)) : 0;
	vec::iterator lam1e = lam1.end(), lam2point = lam2.begin(), lam1point;
	mat::col_iterator cvnb, cvlogib;
	for (int i = 0; i < k2; ++i) {
		lam1point = lam1.begin();
		mmm = logiaagg(b0, b, *lam1point, *lam2point, w, xtr, ytr, l, dl, intercept, maxiter, cri);
		m0i = mmm(0);
		m0 = mmm.tail(p);
		score = m0i + xva*m0;
		cvnb = cvn.begin_col(i);
		cvlogib = cvlogi.begin_col(i);
		*cvnb = sum(m0 != 0);
		*cvlogib++ = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
		if (*cvnb++ < dfmax) {
			while (lam1point != lam1e) {
				lpie = logipieabs(ytr, xtr, m0i, m0, ntr, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
					m01 = zeros<vec>(p);
					m01i = mmm(0);
					m01(vvind) = mmm(span(1, vvind.n_elem));
					lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
							m01 = zeros<vec>(p);
							m01i = mmm(0);
							m01(vvind) = mmm(span(1, vvind.n_elem));
							lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
							uind = ((1 - vind) && (lpie > apie));
							if (!any(uind)) {
								m0i = m01i;
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0i = m01i;
						m0 = m01;
					}
				}
				else {
					m01i = b0;
					m01 = zeros<vec>(p);
					lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							mmm = logiaagg(b0, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), ytr, l.submat(vvind, vvind), dl(vvind), intercept, maxiter, cri);
							m01 = zeros<vec>(p);
							m01i = mmm(0);
							m01(vvind) = mmm(span(1, vvind.n_elem));
							lpie = logipieabs(ytr, xtr, m01i, m01, ntr, *lam2point, ll);
							uind = (1 - vind) && (lpie > apie);
							if (!any(uind)) {
								m0i = m01i;
								m0 = m01;
								break;
							}
						}
					}
					else {
						m0 = zeros<vec>(p);
						m0i = b0;
					}
				}
				score = m0i + xva*m0;
				*cvnb = sum(m0 != 0);
				*cvlogib++ = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
				if (*cvnb++ >= dfmax) {
					break;
				}
			}
		}
		++lam2point;
	}
	List res;
	res["CV"] = cvlogi;
	res["npar"] = cvn;
	return res;
}


//[[Rcpp::export]]
arma::cube cvloginet1(int nf, const arma::mat& x, const arma::vec& y, const arma::vec& b, const arma::mat& l1, const arma::mat& l2, arma::vec& lam2, arma::vec& bets, bool intercept, const arma::uvec& cvwhich, bool meas, int maxiter, double cri) {
	int k1 = bets.n_elem, k2 = lam2.n_elem, p = x.n_cols, nva;
	double b0;
	mat xtr, xva;
	cube cvlogi(k1, k2, nf, fill::zeros);
	uvec trind, vaind;
	vec score, ytr, yva, m0;
	vec::iterator bete = bets.end(), lame = lam2.end(), betpoint, lampoint;
	cube::iterator cvlogib = cvlogi.begin();
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		xtr = x.rows(trind);
		ytr = y(trind);
		xva = x.rows(vaind);
		yva = y(vaind);
		nva = vaind.n_elem;
		b0 = intercept ? logi1(mean(ytr)) : 0;
		lampoint = lam2.begin();
		while (lampoint != lame) {
			betpoint = bets.begin();
			while (betpoint != bete) {
				m0 = loginet(b0, b, xtr, ytr, *lampoint*(*betpoint*l1 + (1 - *betpoint)*l2), intercept, maxiter, cri);
				score = m0(0) + xva*m0(span(1, p));
				*cvlogib++ = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
				++betpoint;
			}
			++lampoint;
		}
	}
	return cvlogi;
}

//[[Rcpp::export]]
arma::mat cvloginet1_pal(arma::mat xtr, arma::mat xva, arma::vec ytr, arma::vec yva, arma::vec b, arma::mat l1, arma::mat l2, arma::vec lam2, arma::vec bets, bool intercept, bool meas, int maxiter, double cri) {
	int k1 = bets.n_elem, k2 = lam2.n_elem, p = xtr.n_cols, nva = yva.n_elem;
	double b0;
	mat cvlogi(k1, k2, fill::zeros);
	vec score, m0;
	vec::iterator lampoint = lam2.begin(), lame = lam2.end(), bete = bets.end(), betpoint;
	mat::iterator cvb = cvlogi.begin();
	b0 = intercept ? logi1(mean(ytr)) : 0;
	while (lampoint != lame) {
		betpoint = bets.begin();
		while (betpoint != bete) {
			m0 = loginet(b0, b, xtr, ytr, *lampoint*(*betpoint*l1 + (1 - *betpoint)*l2), intercept, maxiter, cri);
			score = m0(0) + xva*m0(span(1, p));
			*cvb++ = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
			++betpoint;
		}
		++lampoint;
	}
	return cvlogi;
}


//[[Rcpp::export]]
arma::mat cvloginet2(int nf, const arma::mat& x, const arma::vec& y, const arma::vec& b, const arma::mat& l, arma::vec& lam2, bool intercept, const arma::uvec& cvwhich, bool meas, int maxiter, double cri) {
	int k = lam2.n_elem, p = x.n_cols, nva;
	double b0;
	mat cvlogi(k, nf, fill::zeros), xtr, xva;
	uvec trind, vaind;
	vec score, ytr, yva, m0;
	vec::iterator lame = lam2.end(), lampoint;
	mat::iterator cvpoint = cvlogi.begin();
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		xtr = x.rows(trind);
		ytr = y(trind);
		xva = x.rows(vaind);
		yva = y(vaind);
		nva = vaind.n_elem;
		b0 = intercept ? logi1(mean(ytr)) : 0;
		lampoint = lam2.begin();
		while (lampoint != lame) {
			m0 = loginet(b0, b, xtr, ytr, (*lampoint++)*l, intercept, maxiter, cri);
			score = m0(0) + xva*m0.tail(p);
			*cvpoint++ = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
		}
	}
	return cvlogi;
}

//[[Rcpp::export]]
arma::vec cvloginet2_pal(arma::mat xtr, arma::mat xva, arma::vec ytr, arma::vec yva, arma::vec b, arma::mat l, arma::vec lam2, bool intercept, bool meas, int maxiter, double cri) {
	int k = lam2.n_elem, p = xtr.n_cols, nva = yva.n_elem;
	double b0;
	vec cvlogi(k, fill::zeros), score, m0;
	b0 = intercept ? logi1(mean(ytr)) : 0;
	vec::iterator lamb = lam2.begin(), cvpoint = cvlogi.begin(), lame = lam2.end();
	while (lamb != lame) {
		m0 = loginet(b0, b, xtr, ytr, (*lamb++)*l, intercept, maxiter, cri);
		score = m0(0) + xva*m0.tail(p);
		*cvpoint++ = meas ? (dot(yva, score) - sum(log(1 + exp(score)))) : auc(nva, score, yva);
	}
	return cvlogi;
}

//[[Rcpp::export]]
List findlogi1se(double b0, const arma::vec& b, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::mat& l, const arma::vec& dl, bool intercept, int maxiter, double cri, const arma::mat& tune) {
	int nmod = tune.n_rows, npara, npara0, p = x.n_cols;
	double parasave0;
	vec m0, parasave, mmm;
	rowvec tunesave;
	mmm = logiaagg(b0, b, tune(0, 0), tune(0, 1), w, x, y, l, dl, intercept, maxiter, cri);
	m0 = mmm.tail(p);
	npara0 = sum(m0 != 0);
	npara = npara0;
	parasave0 = mmm(0);
	parasave = m0;
	tunesave = tune.row(0);
	List res;
	for (int i = 1; i < nmod; ++i) {
		mmm = logiaagg(b0, b, tune(i, 0), tune(i, 1), w, x, y, l, dl, intercept, maxiter, cri);		
		m0 = mmm.tail(p);
		npara = sum(m0 != 0);
		if (npara < npara0) {
			parasave0 = mmm(0);
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
	res["b0"] = parasave0;
	res["b"] = parasave;
	return res;
}
