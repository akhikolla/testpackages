#include "aag.h"
#include "Coxfunc.h"
using namespace arma;
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
arma::vec Coxnet(arma::vec b, arma::mat& x, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, arma::vec& del, arma::mat& l, int maxiter, double cri) {
	int n = x.n_rows, p = b.n_elem, ntie;
	vec dd;
	if (tie) {
		ntie = tie1.n_elem;
		dd = getdd(del, tie1, tie2, tie3);
	}
	else {
		ntie = 0;
		dd = del;
	}
	return Coxnet0(n, p, ntie, b, x, tie, tie1, tie2, del, dd, l, maxiter, cri);
}

//[[Rcpp::export]]
arma::vec Coxaagg(arma::vec b, double lam1, double lam2, int ntie, arma::vec w, arma::mat x, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, arma::vec& del, arma::mat l, arma::vec dl, int maxiter, double cri) {
	int n = del.n_elem, p = b.n_elem;
	vec dd;
	if (tie) {
		dd = getdd(del, tie1, tie2, tie3);
	}
	else {
		dd = del;
	}
	return Coxaagg0(n, p, ntie, b, lam1, lam2, w, x, tie, tie1, tie2, del, dd, l, dl, maxiter, cri);
}

//[[Rcpp::export]]
arma::cube cvCoxnet1(int nf, arma::vec b, arma::mat& x, arma::vec& y, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, arma::vec& del, arma::mat& l1, arma::mat& l2, arma::vec& lam2, arma::vec& bets, arma::uvec& cvwhich, bool meas, int maxiter, double cri) {
	bool tietr;
	int k1 = bets.n_elem, k2 = lam2.n_elem, n = y.n_elem, p = b.n_elem, ntie = tie1.n_elem, nutie = 0, ntr, nva, ntietr, j0, j2;
	mat xtr, xva;
	cube cvCox(k1, k2, nf, fill::zeros);
	uvec utie2, trind, vaind, subtie1, gettie, tie1tr0(n, fill::zeros), tie2tr0(n, fill::zeros), tie3tr0(n, fill::zeros), tie1tr, tie2tr, tie3tr;
	vec ytr, yva, deltr, delva, m0, ddtr, time, ytie;
	vec::iterator bete = bets.end(), lame = lam2.end(), lampoint, betpoint;
	cube::iterator cvCoxb = cvCox.begin();
	if (tie) {
		utie2 = unique(tie2);
		nutie = utie2.n_elem;
		ytie = y(utie2);
	}
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		ntr = trind.n_elem;
		xtr = x.rows(trind);
		ytr = y(trind);
		deltr = del.elem(trind);
		nva = vaind.n_elem;
		xva = x.rows(vaind);
		yva = y(vaind);
		delva = del(vaind);
		if (tie) {
			j2 = 0;
			for (int j = 0; j < nutie; ++j) {
				gettie = find(ytr == ytie(j));
				j0 = gettie.n_elem;
				if (j0 > 1) {
					tie1tr0(span(j2, j2 + j0 - 1)) = gettie;
					tie2tr0(span(j2, j2 + j0 - 1)) = gettie(0)*ones<uvec>(j0);
					tie3tr0(span(j2, j2 + j0 - 1)) = gettie(j0 - 1)*ones<uvec>(j0);
					j2 += j0;
				}
				
			}
			tietr = j2 > 1;
			ntietr = tietr ? j2 : 0;
			tie1tr = tie1tr0.head(ntietr);
			tie2tr = tie2tr0.head(ntietr);
			tie3tr = tie3tr0.head(ntietr);
		}
		else {
			tietr = tie;
			ntietr = 0;
			tie1tr = tie1;
			tie2tr = tie2;
			tie3tr = tie3;
		}
		if (tietr) {
			ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
		}
		else {
			ddtr = deltr;
		}
		lampoint = lam2.begin();
		while (lampoint != lame) {
			betpoint = bets.begin();
			while (betpoint != bete) {
				m0 = Coxnet0(ntr, p, ntietr, b, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, *lampoint*(*betpoint*l1 + (1 - *betpoint)*l2), maxiter, cri);
				*cvCoxb++ = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
				++betpoint;
			}
			++lampoint;
		}
	}
	return cvCox;
}


//[[Rcpp::export]]
arma::mat cvCoxnet1_pal(arma::vec b, arma::mat x, arma::mat xtr, arma::mat xva, arma::vec yva, bool tie, arma::uvec tie1, arma::uvec tie2, arma::uvec tie3, bool tietr, arma::uvec tie1tr, arma::uvec tie2tr, arma::uvec tie3tr, arma::vec del, arma::vec deltr, arma::vec delva, arma::mat l1, arma::mat l2, arma::vec lam2, arma::vec bets, bool meas, int maxiter, double cri) {
	int k1 = bets.n_elem, k2 = lam2.n_elem, n = x.n_rows, ntr = xtr.n_rows, nva = yva.n_elem, p = b.n_elem, ntie = tie1.n_elem, ntietr = tie1tr.n_elem;
	mat cvCox(k1, k2, fill::zeros);
	vec m0, ddtr, time;
	vec::iterator lampoint = lam2.begin(), lame = lam2.end(), bete = bets.end(), betpoint;
	mat::iterator cvb = cvCox.begin();
	if (tietr) {
		ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
	}
	else {
		ddtr = deltr;
	}
	while (lampoint != lame) {
		betpoint = bets.begin();
		while (betpoint != bete) {
			m0 = Coxnet0(ntr, p, ntietr, b, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, *lampoint*(*betpoint*l1 + (1 - *betpoint)*l2), maxiter, cri);
			*cvb++ = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
			++betpoint;
		}
		++lampoint;
	}
	return cvCox;
}

//[[Rcpp::export]]
arma::mat cvCoxnet2(int nf, arma::vec& b, arma::mat& x, arma::vec& y, bool tie, int ntie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, arma::vec& del, arma::mat& l, arma::vec& lam, arma::uvec& cvwhich, bool meas, int maxiter, double cri) {
	int k = lam.n_elem, n = y.n_elem, p = b.n_elem, ntietr = 0, nutie = 0, ntr, nva, j0, j2;
	mat cvCox(k, nf, fill::zeros), xtr, xva;
	bool tietr;
	vec ytr, yva, deltr, delva, m0, ddtr, ytie;
	vec::iterator lame = lam.end(), lampoint;
	mat::iterator cvpoint = cvCox.begin();
	uvec utie2, trind, vaind, subtie1, gettie, tie1tr0(n, fill::zeros), tie2tr0(n, fill::zeros), tie3tr0(n, fill::zeros), tie1tr, tie2tr, tie3tr;
	if (tie) {
		utie2 = unique(tie2);
		nutie = utie2.n_elem;
		ytie = y(utie2);
	}
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		vaind = find(cvwhich == i);
		ntr = trind.n_elem;
		xtr = x.rows(trind);
		ytr = y(trind);
		deltr = del(trind);
		nva = vaind.n_elem;
		xva = x.rows(vaind);
		yva = y(vaind);
		delva = del(vaind);
		if (tie) {
			j2 = 0;
			for (int j = 0; j < nutie; ++j) {
				gettie = find(ytr == ytie(j));
				j0 = gettie.n_elem;
				if (j0 > 1) {
					tie1tr0(span(j2, j2 + j0 - 1)) = gettie;
					tie2tr0(span(j2, j2 + j0 - 1)) = gettie(0)*ones<uvec>(j0);
					tie3tr0(span(j2, j2 + j0 - 1)) = gettie(j0 - 1)*ones<uvec>(j0);
					j2 += j0;
				}
				
			}
			tietr = j2 > 1;
			ntietr = tietr ? j2 : 0;
			tie1tr = tie1tr0.head(ntietr);
			tie2tr = tie2tr0.head(ntietr);
			tie3tr = tie3tr0.head(ntietr);
		}
		else {
			tietr = tie;
			ntietr = 0;
			tie1tr = tie1;
			tie2tr = tie2;
			tie3tr = tie3;
		}
		if (tietr) {
			ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
		}
		else {
			ddtr = deltr;
		}
		lampoint = lam.begin();
		while (lampoint != lame) {
			m0 = Coxnet0(ntr, p, ntietr, b, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, (*lampoint++)*l, maxiter, cri);
			*cvpoint++ = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
		}
	}
	return cvCox;
}

//[[Rcpp::export]]
arma::vec cvCoxnet2_pal(arma::vec& b, arma::mat& x, arma::mat& xtr, arma::mat& xva, arma::vec& yva, bool tie, int ntie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, bool tietr, int ntietr, arma::uvec& tie1tr, arma::uvec& tie2tr, arma::uvec& tie3tr, arma::vec& del, arma::vec& deltr, arma::vec& delva, arma::mat& l, arma::vec& lam, bool meas, int maxiter, double cri) {
	int k = lam.n_elem, n = x.n_rows, ntr = xtr.n_rows, nva = yva.n_elem, p = b.n_elem;
	vec cvCox(k, fill::zeros), m0, ddtr;
	vec::iterator lamb = lam.begin(), cvpoint = cvCox.begin(), lame = lam.end();
	if (tietr) {
		ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
	}
	else {
		ddtr = deltr;
	}
	while (lamb != lame) {
		m0 = Coxnet0(ntr, p, ntietr, b, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, (*lamb++)*l, maxiter, cri);
		*cvpoint++ = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
	}
	return cvCox;
}


//[[Rcpp::export]]
List cvCoxaagg(int nf, int dfmax, int ntie, arma::vec& b, arma::vec& w, arma::mat& x, arma::vec& y, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, arma::vec& del, arma::mat& l, arma::vec& dl, arma::vec& lam1, arma::vec& lam2, arma::uvec& cvwhich, bool meas, int maxiter, double cri) {
	bool tietr;
	int k1 = lam1.n_elem, k2 = lam2.n_elem, n = y.n_elem, p = b.n_elem, nutie = 0, ntr = 0, nva = 0, pv, ntietr, j0, j2;
	uvec utie2, trind, vaind, vind, uind, vvind, subtie1, gettie, tie1tr0(n, fill::zeros), tie2tr0(n, fill::zeros), tie3tr0(n, fill::zeros), tie1tr, tie2tr, tie3tr;
	mat ll = l + diagmat(dl), xtr, xva;
	vec ytr, yva, deltr, delva, lpie, apie, m0(p, fill::zeros), m01(p, fill::zeros), mmm, ddtr, time, ytie;
	cube cvn = p*ones<cube>(k1, k2, nf), cvCox(k1, k2, nf, fill::zeros);
	vec::iterator lam2point, lam1point;
	if (tie) {
		utie2 = unique(tie2);
		nutie = utie2.n_elem;
		ytie = y(utie2);
	}
	for (int i = 0; i < nf; ++i) {
		trind = find(cvwhich != i);
		ntr = trind.n_elem;
		xtr = x.rows(trind);
		ytr = y(trind);
		deltr = del(trind);
		vaind = find(cvwhich == i);
		nva = vaind.n_elem;
		xva = x.rows(vaind);
		yva = y(vaind);
		delva = del(vaind);
		if (tie) {
			j2 = 0;
			for (int j = 0; j < nutie; ++j) {
				gettie = find(ytr == ytie(j));
				j0 = gettie.n_elem;
				if (j0 > 1) {
					tie1tr0(span(j2, j2 + j0 - 1)) = gettie;
					tie2tr0(span(j2, j2 + j0 - 1)) = gettie(0)*ones<uvec>(j0);
					tie3tr0(span(j2, j2 + j0 - 1)) = gettie(j0 - 1)*ones<uvec>(j0);
					j2 += j0;
				}
			}
			tietr = j2 > 1;
			ntietr = tietr ? j2 : 0;
			tie1tr = tie1tr0.head(ntietr);
			tie2tr = tie2tr0.head(ntietr);
			tie3tr = tie3tr0.head(ntietr);
		}
		else {
			tietr = tie;
			ntietr = 0;
			tie1tr = tie1;
			tie2tr = tie2;
			tie3tr = tie3;
		}
		if (tietr) {
			ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
		}
		else {
			ddtr = deltr;
		}
		lam2point = lam2.begin();
		for (int j = 0; j < k2; j++) {
			lam1point = lam1.begin();
			m0 = Coxaagg0(ntr, p, ntietr, b, *lam1point, *lam2point, w, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, l, dl, maxiter, cri);
			cvn(0, j, i) = sum(m0 != 0);
			cvCox(0, j, i) = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
			if (cvn(0, j, i) < dfmax) {
				for (int q = 1; q < k1; q++) {
					lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m0, ntr, *lam2point, ll);
					apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
					++lam1point;
					vind = (lpie >= apie);
					apie = *lam1point*w;
					vvind = find(vind);
					if (any(vind)) {
						pv = vvind.n_elem;
						mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01(vvind) = mmm;
						lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (any(uind)) {
							while (1) {
								vind = (vind || uind);
								vvind = find(vind);
								pv = vvind.n_elem;
								mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
								m01 = zeros<vec>(p);
								m01(vvind) = mmm;
								lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
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
						lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
						uind = ((1 - vind) && (lpie > apie));
						if (any(uind)) {
							while (1) {
								vind = (vind || uind);
								vvind = find(vind);
								pv = vvind.n_elem;
								mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
								m01 = zeros<vec>(p);
								m01(vvind) = mmm;
								lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
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
					cvCox(q, j, i) = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
					if (cvn(q, j, i) >= dfmax) {
						break;
					}
				}
			}
			++lam2point;
		}
	}
	List res;
	res["CV"] = cvCox;
	res["npar"] = cvn;
	return res;
}


//[[Rcpp::export]]
arma::cube ssCoxaagg(int ntie, arma::vec& b, arma::vec& w, arma::mat& x, arma::vec& y, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, arma::vec& del, arma::mat& l, arma::vec& dl, arma::vec& lam1, arma::vec& lam2, arma::umat& sswhich, int maxiter, double cri) {
	bool tietr;
	int k1 = lam1.n_elem, k2 = lam2.n_elem, n = y.n_elem, p = b.n_elem, nutie = 0, ntr = sswhich.n_rows, nsam = sswhich.n_cols, pv, ntietr, j0, j2, ntr2;
	uvec utie2, vind, uind, vvind, subtie1, gettie, tie1tr0(n, fill::zeros), tie2tr0(n, fill::zeros), tie3tr0(n, fill::zeros), tie1tr, tie2tr, tie3tr, ss0, ss2;
	mat ll = l + diagmat(dl), xtr;
	vec ytr, deltr, lpie, apie, m0(p, fill::zeros), m01(p, fill::zeros), mmm, ddtr, time, ytie;
	cube ssn(p, k1, k2, fill::zeros);
	vec::iterator lam2point, lam1point;
	if (tie) {
		utie2 = unique(tie2);
		nutie = utie2.n_elem;
		ytie = y(utie2);
	}
	for (int i = 0; i < nsam; ++i) {
		xtr = x.rows(sswhich.col(i));
		ytr = y(sswhich.col(i));
		deltr = del(sswhich.col(i));
		if (deltr(0) == 0) {
			ss0 = find(deltr == 1, 1);
			ss2 = sswhich(span(ss0(0), ntr - 1), i);
			xtr = x.rows(ss2);
			ytr = y(ss2);
			deltr = del(ss2);
			ntr2 = ss2.n_elem;
		}
		else {
			ntr2 = ntr;
		}
		if (tie) {
			j2 = 0;
			for (int j = 0; j < nutie; ++j) {
				gettie = find(ytr == ytie(j));
				j0 = gettie.n_elem;
				if (j0 > 1) {
					tie1tr0(span(j2, j2 + j0 - 1)) = gettie;
					tie2tr0(span(j2, j2 + j0 - 1)) = gettie(0)*ones<uvec>(j0);
					tie3tr0(span(j2, j2 + j0 - 1)) = gettie(j0 - 1)*ones<uvec>(j0);
					j2 += j0;
				}
			}
			tietr = j2 > 1;
			ntietr = tietr ? j2 : 0;
			tie1tr = tie1tr0.head(ntietr);
			tie2tr = tie2tr0.head(ntietr);
			tie3tr = tie3tr0.head(ntietr);
		}
		else {
			tietr = tie;
			ntietr = 0;
			tie1tr = tie1;
			tie2tr = tie2;
			tie3tr = tie3;
		}
		if (tietr) {
			ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
		}
		else {
			ddtr = deltr;
		}
		lam2point = lam2.begin();
		for (int j = 0; j < k2; ++j) {
			lam1point = lam1.begin();
			m0 = Coxaagg0(ntr2, p, ntietr, b, *lam1point, *lam2point, w, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, l, dl, maxiter, cri);
			ssn.slice(j).col(0) += conv_to<vec>::from(m0 != 0);
			for (int q = 1; q < k1; ++q) {
				lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m0, ntr2, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					pv = vvind.n_elem;
					mmm = Coxaagg0(ntr2, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
					m01 = zeros<vec>(p);
					m01(vvind) = mmm;
					lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr2, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							pv = vvind.n_elem;
							mmm = Coxaagg0(ntr2, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01(vvind) = mmm;
							lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr2, *lam2point, ll);
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
					lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							pv = vvind.n_elem;
							mmm = Coxaagg0(ntr2, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01(vvind) = mmm;
							lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr2, *lam2point, ll);
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
arma::cube ssCoxaagg_pal(int ntietr, arma::vec b, arma::vec w, arma::mat xtr, bool tietr, arma::uvec tie1tr, arma::uvec tie2tr, arma::uvec tie3tr, arma::vec deltr, arma::mat l, arma::vec dl, arma::vec lam1, arma::vec lam2, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, ntr = xtr.n_rows, p = b.n_elem, pv;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl);
	vec lpie, apie, m0(p, fill::zeros), m01(p, fill::zeros), mmm, ddtr, time;
	cube ssn(p, k1, k2, fill::zeros);
	vec::iterator lam2point = lam2.begin(), lam1point;
	if (tietr) {
		ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
	}
	else {
		ddtr = deltr;
	}
	for (int j = 0; j < k2; ++j) {
		lam1point = lam1.begin();
		m0 = Coxaagg0(ntr, p, ntietr, b, *lam1point, *lam2point, w, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, l, dl, maxiter, cri);
		ssn.slice(j).col(0) += conv_to<vec>::from(m0 != 0);
		for (int q = 1; q < k1; ++q) {
			lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m0, ntr, *lam2point, ll);
			apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
			++lam1point;
			vind = (lpie >= apie);
			apie = *lam1point*w;
			vvind = find(vind);
			if (any(vind)) {
				pv = vvind.n_elem;
				mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
				m01 = zeros<vec>(p);
				m01(vvind) = mmm;
				lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						pv = vvind.n_elem;
						mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01(vvind) = mmm;
						lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
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
				lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						pv = vvind.n_elem;
						mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01(vvind) = mmm;
						lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
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
List Coxaagg_search(int ntie, int dfmax, arma::vec b, arma::vec w, arma::mat x, bool tie, arma::uvec tie1, arma::uvec tie2, arma::uvec tie3, arma::vec del, arma::mat l, arma::vec dl, arma::vec lam1, arma::vec lam2, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, n = x.n_rows, p = b.n_elem, pv;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), searchn = p*ones<mat>(k1, k2), searchlik(k1, k2, fill::zeros);
	vec lpie, apie, m0(p, fill::zeros), m01(p, fill::zeros), mmm, dd, time;
	cube searchb(p, k1, k2, fill::zeros);
	List res;
	vec::iterator lam2point = lam2.begin(), lam1point;
	mat::col_iterator sb, sl;
	if (tie) {
		dd = getdd(del, tie1, tie2, tie3);
	}
	else {
		dd = del;
	}
	for (int j = 0; j < k2; ++j) {
		sb = searchn.begin_col(j);
		sl = searchlik.begin_col(j);
		lam1point = lam1.begin();
		m0 = Coxaagg0(n, p, ntie, b, *lam1point, *lam2point, w, x, tie, tie1, tie2, del, dd, l, dl, maxiter, cri);
		searchb.slice(j).col(0) = m0;
		*sb = sum(m0 != 0);
		*sl++ = devianceCox(n, x*m0, tie, ntie, tie1, tie2, del);
		if (*sb++ < dfmax) {
			for (int q = 1; q < k1; ++q) {
				lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m0, n, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					pv = vvind.n_elem;
					mmm = Coxaagg0(n, pv, ntie, b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), tie, tie1, tie2, del, dd, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
					m01 = zeros<vec>(p);
					m01(vvind) = mmm;
					lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							pv = vvind.n_elem;
							mmm = Coxaagg0(n, pv, ntie, b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), tie, tie1, tie2, del, dd, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01(vvind) = mmm;
							lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, *lam2point, ll);
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
					lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							pv = vvind.n_elem;
							mmm = Coxaagg0(n, pv, ntie, b(vvind), *lam1point, *lam2point, w(vvind), x.cols(vvind), tie, tie1, tie2, del, dd, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01(vvind) = mmm;
							lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, *lam2point, ll);
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
				*sl++ = devianceCox(n, x*m0, tie, ntie, tie1, tie2, del);
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
List Coxaagg_search_pal(int ntie, int dfmax, arma::vec b, arma::vec w, arma::mat x, bool tie, arma::uvec tie1, arma::uvec tie2, arma::uvec tie3, arma::vec del, arma::mat l, arma::vec dl, arma::vec lam1, double lam2, int maxiter, double cri) {
	int k1 = lam1.n_elem, n = x.n_rows, p = b.n_elem, pv;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), searchb(p, k1, fill::zeros);
	vec searchn = p*ones<vec>(k1), searchlik(k1, fill::zeros), lpie, apie, m0(p, fill::zeros), m01(p, fill::zeros), mmm, dd, time;
	List res;
	vec::iterator lam1point = lam1.begin(), sb = searchn.begin(), sl = searchlik.begin();
	if (tie) {
		dd = getdd(del, tie1, tie2, tie3);
	}
	else {
		dd = del;
	}
	m0 = Coxaagg0(n, p, ntie, b, *lam1point, lam2, w, x, tie, tie1, tie2, del, dd, l, dl, maxiter, cri);
	searchb.col(0) = m0;
	*sb = sum(m0 != 0);
	*sl++ = devianceCox(n, x*m0, tie, ntie, tie1, tie2, del);
	if (*sb++ < dfmax) {
		for (int q = 1; q < k1; ++q) {
			lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m0, n, lam2, ll);
			apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
			++lam1point;
			vind = (lpie >= apie);
			apie = *lam1point*w;
			vvind = find(vind);
			if (any(vind)) {
				pv = vvind.n_elem;
				mmm = Coxaagg0(n, pv, ntie, b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), tie, tie1, tie2, del, dd, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
				m01 = zeros<vec>(p);
				m01(vvind) = mmm;
				lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, lam2, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						pv = vvind.n_elem;
						mmm = Coxaagg0(n, pv, ntie, b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), tie, tie1, tie2, del, dd, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01(vvind) = mmm;
						lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, lam2, ll);
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
				lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, lam2, ll);
				uind = ((1 - vind) && (lpie > apie));
				if (any(uind)) {
					while (1) {
						vind = (vind || uind);
						vvind = find(vind);
						pv = vvind.n_elem;
						mmm = Coxaagg0(n, pv, ntie, b(vvind), *lam1point, lam2, w(vvind), x.cols(vvind), tie, tie1, tie2, del, dd, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
						m01 = zeros<vec>(p);
						m01(vvind) = mmm;
						lpie = Coxpieabs(tie, ntie, tie1, tie2, del, dd, x, m01, n, lam2, ll);
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
			*sl++ = devianceCox(n, x*m0, tie, ntie, tie1, tie2, del);
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
List cvCoxaagg_pal(int dfmax, int ntie, int ntietr, arma::vec b, arma::vec w, arma::mat x, arma::mat xtr, arma::mat xva, arma::vec yva, bool tie, arma::uvec tie1, arma::uvec tie2, arma::uvec tie3, bool tietr, arma::uvec tie1tr, arma::uvec tie2tr, arma::uvec tie3tr, arma::vec del, arma::vec deltr, arma::vec delva, arma::mat l, arma::vec dl, arma::vec lam1, arma::vec lam2, bool meas, int maxiter, double cri) {
	int k1 = lam1.n_elem, k2 = lam2.n_elem, n = x.n_rows, ntr = xtr.n_rows, nva = yva.n_elem, p = b.n_elem, pv;
	uvec vind, uind, vvind;
	mat ll = l + diagmat(dl), cvn = p*ones<mat>(k1, k2), cvCox(k1, k2, fill::zeros);
	vec lpie, apie, m0(p, fill::zeros), m01(p, fill::zeros), mmm, ddtr, time;
	vec::iterator lam1e = lam1.end(), lam2point = lam2.begin(), lam1point;
	mat::col_iterator cvnb, cvCoxb;
	if (tietr) {
		ddtr = getdd(deltr, tie1tr, tie2tr, tie3tr);
	}
	else {
		ddtr = deltr;
	}
	for (int i = 0; i < k2; ++i) {
		lam1point = lam1.begin();
		m0 = Coxaagg0(ntr, p, ntietr, b, *lam1point, *lam2point, w, xtr, tietr, tie1tr, tie2tr, deltr, ddtr, l, dl, maxiter, cri);
		cvnb = cvn.begin_col(i);
		cvCoxb = cvCox.begin_col(i);
		*cvnb = sum(m0 != 0);
		*cvCoxb = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
		if (*cvnb++ < dfmax) {
			while (lam1point != lam1e) {
				lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m0, ntr, *lam2point, ll);
				apie = (2 * (*(lam1point + 1)) - *lam1point)*w;
				++lam1point;
				vind = (lpie >= apie);
				apie = *lam1point*w;
				vvind = find(vind);
				if (any(vind)) {
					pv = vvind.n_elem;
					mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
					m01 = zeros<vec>(p);
					m01.elem(vvind) = mmm;
					lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							pv = vvind.n_elem;
							mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01(vvind) = mmm;
							lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
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
					lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
					uind = ((1 - vind) && (lpie > apie));
					if (any(uind)) {
						while (1) {
							vind = (vind || uind);
							vvind = find(vind);
							pv = vvind.n_elem;
							mmm = Coxaagg0(ntr, pv, ntietr, b(vvind), *lam1point, *lam2point, w(vvind), xtr.cols(vvind), tietr, tie1tr, tie2tr, deltr, ddtr, l.submat(vvind, vvind), dl(vvind), maxiter, cri);
							m01 = zeros<vec>(p);
							m01.elem(vvind) = mmm;
							lpie = Coxpieabs(tietr, ntietr, tie1tr, tie2tr, deltr, ddtr, xtr, m01, ntr, *lam2point, ll);
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
				*cvCoxb++ = meas ? (devianceCox(n, x*m0, tie, ntie, tie1, tie2, del) - devianceCox(ntr, xtr*m0, tietr, ntietr, tie1tr, tie2tr, deltr)) : cidx(nva, xva*m0, yva, delva);
				if (*cvnb++ >= dfmax) {
					break;
				}
			}
		}
		++lam2point;
	}
	List res;
	res["CV"] = cvCox;
	res["npar"] = cvn;
	return res;
}

//[[Rcpp::export]]
List findCox1se(int n, int p, int ntie, arma::vec& b, arma::vec& w, arma::mat& x, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::uvec& tie3, arma::vec& del, arma::mat& l, arma::vec& dl, int maxiter, double cri, arma::mat& tune) {
	int nmod = tune.n_rows, npara, npara0;
	vec m0, parasave, dd;
	rowvec tunesave;
	if (tie) {
		dd = getdd(del, tie1, tie2, tie3);
	}
	else {
		dd = del;
	}
	m0 = Coxaagg0(n, p, ntie, b, tune(0, 0), tune(0, 1), w, x, tie, tie1, tie2, del, dd, l, dl, maxiter, cri);
	npara0 = sum(m0 != 0);
	npara = npara0;
	parasave = m0;
	tunesave = tune.row(0);
	List res;
	for (int i = 1; i < nmod; ++i) {
		m0 = Coxaagg0(n, p, ntie, b, tune(i, 0), tune(i, 1), w, x, tie, tie1, tie2, del, dd, l, dl, maxiter, cri);
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
