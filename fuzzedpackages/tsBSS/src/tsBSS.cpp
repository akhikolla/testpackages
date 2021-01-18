#include "tsBSS.h"
#include <Rcpp.h>
#include <cmath>
//#include <RcppArmadillo.h>
//#include <R.h>
//#include <Rdefines.h>
//#include <Rinternals.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP CCK(SEXP Y, SEXP k) {
  mat X = as<arma::mat>(Y);
  vec tau = as<vec>(k);
  int p = X.n_cols;
  int n = X.n_rows;
  int nk = tau.n_elem;
  mat Ip = eye(p, p);
  mat CK;
  mat C23;
  cube CCK(p, p, p * p * nk);
  cube CC1(p, p, p * p * nk);
  cube CC23(p, p, p * p * nk);
  int ind = 0;

  for (int kk = 0; kk < nk; kk = kk + 1) {
    int K = tau(kk);
    int nK = n - K;
    mat LAMBDA_K = trans(X.rows(0, nK - 1)) * X.rows(K, n - 1) / nK;
    cube Hkt(p,p,nK);

    for (int hh = 0; hh < nK; hh = hh + 1) {
        Hkt.slice(hh) = trans(X.row(hh + K)) * X.row(hh);
    }
    for (int ii = 0; ii < p; ii = ii + 1) {
      for (int jj = 0; jj < p; jj = jj + 1) {
        mat Eij = zeros(p, p);
        Eij(ii, jj) = 1;
        mat Eij2 = Eij + Eij.t();
        mat BK = zeros(p, p);

        for (int bb = 0; bb < nK; bb = bb + 1) {
          BK = BK + trans(Hkt.slice(bb)) * Eij * Hkt.slice(bb);
        }
        BK = BK/nK;

        if(ii == jj) {
          C23 = trans(LAMBDA_K) * Eij2 * LAMBDA_K + Ip;
          CK = BK - C23;
        }
        else {
          C23 = trans(LAMBDA_K) * Eij2 * LAMBDA_K;
          CK = BK - C23;
        }
        CCK.slice(ind) = CK;
        CC1.slice(ind) = BK;
        CC23.slice(ind) = C23;
        ind = ind + 1;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("k") = tau,
                            Rcpp::Named("CCK") = CCK
                            );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP PVCk(SEXP Y, SEXP k) {
  mat X = as<arma::mat>(Y);
  vec tau = as<vec>(k);
  int nk = tau.n_elem;
  int p = X.n_cols;
  int n = X.n_rows;
  mat R = zeros(p, p);
  mat XXcc = zeros(p, p);
  mat Gij = zeros(p, p);
  
  for (int m = 0; m < nk; m = m + 1) {
    int K = tau(m);
    int nK = n - K;
    mat Xt = X.rows(0, nK - 1);
    mat Xttau = X.rows(K, n - 1);
    
    for (int i = 0; i < p; i = i + 1) {
      for (int j = 0; j < p; j = j + 1) {
        vec xij = Xttau.col(i) % Xttau.col(j);
        vec xijc = xij - mean(xij);
        
        cube XX(p, p, nK);
        for (int kk = 0; kk < nK; kk = kk + 1) {
          XX.slice(kk) = trans(Xt.row(kk)) * Xt.row(kk);
        }
        XXcc = mean(XX, 2);
        
        cube XXy(p, p, nK);
        for (int kk = 0; kk < nK; kk = kk + 1) {
          XXy.slice(kk) = (XX.slice(kk) - XXcc) * xijc(kk);
        }
        Gij = mean(XXy, 2);
        R = R + trans(Gij) * Gij;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("k") = tau,
                            Rcpp::Named("R") = R
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP TIK(SEXP Y, SEXP U, SEXP k, SEXP method) {
  mat X = as<arma::mat>(Y);
  mat O = as<arma::mat>(U);
  int tau = as<int>(k);
  int met = as<int>(method);
  int p = X.n_cols;
  int n = X.n_rows;
  mat Xt = X.rows(0, n - 1 - tau);
  mat Xttau = X.rows(tau, n - 1);
  mat Tik = zeros(p, p);

  for (int ii = 0; ii < p; ii = ii + 1) {
    vec oy = Xt * O.col(ii);
    vec oytau = Xttau * O.col(ii);
    double Tik1 = mean(pow(oy % oytau, 2));
    mat Tik2 = mean((Xt % repmat(2*(oy % oytau) % oytau, 1, p)), 0);
    mat Tik3 = mean((Xttau % repmat(2*(oy % oytau) % oy, 1, p)), 0);
	  
	  if (met == 1) {
      Tik.col(ii) = trans(Tik2 + Tik3);
	  }
	  else if (met == 2) {
      Tik.col(ii) = trans((copysign(1.0, Tik1 - 1)) * (Tik2 + Tik3));
	  }
	  else { //if(met == 3)
      Tik.col(ii) = trans((Tik1-1) * (Tik2 + Tik3));
	  }
  }
  return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("k") = tau,
                            Rcpp::Named("Tik") = Tik
                            );
}
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP TIKlc(SEXP Y, SEXP U, SEXP k, SEXP method) {
  mat X = as<arma::mat>(Y);
  mat O = as<arma::mat>(U);
  int tau = as<int>(k);
  int met = as<int>(method);
  int p = X.n_cols;
  int n = X.n_rows;
  mat Xt = X.rows(0, n - 1 - tau);
  mat Xttau = X.rows(tau, n - 1);
  mat Tik = zeros(p, p);

  for (int ii = 0; ii < p; ii = ii + 1) {
    vec oy = Xt * O.col(ii);
    vec oytau = Xttau * O.col(ii);
    vec coy = cosh(oy);
    vec coytau = cosh(oytau);
    vec lcoy = log(coy);
    vec lcoytau = log(coytau);
		double Tik1a = mean(lcoy % lcoytau);
		double Tik1b = mean(lcoy);
		double Tik1c = mean(lcoytau);
    mat Tik2 = mean((Xt % repmat(tanh(oy) % lcoytau, 1, p)), 0);
    mat Tik3 = mean((Xttau % repmat(tanh(oytau) % lcoy, 1, p)), 0);
    mat Tik4 = mean(lcoy)*mean((Xttau % repmat(tanh(oytau), 1, p)), 0);
    mat Tik5 = mean(lcoytau)*mean((Xt % repmat(tanh(oy), 1, p)), 0);

    if (met == 1) {
      Tik.col(ii) = trans(Tik2 + Tik3);
		}
		else if (met == 2) {
      Tik.col(ii) = trans((copysign(1.0, Tik1a - Tik1b*Tik1c)) * (Tik2 + Tik3 - Tik4 - Tik5));
		}
		else {//if(met == 3)
		  Tik.col(ii) = trans((Tik1a - Tik1b*Tik1c) * (Tik2 + Tik3 - Tik4 - Tik5));
		}
  }
  return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("k") = tau,
                            Rcpp::Named("Tik") = Tik
                            );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP TIK1(SEXP Y, SEXP U, SEXP k) {
  mat X = as<arma::mat>(Y);
  mat O = as<arma::mat>(U);
  int tau = as<int>(k);
  int p = X.n_cols;
  int n = X.n_rows;
  mat Xt = X.rows(0, n - 1 - tau);
  mat Xttau = X.rows(tau, n - 1);
  mat Tik = zeros(p,p);
  
  for (int ii = 0; ii < p; ii = ii + 1) {
    vec oy = Xt * O.col(ii);
    vec oytau = Xttau * O.col(ii);
    double Tik1 = mean(oy % oytau);
    mat Tik2 = mean((Xt % repmat(oytau, 1, p)), 0);
    mat Tik3 = mean((Xttau % repmat(oy, 1, p)), 0);
    
    Tik.col(ii) = trans(Tik1*(Tik2 + Tik3));
  }
  return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("k") = tau,
                            Rcpp::Named("Tik") = Tik
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP TSIR(SEXP X, SEXP slices, SEXP k, SEXP h) {
  mat XX = as<arma::mat>(X);
  vec sli = as<vec>(slices);
  int hh = as<int>(h);
  int kk = as<int>(k);
  int n = XX.n_rows;
  int p = XX.n_cols;
  mat sl = sli.subvec(kk, n - 1);
  mat Xk = XX.rows(0, n - kk - 1);
  mat Colm = zeros(hh, p);
  
  for (int i = 0; i < hh; i = i + 1) {
    uvec ind = find(sl == i + 1);
    mat Xksl = Xk.rows(ind);
    Colm.row(i) = mean(Xksl, 0);
  }
  mat RES = cov(Colm);
  return Rcpp::List::create(Rcpp::Named("RES") = RES);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP TSAVE(SEXP X, SEXP slices, SEXP k, SEXP h) {
  mat XX = as<arma::mat>(X);
  vec sli = as<vec>(slices);
  int hh = as<int>(h);
  int kk = as<int>(k);
  int n = XX.n_rows;
  int p = XX.n_cols;
  mat sl = sli.subvec(kk, n - 1);
  mat Xk = XX.rows(0, n - kk - 1);
  mat Ip = eye(p, p);
  cube COV(p, p, hh);
  
  for (int i = 0; i < hh; i = i + 1) {
    uvec ind = find(sl == i + 1);
    mat Xksl = Xk.rows(ind);
    mat CovSl = cov(Xksl);
    mat subt = Ip - CovSl;
    COV.slice(i) = subt * subt.t();
  }
  mat RES = mean(COV, 2);
  return Rcpp::List::create(Rcpp::Named("RES") = RES);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat varExx_k(arma::mat X, arma::vec k) {
  int n = X.n_rows;
  int p = X.n_cols;
  int K = k.n_elem;
  mat vExx_tau = zeros(K, p);
  
  for (int j = 0; j < p; j = j + 1) {
    X.col(j) = X.col(j)/stddev(X.col(j)); //variances scaled to 1;
    for (int i = 0; i < K; i = i + 1) {
      vec vec1 = pow(X(span(0, n - 1 - k(i)), j), 2);
      vec vec2 = pow(X(span(k(i), n - 1), j), 2);
      vExx_tau(i, j) = mean(vec1 % vec2);
      for (int m = 1; m < 21; m = m + 1) {
        vExx_tau(i, j) = vExx_tau(i, j) + 2*((n - m)/n)*(mean(X(span(0, n - m - 1 - k(i)), j) % 
                                                              X(span(k(i), n - m - 1), j) % X(span(m, n - 1 - k(i)), j) %
                                                              X(span(m + k(i), n - 1), j)));
      }
    }
  }
  return vExx_tau;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP lblinM(SEXP X, SEXP k) {
  mat XX = as<arma::mat>(X);
  vec kk = as<vec>(k);
  int n = XX.n_rows;
  int p = XX.n_cols;
  int K = kk.n_elem;
  mat vars = varExx_k(XX, kk);
  vec TS = zeros<vec>(p);
  
  for (int j = 0; j < p; j = j + 1) {
    XX.col(j) = XX.col(j)/stddev(XX.col(j)); //variances scaled to 1;
    for (int i = 0; i < K; i = i + 1) {
      TS(j) = TS(j) + n*pow(mean(XX(span(0, n - 1 - kk(i)), j) % XX(span(kk(i), n - 1), j)), 2)/vars(i, j);
    }
  }
  return Rcpp::List::create(Rcpp::Named("RES") = TS);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP lbsqM(SEXP X, SEXP k) {
  mat XX = as<arma::mat>(X);
  vec kk = as<vec>(k);
  int n = XX.n_rows;
  int p = XX.n_cols;
  int K = kk.n_elem;
  vec TS = zeros<vec>(p);
  
  for (int j = 0; j < p; j = j + 1) {
    XX.col(j) = XX.col(j)/stddev(XX.col(j)); //variances scaled to 1;
    for (int i = 0; i < K; i = i + 1) {
      vec vec1 = pow(XX(span(0, n - 1 - kk(i)), j), 2);
      vec vec2 = pow(XX(span(kk(i), n - 1), j), 2);
      TS(j) = TS(j) + n*pow(mean(vec1 % vec2 - 1), 2)/4;
    }
  }
  return Rcpp::List::create(Rcpp::Named("RES") = TS);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP EIGEN(SEXP X) {  
  arma::mat S = as<arma::mat>(X);
  arma::vec eigvalS;
  arma::mat eigvecS;
  eig_sym(eigvalS,eigvecS,S);
  return Rcpp::List::create(Rcpp::Named("values") = reverse(eigvalS),
                            Rcpp::Named("vectors") = fliplr(eigvecS)
  );
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP PREPBSS(SEXP X, SEXP n){  
  mat x = as<arma::mat>(X); 
  int N = as<int>(n);

  arma::rowvec MEAN = mean(x,0);
  mat xc = x.each_row() - MEAN;
  mat S = xc.t()*xc/(N-1);
  vec eigvalS;
  mat eigvecS;
  eig_sym(eigvalS,eigvecS,S);
  vec SqrtEigvalSI = 1/sqrt(eigvalS);
  mat EVS = diagmat( SqrtEigvalSI );
  mat SqrtI = eigvecS * EVS * eigvecS.t();
  mat y = xc * SqrtI;
  return Rcpp::List::create(Rcpp::Named("Y") = y,
                              Rcpp::Named("X.C") = xc,
                              Rcpp::Named("COV.sqrt.i") = SqrtI,
                              Rcpp::Named("MEAN") = MEAN
    );
}
