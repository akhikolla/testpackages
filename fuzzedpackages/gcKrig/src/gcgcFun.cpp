#include <Rmath.h>
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

/*

The Function ghkgcmr is inherited from package gcmr written by Masarotto and Varin
with little modification (antithetic variable).

Function ghkgcmr2 is the function that provides a reordering of variables.

Function ldgc is to calculate gaussian copula log-density.

mvnintGHKOcpp is re-ordered version function similar as
 ghkgcmr2 (except that arguments are different)
 that is used to compute integral values.

mvnintGHKcpp is original version without reordering to compute integral values.
*/

RcppExport SEXP ghkgcmr(SEXP R_, SEXP pdf_, SEXP cdf_, SEXP nrep_){

  arma::mat R = as<arma::mat>(R_);
  arma::vec pdf = as<arma::vec>(pdf_);
  arma::vec cdf = as<arma::vec>(cdf_);
  int nrep = as<int>(nrep_);
  int  r , i, j, dim = R.n_cols;
  arma::mat L = arma::chol(R, "lower");
  arma::vec llik(dim); arma::vec lower(dim); arma::vec upper(dim);
  double yi, lk , plow, pup, mw, s2, mwold, biasold, xr1, xr2, randu,
  EPS = sqrt(DOUBLE_EPS), EPS1=1-EPS;
  arma::vec xr(dim);
  arma::vec w(nrep); w.ones();
  mwold = 1 ; biasold = 0 ;
  GetRNGstate();
  for ( i=0 ; i<dim ;  i++) {
    lower(i) = R::qnorm(R::fmax2(EPS,R::fmin2(EPS1,cdf(i) - pdf(i))),0,1,1,0);
    upper(i) = R::qnorm(R::fmax2(EPS, R::fmin2(EPS1,cdf(i))), 0,1,1,0);
    mw = s2  = 0.0 ;
    for ( r=0 ; r < nrep ; r++ ) {
      yi=0;
      for ( j=0;  j < i ; j++ ) {
        yi += L(i,j)*(xr(j)) ;
      }
      plow = Rcpp::stats::pnorm_0( (lower(i)-yi)/L(i,i), 1 , 0 ) ;
      pup = Rcpp::stats::pnorm_0( (upper(i)-yi)/L(i,i), 1 , 0 ) ;
      randu = unif_rand();
      xr1 = Rcpp::stats::qnorm_0(R::fmax2(EPS,R::fmin2(EPS1,plow+randu*(pup-plow))), 1, 0);
      xr2 = Rcpp::stats::qnorm_0(R::fmax2(EPS,R::fmin2(EPS1,pup+randu*(plow-pup))), 1, 0);
      xr(i) = (xr1+xr2)/2;
      w(r) /= mwold ;
      lk = R::fmax2(EPS,R::fmin2(EPS1,pup-plow));
      w(r) *= lk ;
      mw += w(r) ;
      s2 += w(r)*w(r);
    }
    mw /= nrep ;
    s2 = ((s2/nrep)-mw*mw)/(2*(nrep-1)*mw*mw) ;
    llik(i) = log(mw) + s2 - biasold ;
    mwold = mw ;
    biasold = s2 ;
  }
  return wrap(arma::accu(llik));
  PutRNGstate();
}






RcppExport SEXP ghkgcmr2(SEXP R_, SEXP pdf_, SEXP cdf_, SEXP nrep_){

  arma::mat R = as<arma::mat>(R_);
  arma::vec pdf = as<arma::vec>(pdf_);
  arma::vec cdf = as<arma::vec>(cdf_);
  int nrep = as<int>(nrep_);
  int dim = R.n_cols;
  int i, j, r ,ii, k;
  arma::vec llik(dim);
  double yi, lk , plow, pup, mw, s2, mwold, biasold, xr1, xr2, randu,
  EPS = sqrt(DOUBLE_EPS), EPS1=1-EPS;
  const double eps = std::numeric_limits<double>::epsilon();
  arma::vec xr(dim), ap(dim), bp(dim), w(nrep); w.ones();
  mwold = 1 ; biasold = 0 ;
  arma::mat L = R;

  arma::vec y;  y.zeros(dim);
  double ctmp_, vtmp_, dem_, am_ = 0, bm_ = 0, de_, ai_, bi_, s_;

  double* ctmp = &ctmp_;  double* vtmp = &vtmp_; double* dem = &dem_;
  double* am = &am_; double* bm = &bm_; double* de = &de_;
  double* ai = &ai_; double* bi = &bi_; double* s = &s_;
  int im; arma::mat tmp;  arma::mat mattmp;

  for(ii = 0; ii<dim; ii++){
    ap(ii) = R::qnorm(R::fmax2(EPS,R::fmin2(EPS1,cdf(ii) - pdf(ii))),0,1,1,0);
    bp(ii) = R::qnorm(R::fmax2(EPS, R::fmin2(EPS1,cdf(ii))), 0,1,1,0);
  }

  for( k=0; k<dim; k++){
    im = k;  *vtmp = 0;  *dem = 1.0; *s = 0;
     for(i=k; i<dim; i++){
       if(L(i,i) > eps){
        *ctmp = sqrt(R::fmax2(L(i,i), 0 ));
        if(k == 0){*s = 0;}
        if(i > 0 && k > 0){
          mattmp = L(i, arma::span(0,k-1))*y(arma::span(0,k-1)); *s = mattmp(0,0);
        }
        *ai = (ap(i)-*s)/(*ctmp);  *bi = (bp(i)-*s)/(*ctmp);
        *de = R::pnorm(*bi,0,1,1,0) - R::pnorm(*ai,0,1,1,0);
        if(*de<=*dem){ *vtmp = *ctmp;  *dem = *de;  *am = *ai;  *bm = *bi; im = i;}
      }
    }

    if(im > k){
      ap.swap_rows(im,k); bp.swap_rows(im,k);
      L(im, im) = L(k,k);
      if(k > 0){
        L(arma::span::all,arma::span(0,k-1)).swap_rows(im,k);
      }
      if( dim-im >= 2){
        L(arma::span(im+1,dim-1),arma::span::all).swap_cols(im,k);
      }
      if(im - k >= 2){
        tmp = L(arma::span(k+1,im-1),k); L(arma::span(k+1,im-1),k) = L(im,arma::span(k+1,im-1)).t();
        L(im,arma::span(k+1,im-1)) = tmp.t();
      }
    }

    if(dim - k >= 2)  L(k,arma::span(k+1,dim-1)).zeros();
    if(*vtmp > eps*k){
      L(k,k) = *vtmp;
      for(int i=k+1; i<dim; i++){
        L(i,k) = L(i,k)/(*vtmp);
        L(i, arma::span(k+1,i)) = L(i, arma::span(k+1,i))-L(i,k)*L(arma::span(k+1,i),k).t();
      }
      if(std::abs(*dem) > eps){
        y(k) = (R::dnorm(*am,0,1,0) - R::dnorm(*bm,0,1,0) )/ (*dem);
      }else{
        if(*am < -10) {y(k) = *bm;}
        else if(*bm > 10) {y(k) = *am;}
        else{y(k) = (*am+*bm)/2;}
      }
    }
    else{L(arma::span(k,dim-1), k).zeros(); y(k) = 0; }
  }

  GetRNGstate();
  for ( i=0 ; i<dim ;  i++) {
    mw = s2  = 0.0 ;
    for ( r=0 ; r < nrep ; r++ ) {
      yi=0;
      for ( j=0;  j < i ; j++ ) {
        yi += L(i,j)*(xr(j)) ;
      }
      plow = Rcpp::stats::pnorm_0( (ap(i)-yi)/L(i,i), 1 , 0 ) ;
      pup = Rcpp::stats::pnorm_0( (bp(i)-yi)/L(i,i), 1 , 0 ) ;
      randu = unif_rand();
      xr1 = Rcpp::stats::qnorm_0(R::fmax2(EPS,R::fmin2(EPS1,plow+randu*(pup-plow))), 1, 0);
      xr2 = Rcpp::stats::qnorm_0(R::fmax2(EPS,R::fmin2(EPS1,pup+randu*(plow-pup))), 1, 0);
      xr(i) = (xr1+xr2)/2;
      w(r) /= mwold ;
      lk = R::fmax2(EPS,R::fmin2(EPS1,pup-plow));
      w(r) *= lk ;
      mw += w(r) ;
      s2 += w(r)*w(r);
    }
    mw /= nrep ;
    s2 = ((s2/nrep)-mw*mw)/(2*(nrep-1)*mw*mw) ;
    llik(i) = log(mw) + s2 - biasold ;
    mwold = mw ;
    biasold = s2 ;
  }
  return wrap(arma::accu(llik));
  PutRNGstate();
}




RcppExport SEXP ldgc(SEXP u_, SEXP R_) {

    arma::rowvec u = as<arma::rowvec>(u_);
    arma::mat R = as<arma::mat>(R_);

    int dimn = R.n_rows;
    arma::mat RI = arma::inv(R)-arma::speye(dimn, dimn);
    arma::rowvec phi(dimn);
    for(int i=0; i<dimn; i++){
        phi[i] = Rcpp::stats::qnorm_0(u[i], 1, 0);
    }
    arma::mat outtmp = phi*RI*phi.t();
    double logdet = sum(arma::log(arma::eig_sym(R)));
    double out = -0.5*(logdet+outtmp(0));
    if(std::isinf(out) == true || std::isnan(out) == true){
        out = -9.532493e+14;
    }
    return wrap(out);
}
