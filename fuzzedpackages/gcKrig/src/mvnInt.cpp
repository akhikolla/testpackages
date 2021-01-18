#include <Rmath.h>
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

RcppExport SEXP mvnintGHKOcpp(SEXP mu_, SEXP R_, SEXP lower_, SEXP upper_, SEXP nrep_){
  arma::vec mu = as<arma::vec>(mu_);
  arma::mat R = as<arma::mat>(R_);
  arma::vec lower = as<arma::vec>(lower_);
  arma::vec upper = as<arma::vec>(upper_);
  int nrep = as<int>(nrep_);
  arma::vec mcsample(nrep);

  const double eps = std::numeric_limits<double>::epsilon();
  const double EPS = std::numeric_limits<double>::min();
  int n = R.n_cols;
  arma::mat L = R; arma::vec ap = lower; arma::vec bp = upper;
  arma::vec mup = mu;
  arma::vec y;  y.zeros(n);
  double ctmp_, vtmp_, dem_, am_ = 0, bm_ = 0, de_, ai_, bi_, s_;
  double res_ = 0, prod, MU_, eta_, gamma_, ans_;
  arma::vec z(n);

  double* ctmp = &ctmp_;  double* vtmp = &vtmp_; double* dem = &dem_;
  double* am = &am_; double* bm = &bm_; double* de = &de_;
  double* ai = &ai_; double* bi = &bi_; double* s = &s_;
  double* eta = & eta_; double* res = & res_; double* ans = &ans_;
  double* gamma = & gamma_; double* MU = & MU_;
  int im; arma::mat tmp;  arma::mat mattmp;
  arma::vec d = arma::sqrt(arma::diagvec(L));

  for(int i = 0; i<n; i++){
    if(d(i) > 0){
      L(arma::span::all, i) = L(arma::span::all, i)/d(i);
      L(i, arma::span::all) = L(i, arma::span::all)/d(i);
      ap(i) = ap(i)/d(i); bp(i) = bp(i)/d(i); mup(i) = mup(i)/d(i);
    }
  }

  for(int k=0; k<n; k++){
    im = k;  *vtmp = 0;  *dem = 1.0; *s = 0;
    for(int i=k; i<n; i++){
      if(L(i,i) > eps){
        *ctmp = sqrt(R::fmax2(L(i,i), 0 ));
        if(k == 0){*s = 0;}
        if(i > 0 && k > 0){
          mattmp = L(i, arma::span(0,k-1))*y(arma::span(0,k-1)); *s = mattmp(0,0);
        }
        *ai = (ap(i)-*s-mup(i))/(*ctmp);  *bi = (bp(i)-*s-mup(i))/(*ctmp);
        *de = R::pnorm(*bi,0,1,1,0) - R::pnorm(*ai,0,1,1,0);
        if(*de<=*dem){ *vtmp = *ctmp;  *dem = *de;  *am = *ai;  *bm = *bi; im = i;}
      }
    }

    if(im > k){
      ap.swap_rows(im,k); bp.swap_rows(im,k);
      mup.swap_rows(im,k);
      L(im, im) = L(k,k);
      if(k > 0){
        L(arma::span::all,arma::span(0,k-1)).swap_rows(im,k);
      }
      if( n-im >= 2){
        L(arma::span(im+1,n-1),arma::span::all).swap_cols(im,k);
      }
      if(im - k >= 2){
        tmp = L(arma::span(k+1,im-1),k); L(arma::span(k+1,im-1),k) = L(im,arma::span(k+1,im-1)).t();
        L(im,arma::span(k+1,im-1)) = tmp.t();
      }
    }

    if(n - k >= 2)  L(k,arma::span(k+1,n-1)).zeros();
    if(*vtmp > eps*k){
      L(k,k) = *vtmp;
      for(int i=k+1; i<n; i++){
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
    else{L(arma::span(k,n-1), k).zeros(); y(k) = 0; }
  }
  GetRNGstate();
  for(int ii=0; ii<nrep; ii++) {
    prod = 1.0;
    for(int jj=0; jj<n; jj++) {
      *MU=0;
      for(int kk=0; kk<jj; kk++)  *MU += L(jj,kk)*z(kk);
      *eta = Rcpp::stats::pnorm_1((ap(jj)-*MU-mup(jj))/L(jj,jj),0,1,0);
      *gamma = Rcpp::stats::pnorm_1((bp(jj)-*MU-mup(jj))/L(jj,jj),0,1,0);
      prod *= *gamma-*eta;
      double u = Rf_runif(0, 1);
      *ans = R::fmin2(R::fmax2(eps, u*(*gamma)+(1-u)*(*eta)), 1-eps);
      z(jj) = Rcpp::stats::qnorm_1(*ans,0,1,0);
    }
    *res += R::fmax2(prod, EPS);
    mcsample(ii) = prod;
  }

  *res /= nrep;
  double mcstd = arma::stddev(mcsample)/std::sqrt((double)nrep);
  if(*res == 0 || std::isnan(*res)){
    *res = EPS;
  }
  PutRNGstate();
  return List::create(Named("value") = log(*res),
                      Named("error") = mcstd );
}



RcppExport SEXP mvnintGHKcpp(SEXP mu_, SEXP R_, SEXP lower_, SEXP upper_, SEXP nrep_){

  arma::vec mu = as<arma::vec>(mu_);
  arma::mat R = as<arma::mat>(R_);
  arma::vec lower = as<arma::vec>(lower_);
  arma::vec upper = as<arma::vec>(upper_);
  int nrep = as<int>(nrep_);
  arma::vec mcsample(nrep);
  int dim = R.n_cols;
  arma::vec z(dim);

  arma::mat L = arma::chol(R).t();
  double res_ = 0, prod, MU_, eta_, gamma_, ans_;
  const double EPS = std::numeric_limits<double>::min();
  const double EPS2 = DOUBLE_EPS;
  double* eta = & eta_; double* res = & res_; double* ans = &ans_;
  double* gamma = & gamma_; double* MU = & MU_;
  GetRNGstate();
  for(int i=0;i<nrep;i++) {
    prod = 1.0;
    for(int j=0;j<dim;j++) {
      *MU=0;
      for(int k=0;k<j;k++)  *MU += L(j,k)*z(k);

      *eta = Rcpp::stats::pnorm_0((lower(j)-*MU-mu(j))/L(j,j),1,0);
      *gamma = Rcpp::stats::pnorm_0((upper(j)-*MU-mu(j))/L(j,j),1,0);
      prod *= *gamma-*eta;
      *ans = R::fmin2(R::fmax2(EPS2,*eta + unif_rand()*(*gamma-*eta)), 1-EPS2);
      z(j) = Rcpp::stats::qnorm_0(*ans,1,0);
    }
    *res += R::fmax2(prod, EPS);
    mcsample(i) = prod;
  }
  *res /= nrep;
  double mcstd = arma::stddev(mcsample)/std::sqrt((double)nrep);
  if(*res == 0 || std::isnan(*res)){
    *res = EPS;
  }
  PutRNGstate();
  return List::create(Named("value") = log(*res),
                      Named("error") = mcstd );
}



