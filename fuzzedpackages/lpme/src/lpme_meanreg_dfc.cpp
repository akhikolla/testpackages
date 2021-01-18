#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// estimate g function JASA when error is laplace
RcppExport SEXP fitjasaLap( SEXP x_, SEXP h_, SEXP W_, SEXP Y_, SEXP sigU_, SEXP dt_, SEXP t_) {
	BEGIN_RCPP
  
	// Transfer R variables into C++;
	NumericVector x(x_);
  NumericVector W(W_);
  NumericVector Y(Y_);
  double sigU = as<double>(sigU_);
  double h = as<double>(h_);
  double dt = as<double>(dt_);
  NumericVector t(t_);
  int nx = x.size();
  
  // results to save 
  NumericVector res(nx);
	
	RNGScope scope;
  
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);
	
  // start estimating 
  gjasaLap(res, x, t, dt, W, Y, sigU, h);
  
  return List::create(Named("ghat")=res);
	END_RCPP
}

// estimate g function JASA when error is Gaussian
RcppExport SEXP fitjasaGau( SEXP x_, SEXP h_, SEXP W_, SEXP Y_, SEXP sigU_, SEXP dt_, SEXP t_) {
  BEGIN_RCPP
  
	// Transfer R variables into C++;
	NumericVector x(x_);
  NumericVector W(W_);
  NumericVector Y(Y_);
  double sigU = as<double>(sigU_);
  double h = as<double>(h_);
  double dt = as<double>(dt_);
  NumericVector t(t_);
  int nx = x.size();
  
  // results to save 
  NumericVector res(nx);
	
	RNGScope scope;
  
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);
	
  // start estimating 
  gjasaGau(res, x, t, dt, W, Y, sigU, h);
  
  return List::create(Named("ghat")=res);
	END_RCPP
}

// SIMEX bandwidth selection when error is laplace
RcppExport SEXP SIMEXjasaLap(SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, 
                  SEXP pW_, SEXP pWs_, SEXP dt_, SEXP t_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector Y(Y_);
  NumericMatrix Ws(Ws_);
  NumericMatrix Wss(Wss_);
  NumericVector pW(pW_);
  NumericMatrix pWs(pWs_);
  const IntegerVector cumfold(cumfold_);
  const double sigU = as<double>(sigU_);
  NumericVector h1(h1_);
  NumericVector h2(h2_);
  NumericVector t(t_);
  double dt = as<double>(dt_);
  int B = Ws.ncol();
  int n = W.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  
  // temp variable
  NumericVector CVh1(nh1);
  NumericVector CVh2(nh2);
  
  RNGScope scope;
	
  // start SIMEX 
  for (int i=0; i<nh1; ++i){
    double h = h1[i];
    Rprintf( "Evaluating CV1: i=%d\n", i+1 );
    //Rprintf( "%f\n", h );
    NumericVector CV(B);
    for (int b=0; b<B; ++b){
      R_CheckUserInterrupt();
      NumericVector Wstar = Ws(_,b);
      NumericVector gW(n);
      for (int j=1; j<cumfold.size(); ++j){
        int ind1 = cumfold[j-1];
        int ind2 = cumfold[j]-1;
        NumericVector xj = W[Range(ind1,ind2)];
        NumericVector res(xj.size());
        int nrem = n - xj.size(); 
        NumericVector w(nrem);
        NumericVector y(nrem);
        subvecij(Wstar, Y, ind1, ind2, w, y);
        gjasaLap(res, xj, t, dt, w, y, sigU, h);
        gW[Range(ind1,ind2)] = res;
      }
      CV[b] = Rcpp::mean( Rcpp::pow((Y-gW),2)*pW );
      //Rprintf( "%f\n", CV[b] );
    }
    CVh1[i] = Rcpp::mean(CV);
    //Rprintf( "%f\n", CVh1[i] );
  }
  
  for (int i=0; i<nh2; ++i){
    double h = h2[i];
    Rprintf( "Evaluating CV2: i=%d\n", i+1 );
    //Rprintf( "%f\n", h );
    NumericVector CV(B);
    for (int b=0; b<B; ++b){
      R_CheckUserInterrupt();
      NumericVector Wstarstar = Wss(_,b);
      NumericVector Wstar = Ws(_,b);
      NumericVector pWsb = pWs(_,b);
      NumericVector gW(n);
      for (int j=1; j<cumfold.size(); ++j){
        int ind1 = cumfold[j-1];
        int ind2 = cumfold[j]-1;
        NumericVector xj = Wstar[Range(ind1,ind2)];
        NumericVector res(xj.size());
        int nrem = n - xj.size(); 
        NumericVector w(nrem);
        NumericVector y(nrem);
        subvecij(Wstarstar, Y, ind1, ind2, w, y);
        gjasaLap(res, xj, t, dt, w, y, sigU, h);
        gW[Range(ind1,ind2)] = res;
      }
      CV[b] = Rcpp::mean( Rcpp::pow((Y-gW),2)*pW );
      //Rprintf( "%f\n", CV[b] );
    }
    CVh2[i] = Rcpp::mean(CV);
    //Rprintf( "%f\n", CVh2[i] );
  }
  
  return List::create(Named("h1")=h1,
                    Named("CVh1")=CVh1,
                    Named("h2")=h2,
                    Named("CVh2")=CVh2);
	END_RCPP
}

// SIMEX bandwidth selection when error is Gaussian
RcppExport SEXP SIMEXjasaGau(SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, 
                  SEXP pW_, SEXP pWs_, SEXP dt_, SEXP t_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector Y(Y_);
  NumericMatrix Ws(Ws_);
  NumericMatrix Wss(Wss_);
  NumericVector pW(pW_);
  NumericMatrix pWs(pWs_);
  const IntegerVector cumfold(cumfold_);
  const double sigU = as<double>(sigU_);
  NumericVector h1(h1_);
  NumericVector h2(h2_);
  NumericVector t(t_);
  double dt = as<double>(dt_);
  int B = Ws.ncol();
  int n = W.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  
  // temp variable
  NumericVector CVh1(nh1);
  NumericVector CVh2(nh2);
  
  RNGScope scope;
  
  // start SIMEX 
  for (int i=0; i<nh1; ++i){
    double h = h1[i];
    Rprintf( "Evaluating CV1: i=%d\n", i+1 );
    //Rprintf( "%f\n", h );
    NumericVector CV(B);
    for (int b=0; b<B; ++b){
      R_CheckUserInterrupt();
      NumericVector Wstar = Ws(_,b);
      NumericVector gW(n);
      for (int j=1; j<cumfold.size(); ++j){
        int ind1 = cumfold[j-1];
        int ind2 = cumfold[j]-1;
        NumericVector xj = W[Range(ind1,ind2)];
        NumericVector res(xj.size());
        int nrem = n - xj.size(); 
        NumericVector w(nrem);
        NumericVector y(nrem);
        subvecij(Wstar, Y, ind1, ind2, w, y);
        gjasaGau(res, xj, t, dt, w, y, sigU, h);
        gW[Range(ind1,ind2)] = res;
      }
      CV[b] = Rcpp::mean( Rcpp::pow((Y-gW),2)*pW );
      //Rprintf( "%f\n", CV[b] );
    }
    CVh1[i] = Rcpp::mean(CV);
    //Rprintf( "%f\n", CVh1[i] );
  }
  
  for (int i=0; i<nh2; ++i){
    double h = h2[i];
    Rprintf( "Evaluating CV2: i=%d\n", i+1 );
    //Rprintf( "%f\n", h );
    NumericVector CV(B);
    for (int b=0; b<B; ++b){
      R_CheckUserInterrupt();
      NumericVector Wstarstar = Wss(_,b);
      NumericVector Wstar = Ws(_,b);
      NumericVector pWsb = pWs(_,b);
      NumericVector gW(n);
      for (int j=1; j<cumfold.size(); ++j){
        int ind1 = cumfold[j-1];
        int ind2 = cumfold[j]-1;
        NumericVector xj = Wstar[Range(ind1,ind2)];
        NumericVector res(xj.size());
        int nrem = n - xj.size(); 
        NumericVector w(nrem);
        NumericVector y(nrem);
        subvecij(Wstarstar, Y, ind1, ind2, w, y);
        gjasaGau(res, xj, t, dt, w, y, sigU, h);
        gW[Range(ind1,ind2)] = res;
      }
      CV[b] = Rcpp::mean( Rcpp::pow((Y-gW),2)*pW );
      //Rprintf( "%f\n", CV[b] );
    }
    CVh2[i] = Rcpp::mean(CV);
    //Rprintf( "%f\n", CVh2[i] );
  }
  
  return List::create(Named("h1")=h1,
                    Named("CVh1")=CVh1,
                    Named("h2")=h2,
                    Named("CVh2")=CVh2);
	END_RCPP
}

