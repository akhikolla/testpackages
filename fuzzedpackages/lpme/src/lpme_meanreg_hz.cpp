#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// estimate g function when error is user-assumed
RcppExport SEXP fitNEW( SEXP x_, SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, 
                        SEXP Kinput_, SEXP W_, SEXP Y_, SEXP FfU_, SEXP h_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector x(x_);
  NumericVector input(input_);
  NumericVector output(output_);
  arma::vec mcon = as<arma::vec>(mconst_);
  double beta = as<double>(beta_);
  double beta2 = as<double>(beta2_);
  NumericVector Kinput(Kinput_);
  NumericVector W(W_);
  NumericVector Y(Y_);
  arma::vec FfU = as<arma::vec>(FfU_);
  double h = as<double>(h_);
  
  RNGScope scope;
  
  // start estimating 
  int m = input.size();
  int m_mid = m/2 +1; 
  int n = W.size();
  arma::vec fWin(m);
  arma::vec gWin(m);
  double nh = n*h;
  for (int i=0; i<m; ++i){
    //R_CheckUserInterrupt();
    double x0 = input[i];
    NumericVector a = (x0-W)/h;
    NumericVector a0(n);
    for (int j=0; j<n; ++j){
      int indx = (int)(round(a[j]/beta+m_mid));
      a0[j] = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0 ;
    }
    NumericVector a1 = a0*a;
    NumericVector a2 = a1*a;
    double S0=Rcpp::sum(a0)/(nh);
    double S1=Rcpp::sum(a1)/(nh);
    double S2=Rcpp::sum(a2)/(nh);
    double T0=Rcpp::sum(Y*a0)/(nh);
    double T1=Rcpp::sum(Y*a1)/(nh);
    gWin[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
    fWin[i] = S0;
  }
  
  // Find the support for CF of KK
  int indexl = (int)(round(-1.0/h/beta2+m_mid)); 
  int indexu = (int)(round(1.0/h/beta2+m_mid)); 
  arma::vec support = arma::ones<vec>(m);
  for (int i=0; i<(indexl-1); ++i) {support[i]=0;}
  for (int i=indexu; i<m; ++i) {support[i]=0;}
  
  // FFT for fW
  arma::cx_vec FfW = beta*mcon%arma::fft( mcon%fWin )%support;
  
  // FFT for gW*fW
  arma::cx_vec FgWfW=beta*mcon%arma::fft( mcon%gWin%fWin )%support;
  
  // inverse FFT to get fX
  arma::cx_vec Fratio=(FfW/FfU%support);
  arma::cx_vec fXF = mcon/beta%arma::ifft( mcon%Fratio);
  
  // inverse FFT to get gX*fX
  Fratio=(FgWfW/FfU%support);
  arma::cx_vec gXfX=mcon/beta%arma::ifft( mcon%Fratio);
  
  // estimate of gX 
  arma::vec ghat = arma::real(gXfX)/arma::real(fXF);
  int nx = x.size();
  NumericVector res(nx);
  for (int i=0; i<nx; ++i){
    int ind = (int)(x[i]/beta+m_mid)-1;
    res[i] = ghat[ind];
  }
  
  return List::create(Named("ghat")=res);
  END_RCPP
}

// estimate g function by ignoring measurement error using local polynomial
RcppExport SEXP fitlocpoly( SEXP x_, SEXP beta_, SEXP Kinput_, SEXP X_, SEXP Y_, SEXP h_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
	NumericVector x(x_);
  double beta = as<double>(beta_);
  NumericVector Kinput(Kinput_);
  NumericVector X(X_);
  NumericVector Y(Y_);
  double h = as<double>(h_);
  int nx = x.size();
  int n = X.size();
  int m = Kinput.size();
  int m_mid = m/2 +1; 
  
  // results to save 
  NumericVector ghat(nx);
  NumericVector fhat(nx);
  double nh = n*h;
  for (int i=0; i<nx; ++i){
    R_CheckUserInterrupt();
    NumericVector a = (x[i]-X)/h;
    NumericVector a0(n);
    for (int j=0; j<n; ++j){
      int indx = (int)(round(a[j]/beta+m_mid));
      a0[j] = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0 ;
    }
    NumericVector a1 = a0*a;
    NumericVector a2 = a1*a;
    double S0=Rcpp::sum(a0)/(nh);
    double S1=Rcpp::sum(a1)/(nh);
    double S2=Rcpp::sum(a2)/(nh);
    double T0=Rcpp::sum(Y*a0)/(nh);
    double T1=Rcpp::sum(Y*a1)/(nh);
    ghat[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
    fhat[i] = S0;
  }
  
  return List::create(Named("ghat")=ghat,
                      Named("fhat")=fhat);
	END_RCPP
}

// estimate g function when error is laplace
RcppExport SEXP fitnewLap( SEXP x_, SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP sigU_, SEXP h_) {
  BEGIN_RCPP
  
	// Transfer R variables into C++;
	NumericVector x(x_);
  NumericVector input(input_);
  NumericVector output(output_);
  NumericVector mconst(mconst_);
  double beta = as<double>(beta_);
  double beta2 = as<double>(beta2_);
  NumericVector Kinput(Kinput_);
  NumericVector W(W_);
  NumericVector Y(Y_);
  double sigU = as<double>(sigU_);
  double h = as<double>(h_);
  int nx = x.size();
  
  // results to save 
  NumericVector res(nx);
	
	RNGScope scope;
  
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);
	
  // start estimating 
  gnewLap(res, x, input, output, beta, beta2, mconst, Kinput, W, Y, sigU, h);
  
  return List::create(Named("ghat")=res);
	END_RCPP
}

// estimate g function when error is Gaussian
RcppExport SEXP fitnewGau( SEXP x_, SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP sigU_, SEXP h_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
	NumericVector x(x_);
  NumericVector input(input_);
  NumericVector output(output_);
  NumericVector mconst(mconst_);
  double beta = as<double>(beta_);
  double beta2 = as<double>(beta2_);
  NumericVector Kinput(Kinput_);
  NumericVector W(W_);
  NumericVector Y(Y_);
  double sigU = as<double>(sigU_);
  double h = as<double>(h_);
  int nx = x.size();
  
  // results to save 
  NumericVector res(nx);
	
	RNGScope scope;
  
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);
	
  // start estimating 
  gnewGau(res, x, input, output, beta, beta2, mconst, Kinput, W, Y, sigU, h);
  
  return List::create(Named("ghat")=res);
	END_RCPP
}


// SIMEX bandwidth selection when error is laplace
RcppExport SEXP SIMEXnewLap( SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, SEXP pW_, SEXP pWs_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector input(input_);
  NumericVector output(output_);
  NumericVector mconst(mconst_);
  double beta = as<double>(beta_);
  double beta2 = as<double>(beta2_);
  NumericVector Kinput(Kinput_);
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
        gnewLap(res, xj, input, output, beta, beta2, mconst, Kinput, w, y, sigU, h);
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
        gnewLap(res, xj, input, output, beta, beta2, mconst, Kinput, w, y, sigU, h);
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
RcppExport SEXP SIMEXnewGau( SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_, SEXP Kinput_, 
      SEXP W_, SEXP Y_, SEXP Ws_, SEXP Wss_, SEXP h1_, SEXP h2_, SEXP sigU_, SEXP cumfold_, SEXP pW_, SEXP pWs_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector input(input_);
  NumericVector output(output_);
  NumericVector mconst(mconst_);
  double beta = as<double>(beta_);
  double beta2 = as<double>(beta2_);
  NumericVector Kinput(Kinput_);
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
        gnewGau(res, xj, input, output, beta, beta2, mconst, Kinput, w, y, sigU, h);
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
        gnewGau(res, xj, input, output, beta, beta2, mconst, Kinput, w, y, sigU, h);
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
