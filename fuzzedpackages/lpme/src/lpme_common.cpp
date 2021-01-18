#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma ;
using namespace Rcpp ;
using namespace std;

// Fourier transform for error U
Rcpp::NumericVector FuLap(Rcpp::NumericVector t, double sigU){
  return(1.0/(1.0+Rcpp::pow((sigU*t),2)*0.5));
}
Rcpp::NumericVector FuGau(Rcpp::NumericVector t, double sigU){
  return( Rcpp::exp(-Rcpp::pow((sigU*t),2)*0.5) );
}
Rcpp::NumericVector FuLapinv(Rcpp::NumericVector t, double sigU){
  return((1.0+Rcpp::pow((sigU*t),2)*0.5));
}
Rcpp::NumericVector FuGauinv(Rcpp::NumericVector t, double sigU){
  return( Rcpp::exp(Rcpp::pow((sigU*t),2)*0.5) );
}

// Fourier transform for Kernel K
Rcpp::NumericVector FK(Rcpp::NumericVector t){
  return(Rcpp::pow((1.0-Rcpp::pow(t,2)), 8) );
}

// first derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK1(Rcpp::NumericVector t){
  return(-16.0*t*Rcpp::pow((1.0-Rcpp::pow(t,2)), 7));
}

// second derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK2(Rcpp::NumericVector t){
  Rcpp::NumericVector t2 = Rcpp::pow(t,2);
  return( Rcpp::pow((1.0-t2), 6)*(240.0*t2-16.0) );
}

// second-order Kernel K
double K_sec_order(double x){
  double res=0;
  double xi = std::abs(x);
  if(xi<0.2){
    res = 0.1455068+0.0000996*xi+ -0.0084387*std::pow(xi,2);
  }else{
    res = 48.0*std::cos(xi)/(PI*std::pow(xi,4))*(1.0-15.0/std::pow(xi,2)) - 
      144*std::sin(xi)/(PI*std::pow(xi,5))*(2.0-5.0/std::pow(xi,2));
  }
  return(res);
}

// Fourier transform for second-order Kernel K
Rcpp::NumericVector FK_sec_order(Rcpp::NumericVector t){
  return(Rcpp::pow((1.0-Rcpp::pow(t,2)), 3) );
}

// first derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK1_sec_order(Rcpp::NumericVector t){
  return(-6.0*t*Rcpp::pow((1.0-Rcpp::pow(t,2)), 2));
}

// second derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK2_sec_order(Rcpp::NumericVector t){
  return( ( 1.0-Rcpp::pow(t,2) )*( 30.0*Rcpp::pow(t,2)-6.0 ) );
}

// function to generate subvectors a[-(ind1:ind2)] and b[-(ind1:ind2)] and save them to w and y respectively
// void allows more than two values can be returned. For example, w and y are returned here.
void subvecij(const Rcpp::NumericVector& a, const Rcpp::NumericVector& b, int ind1, int ind2, Rcpp::NumericVector& w, Rcpp::NumericVector& y){
  int ny = y.size();
  for (int i=0; i<ny; ++i){
    if(i<ind1) {
      w[i] = a[i];
      y[i] = b[i];
    } else {
      w[i] = a[i+ind2-ind1+1];
      y[i] = b[i+ind2-ind1+1];
    }
  }
}

// function to estimate ghat of x using JASA when U is laplace
void gjasaLap(Rcpp::NumericVector& res, const Rcpp::NumericVector& x, const Rcpp::NumericVector& t, double dt, const Rcpp::NumericVector& W, 
      const Rcpp::NumericVector& Y, double sigU, double h){
  int nt = t.size();
  int n = W.size();
  int nx = x.size();
  Rcpp::NumericVector rehatFW(nt); 
  Rcpp::NumericVector imhatFW(nt);
  Rcpp::NumericVector reYhatFW(nt);
  Rcpp::NumericVector imYhatFW(nt);
  Rcpp::NumericVector FKt = FK(t);
  Rcpp::NumericVector FKt1= FK1(t);
  Rcpp::NumericVector FKt2= FK2(t);
  Rcpp::NumericVector FUt = FuLap(t/h, sigU);
  for (int i=0; i<nt; i++){
    R_CheckUserInterrupt();
    double ti = t[i];
    Rcpp::NumericVector csW = Rcpp::cos(W*ti/h);
    Rcpp::NumericVector snW = Rcpp::sin(W*ti/h);
    rehatFW[i] = Rcpp::sum(csW);
    imhatFW[i] = Rcpp::sum(snW);
    reYhatFW[i] = Rcpp::sum(Y*csW);
    imYhatFW[i] = Rcpp::sum(Y*snW);
  }
  for (int i=0; i<nx; i++){
    R_CheckUserInterrupt();
    Rcpp::NumericVector rext = Rcpp::cos(x[i]*t/h);
    Rcpp::NumericVector imxt = Rcpp::sin(x[i]*t/h);
    Rcpp::NumericVector tmp0 = rehatFW*rext + imhatFW*imxt;
    Rcpp::NumericVector tmp1 = rehatFW*imxt - imhatFW*rext;
    double S0 = Rcpp::sum(tmp0*FKt/FUt)*dt/(n*h*2*PI);
    double S1 = Rcpp::sum(tmp1*FKt1/FUt)*dt/(n*h*2*PI);
    double S2 = -Rcpp::sum(tmp0*FKt2/FUt)*dt/(n*h*2*PI);
    double T0 = Rcpp::sum((reYhatFW*rext + imYhatFW*imxt)*FKt/FUt)*dt/(n*h*2*PI);
    double T1 = Rcpp::sum((reYhatFW*imxt - imYhatFW*rext)*FKt1/FUt)*dt/(n*h*2*PI);
    res[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
  }
}

// function to estimate ghat of x using JASA when U is Gaussian
void gjasaGau(Rcpp::NumericVector& res, const Rcpp::NumericVector& x, const Rcpp::NumericVector& t, double dt, const Rcpp::NumericVector& W, 
      const Rcpp::NumericVector& Y, double sigU, double h){
  int nt = t.size();
  int n = W.size();
  int nx = x.size();
  Rcpp::NumericVector rehatFW(nt); 
  Rcpp::NumericVector imhatFW(nt);
  Rcpp::NumericVector reYhatFW(nt);
  Rcpp::NumericVector imYhatFW(nt);
  Rcpp::NumericVector FKt = FK(t);
  Rcpp::NumericVector FKt1= FK1(t);
  Rcpp::NumericVector FKt2= FK2(t);
  Rcpp::NumericVector FUt = FuGau(t/h, sigU);
  for (int i=0; i<nt; i++){
    R_CheckUserInterrupt();
    double ti = t[i];
    Rcpp::NumericVector csW = Rcpp::cos(W*ti/h);
    Rcpp::NumericVector snW = Rcpp::sin(W*ti/h);
    rehatFW[i] = Rcpp::sum(csW);
    imhatFW[i] = Rcpp::sum(snW);
    reYhatFW[i] = Rcpp::sum(Y*csW);
    imYhatFW[i] = Rcpp::sum(Y*snW);
  }
  for (int i=0; i<nx; i++){
    R_CheckUserInterrupt();
    Rcpp::NumericVector rext = Rcpp::cos(x[i]*t/h);
    Rcpp::NumericVector imxt = Rcpp::sin(x[i]*t/h);
    Rcpp::NumericVector tmp0 = rehatFW*rext + imhatFW*imxt;
    Rcpp::NumericVector tmp1 = -rehatFW*imxt + imhatFW*rext;
    double S0 = Rcpp::sum(tmp0*FKt/FUt)*dt/(n*h*2*PI);
    double S1 = Rcpp::sum(tmp1*FKt1/FUt)*dt/(n*h*2*PI);
    double S2 = -Rcpp::sum(tmp0*FKt2/FUt)*dt/(n*h*2*PI);
    double T0 = Rcpp::sum((reYhatFW*rext + imYhatFW*imxt)*FKt/FUt)*dt/(n*h*2*PI);
    double T1 = Rcpp::sum((-reYhatFW*imxt + imYhatFW*rext)*FKt1/FUt)*dt/(n*h*2*PI);
    res[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
  }
}

// function to estimate ghat of x using OURS when error is Laplace
void gnewLap(Rcpp::NumericVector& ghatofx, const Rcpp::NumericVector& x, const Rcpp::NumericVector& input, const Rcpp::NumericVector& output, 
        double beta, double beta2, const Rcpp::NumericVector& mconst, const Rcpp::NumericVector& Kinput, 
        const Rcpp::NumericVector& W, const Rcpp::NumericVector& Y, double sigU, double h){
  int m = input.size();
  int m_mid = m/2 +1; 
  int n = W.size();
  Rcpp::NumericVector gWinput(m);
  Rcpp::NumericVector fWinput(m);
  double nh = n*h;
  for (int i=0; i<m; ++i){
    R_CheckUserInterrupt();
    double x0 = input[i];
    Rcpp::NumericVector a = (x0-W)/h;
    Rcpp::NumericVector a0(n);
    for (int j=0; j<n; ++j){
      int indx = (int)(round(a[j]/beta+m_mid));
      a0[j] = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0 ;
    }
    Rcpp::NumericVector a1 = a0*a;
    Rcpp::NumericVector a2 = a1*a;
    double S0=Rcpp::sum(a0)/(nh);
    double S1=Rcpp::sum(a1)/(nh);
    double S2=Rcpp::sum(a2)/(nh);
    double T0=Rcpp::sum(Y*a0)/(nh);
    double T1=Rcpp::sum(Y*a1)/(nh);
    gWinput[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
    fWinput[i] = S0;
  }
  
  // Find the support for CF of KK
  int indexl = (int)(round(-1.0/h/beta2+m_mid)); 
  int indexu = (int)(round(1.0/h/beta2+m_mid)); 
  arma::vec support = arma::ones<vec>(m);
  for (int i=0; i<(indexl-1); ++i) {support[i]=0;}
  for (int i=indexu; i<m; ++i) {support[i]=0;}
  
  // Make arma objects
  arma::vec gWin(gWinput.begin(), m, false);
  arma::vec fWin(fWinput.begin(), m, false);
  arma::vec mcon(const_cast<Rcpp::NumericVector&>(mconst).begin(), m, false);
  
  // FFT for fW
  arma::vec Xin=mcon%fWin;
  arma::cx_vec FfW = beta*mcon%arma::fft( Xin )%support;
  
  // FFT for gW*fW
  arma::cx_vec FgWfW=beta*mcon%arma::fft( mcon%gWin%fWin )%support;
  
  // FfU
  Rcpp::NumericVector FfU=FuLap(output, sigU);
  arma::vec FfU2(FfU.begin(), FfU.size(), false);
  
  // inverse FFT to get fX
  arma::cx_vec Fratio=(FfW/FfU2%support);
  arma::cx_vec fXF = mcon/beta%arma::ifft( mcon%Fratio);
  
  // inverse FFT to get gX*fX
  Fratio=(FgWfW/FfU2%support);
  arma::cx_vec gXfX=mcon/beta%arma::ifft( mcon%Fratio);
  
  // estimate of gX 
  arma::vec ghat = arma::real(gXfX)/arma::real(fXF);
  int nx = x.size();
  for (int i=0; i<nx; ++i){
    int ind = (int)(x[i]/beta+m_mid)-1;
    ghatofx[i] = ghat[ind];
  }
}

// function to estimate ghat of x using OURS when error is Gaussian
void gnewGau(Rcpp::NumericVector& ghatofx, const Rcpp::NumericVector& x, const Rcpp::NumericVector& input, const Rcpp::NumericVector& output, 
        double beta, double beta2, const Rcpp::NumericVector& mconst, const Rcpp::NumericVector& Kinput, 
        const Rcpp::NumericVector& W, const Rcpp::NumericVector& Y, double sigU, double h){
  int m = input.size();
  int m_mid = m/2 +1; 
  int n = W.size();
  Rcpp::NumericVector gWinput(m);
  Rcpp::NumericVector fWinput(m);
  double nh = n*h;
  for (int i=0; i<m; ++i){
    R_CheckUserInterrupt();
    double x0 = input[i];
    Rcpp::NumericVector a = (x0-W)/h;
    Rcpp::NumericVector a0(n);
    for (int j=0; j<n; ++j){
      int indx = (int)(round(a[j]/beta+m_mid));
      a0[j] = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0 ;
    }
    Rcpp::NumericVector a1 = a0*a;
    Rcpp::NumericVector a2 = a1*a;
    double S0=Rcpp::sum(a0)/(nh);
    double S1=Rcpp::sum(a1)/(nh);
    double S2=Rcpp::sum(a2)/(nh);
    double T0=Rcpp::sum(Y*a0)/(nh);
    double T1=Rcpp::sum(Y*a1)/(nh);
    gWinput[i] = (T0*S2-T1*S1)/(S0*S2-S1*S1+1e-30);
    fWinput[i] = S0;
  }
  
  // Find the support for CF of KK
  int indexl = (int)(round(-1.0/h/beta2+m_mid)); 
  int indexu = (int)(round(1.0/h/beta2+m_mid)); 
  
  // Make arma objects
  arma::vec gWin(gWinput.begin(), m, false);
  arma::vec fWin(fWinput.begin(), m, false);
  arma::vec mcon(const_cast<Rcpp::NumericVector&>(mconst).begin(), m, false);
  
  // FFT for fW
  arma::vec Xin=mcon%fWin;
  arma::cx_vec FfW = beta*mcon%arma::fft( Xin );
  
  // FFT for gW*fW
  arma::cx_vec FgWfW=beta*mcon%arma::fft( mcon%gWin%fWin );
  
  // FfU
  Rcpp::NumericVector FfU=FuGau(output, sigU);
  arma::vec FfU2(FfU.begin(), FfU.size(), false);
  //(FfU2.elem( arma::find_nonfinite(FfU2>1e-30) )).fill(0);
  
  // inverse FFT to get fX
  arma::cx_vec Fratio=(FfW/FfU2);
  for (int i=0; i<(indexl-1); ++i) {Fratio[i]=0;}
  for (int i=indexu; i<m; ++i) {Fratio[i]=0;}
  arma::cx_vec fXF = mcon/beta%arma::ifft( mcon%Fratio);
  
  // inverse FFT to get gX*fX
  Fratio=(FgWfW/FfU2);
  for (int i=0; i<(indexl-1); ++i) {Fratio[i]=0;}
  for (int i=indexu; i<m; ++i) {Fratio[i]=0;}
  arma::cx_vec gXfX=mcon/beta%arma::ifft( mcon%Fratio);
  
  // estimate of gX 
  arma::vec ghat = arma::real(gXfX)/arma::real(fXF);
  int nx = x.size();
  for (int i=0; i<nx; ++i){
    int ind = (int)(x[i]/beta+m_mid)-1;
    ghatofx[i] = ghat[ind];
  }
}
