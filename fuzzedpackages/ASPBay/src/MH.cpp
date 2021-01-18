// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;

inline double xlogy(double x, double y) {
  return ((x == 0)?0:x*log(y));
}

inline double likelihoodcpp(arma::rowvec f, double psi, IntegerVector S, IntegerVector R) {
  const double f_ab = f(0);
  const double f_Ab = f(1);
  const double f_aB = f(2);
  const double f_AB = f(3);
  const double f_a = f_ab + f_aB;
  const double f_A = 1-f_a;
  const double f_b = f_ab + f_Ab;
  const double f_B = 1-f_b;

  arma::rowvec IBD(3);
  arma::vec GA(3);
  arma::vec rr(3);
  IBD << 0.25 << 0.5 << 0.25 ; 
  GA << f_A*f_A << 2*f_A*f_a << f_a*f_a;
  rr << 1 << psi << psi*psi;

  arma::mat G(3,3);
  G(0,0) = f_AB*f_AB;   G(1,0) = 2*f_AB*f_aB;               G(2,0) = f_aB*f_aB;   
  G(0,1) = 2*f_AB*f_Ab; G(1,1) = 2*(f_aB*f_Ab + f_AB*f_ab); G(2,1) = 2*f_aB*f_ab;
  G(0,2) = f_Ab*f_Ab;   G(1,2) = 2*f_Ab*f_ab;               G(2,2) = f_ab*f_ab;

  //P(GA1 = k, IBD = i | ASP)
  arma::mat HA = (GA * IBD);
  HA.col(0) = HA.col(0) % rr;
  HA.col(1) = HA.col(1) % rr;
  HA.col(2) = HA.col(2) % rr;

  arma::mat L(3,3);
  L(0,0) = sum(GA%rr);  L(0,1) = f_A + psi*f_a;               L(0,2) = 1;
  L(1,0) = L(0,0);      L(1,1) = 0.5*(f_A + psi*(1+psi*f_a)); L(1,2) = psi;
  L(2,0) = L(0,0);      L(2,1) = psi*(f_A + psi*f_a);         L(2,2) = psi*psi;
 
  HA = HA % L;
  HA = HA/sum(sum(HA)); 
  
  //P(GB1 = k, IBD = i | ASP)
  arma::mat HB(3,3);
  arma::vec sG = sum(G,1); 
  HB(0,0) = sum( HA.col(0) % G.col(0) / sG);
  HB(0,1) = sum( HA.col(1) % G.col(0) / sG);
  HB(0,2) = sum( HA.col(2) % G.col(0) / sG);

  HB(1,0) = sum( HA.col(0) % G.col(1) / sG);
  HB(1,1) = sum( HA.col(1) % G.col(1) / sG);
  HB(1,2) = sum( HA.col(2) % G.col(1) / sG);

  HB(2,0) = sum( HA.col(0) % G.col(2) / sG);
  HB(2,1) = sum( HA.col(1) % G.col(2) / sG);
  HB(2,2) = sum( HA.col(2) % G.col(2) / sG);
 
  HB = HB / sum(sum(HB));

  return 2*xlogy(S(0),f_B) + xlogy(S(1), 2*f_B*f_b) + 2*xlogy(S(2), f_b) 
       + xlogy(R(0), HB(0,0)) + xlogy(R(1), HB(0,1)) + xlogy(R(2), HB(0,2)) 
       + xlogy(R(3), HB(1,0)) + xlogy(R(4), HB(1,1)) + xlogy(R(5), HB(1,2))
       + xlogy(R(6), HB(2,0)) + xlogy(R(7), HB(2,1)) + xlogy(R(8), HB(2,2)) ;
         
}

// [[Rcpp::export]]
arma::mat MHcpp(int N, int thin, IntegerVector S, IntegerVector R, double sd_freq, double sd_psi, arma::rowvec p0, double psi_prior) {

  arma::mat X(1 + (N-1)/thin,5);
  // initial values
  X.row(0) = p0;
  arma::rowvec f = X.row(0).cols(0,3);
  double psi = X(0,4);

  // initial likelihood
  double slog = sum(log(f));
  double L = likelihoodcpp(f, psi, S, R) - 0.5*psi_prior*log(psi)*log(psi) - 0.5*slog;

  for(int i = 1; i < N; i++) {
    // candidate values
    arma::rowvec f1 = f;
    f1(0) *= exp(R::rnorm(0,sd_freq));
    f1(1) *= exp(R::rnorm(0,sd_freq));
    f1(2) *= exp(R::rnorm(0,sd_freq));
    f1(3) *= exp(R::rnorm(0,sd_freq));
    f1 = f1/sum(f1);
    double psi1 = psi*exp(R::rnorm(0,sd_psi));
    // likelihood
    double slog1 = sum(log(f1));
    double L1 = likelihoodcpp(f1, psi1, S, R) - 0.5*psi_prior*log(psi1)*log(psi1) - 0.5*slog1;
    // accept
    if( log(R::runif(0,1)) < L1 - L + slog1 - slog) {
      f = f1; psi = psi1; L = L1; slog = slog1;
    }
    // fill X
    if(i % thin == 0) {
      X.row(i/thin).cols(0,3) = f;
      X(i/thin,4) = psi;
    }
  }
  return X;
}
