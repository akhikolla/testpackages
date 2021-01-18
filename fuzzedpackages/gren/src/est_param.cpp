#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

/* This function estimates new variational parameters*/
// [[Rcpp::export]]
List est_param(arma::mat xr, arma::mat xu, arma::vec kappa, arma::vec m, int n, int p, arma::vec ciold, 
               double phi, arma::vec chiold, double lambda2, arma::vec lambdag, 
               arma::vec lambdagold, bool intercept, bool unpen, bool posterior, 
               bool elbo, bool start) {
  
  arma::vec lambda2vec = lambda2*lambdag;
  int u = 0;
  int r = xr.n_cols;
  int nvars = xr.n_cols;
  arma::mat x = xr;
  if(intercept==true && unpen==false) {
    u = 1;
    nvars += 1;
    arma::vec evec(n,arma::fill::ones);
    x = join_rows(evec, xr);
  } else if(unpen==true){
    u = xu.n_cols;
    nvars += xu.n_cols;
    x = join_rows(xu,xr);
  }
  
  // initialize c, chi, dsigma and mu
  double elboval;
  arma::vec ci(n);
  arma::vec chi(r);
  arma::vec dsigma(nvars);
  arma::mat sigma(nvars,nvars);
  arma::vec mu(nvars);
  
  // used in all situations
  arma::mat xrt = xr.t();
  arma::vec h(r);
  arma::vec om(n);
  if(start) {
    // if we are calculating starting values we need the diagonal weight 
    // matrix with wi = pi/(1-pi) and h=2*lambda
    om = ciold/(1 - ciold);
    h = 2*lambda2vec;
  } else {
    om = 0.5*m/ciold%tanh(ciold/2);
    h = lambda2vec + lambda2vec%sqrt(phi/chiold);
  }
  arma::vec hinv = 1/h; 
  arma::mat C = xr.each_row() % hinv.t();
  arma::mat Ct = C.t();
  arma::mat D = C*xrt;
  arma::mat E(n,n), F(r,n), N(nvars,n);
  
  double ldetsig;
  
  if(p <= n) {
    sigma = (x.each_col() % om).t()*x;
    arma::vec happend = h;
    happend.insert_rows(0, nvars - r, true);
    sigma.diag() += happend;
    if(elbo) {
      ldetsig = - real(log_det(sigma));
    }
    sigma = sigma.i();
    dsigma = sigma.diag();
    N = sigma*x.t();
  } else {
    if(intercept==false && unpen==false) {
      E = D;
      E.diag() += 1/om;
      E = E.i();
      F = Ct*E;
      N = Ct - F*D;
      if(posterior) {
        sigma = diagmat(hinv) - F*C;
        dsigma = sigma.diag();
      } else if(start==true){
        arma::mat NW = N.each_row() % om.t();
        dsigma = sum(NW % N,1);
        N = NW*(N.t()*x.t());
      } else {
        dsigma = hinv - sum(F%Ct,1);
      }
      if(elbo) {
        ldetsig = -real(log_det(diagmat(1/om) + D)) - sum(log(om)) - 
          sum(log(h));
      }
    } else if(intercept==true && unpen==false){ 
      double A = sum(om);
      double Ainv = 1/A;
      arma::mat Om = diagmat(om);
      Om = Om - Ainv*om*om.t();
      E = D*Om;
      E.diag() += 1;
      E = E.i();
      F = Ct*Om*E;
      arma::mat G = Ainv*om.t()*xr;
      arma::mat J = C*G.t();
      arma::mat K = G*F;
      arma::mat M = G.each_row() % hinv.t();
      N.submat(0,0,0,n-1).fill(Ainv + arma::conv_to<double>::from(M*G.t() - K*J));
      N.submat(0,0,0,n-1) += (G*F*D - G*Ct);
      N.submat(1,0,nvars-1,n-1).each_col() = F*J - M.t();
      N.submat(1,0,nvars-1,n-1) += Ct - F*D;
      if(posterior) {
        sigma.submat(0,0,0,0) = Ainv + arma::conv_to<double>::from(M*G.t() - K*J);
        sigma.submat(0,1,0,nvars-1) = K*C - M;
        sigma.submat(1,0,nvars-1,0) = F*J - M.t();
        sigma.submat(1,1,nvars-1,nvars-1) = diagmat(hinv) - F*C;
        dsigma = sigma.diag();
      } else if(start==true){
        arma::mat NW = N.each_row() % om.t();
        dsigma = sum(NW % N,1);
        N = NW*(N.t()*x.t());
      } else {
        dsigma(0) = Ainv + arma::conv_to<double>::from(M*G.t() - K*J);
        dsigma.subvec(1,nvars-1) = hinv - sum(F % Ct,1);
      }
      if(elbo) {
        arma::vec Dom = D*om;
        arma::mat OmD = diagmat(1/om) + D;
        ldetsig = -real(log_det(OmD)) - sum(log(om)) - sum(log(h)) -
          log(A - arma::conv_to<double>::from(om.t()*Dom + Dom.t()*OmD.i()*Dom));
      }
    } else { // checked (compared with calculations in R)
      arma::mat xut = xu.t();
      arma::mat Om = diagmat(om);
      arma::mat xutom = (xut.each_row() % om.t());
      arma::mat A = xutom*xu;
      arma::mat Ainv = A.i();
      Om = Om - xutom.t()*Ainv*xutom;
      arma::mat B = xutom*xr;
      E = D*Om;
      E.diag() += 1;
      E = E.i();
      F = Ct*Om*E;
      arma::mat G = Ainv*B;
      arma::mat J = G.t()*xut;
      arma::mat K = G*F;
      arma::mat L = K*C;
      arma::mat M = G.each_row() % hinv.t(); // problem
      N.submat(0,0,u-1,n-1) = Ainv*xut + M*J - L*J - G*Ct + K*D;
      N.submat(u,0,nvars-1,n-1) = F*C*J - J.each_col() % hinv + Ct - F*D;
      if(posterior) {
        arma::mat O = F*C;
        sigma.submat(0,0,u-1,u-1) = Ainv + (M - L)*G.t();
        sigma.submat(0,u,u-1,nvars-1) = -M + K*C;
        sigma.submat(u,0,nvars-1,u-1) = -M.t() + O*G.t();
        sigma.submat(u,u,nvars-1,nvars-1) = diagmat(hinv) - O;
        dsigma = sigma.diag();
      } else if(start==true){
        arma::mat NW = N.each_row() % om.t();
        dsigma = sum(NW % N,1);
        N = NW*(N.t()*x.t());
      } else {
        dsigma.subvec(0,u-1) = Ainv.diag() + sum((M - L) % G, 1);
        dsigma.subvec(u,u+r-1) = hinv - sum(F % Ct, 1);
      }
      if(elbo) {
        arma::mat BCt = B*Ct;
        arma::mat OmD = diagmat(1/om) + D;
        ldetsig = -real(log_det(OmD)) - sum(log(om)) - sum(log(h)) -
          real(log_det(A - (B.each_row() % hinv.t())*B.t() + BCt*OmD.i()*
          BCt.t()));
      }
    }
  }
  
  mu = N*kappa;
  ci = sqrt(sum(x.t() % N,0).t() + square(x*mu));
  chi = lambda2vec%(dsigma.subvec(nvars-r,nvars-1) + 
    square(mu.subvec(nvars-r,nvars-1)));
  
  if(elbo) {
    elboval = arma::conv_to<double>::from(kappa.t()*x*mu + m.t()*(0.5*ci - 
      log(exp(ci) + 1))) + 0.5*sum(log(lambdag)) - 
      0.5*sqrt(phi)*sum(sqrt(chi)) -
      0.5*lambda2*sum(lambdag%(1 + sqrt(phi/chi))%
      (dsigma.subvec(u,u+r-1) + square(mu.subvec(u,u+r-1)))) + 0.5*ldetsig;
  }
  
  if(posterior && elbo) {
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("sigma") = sigma,
                        Named("mu") = mu,
                        Named("elbo") = elboval);
  } else if(elbo){
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("dsigma") = dsigma,
                        Named("mu") = mu,
                        Named("elbo") = elboval);
  } else if(posterior) {
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("sigma") = sigma,
                        Named("mu") = mu);
  } else {
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("dsigma") = dsigma,
                        Named("mu") = mu);
  }
  
}
