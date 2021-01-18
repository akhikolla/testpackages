#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

vec F(vec x, int N, mat Rx) {
  
  mat Mx(N,N,fill::zeros);
  mat Mx2(N,N,fill::zeros);
  vec M(N*N,fill::zeros);
  
  Mx = reshape(x,N,N);
  Mx = Mx.t();
  
  Mx2 = sum(Mx,0);
  Mx2 = repmat(Mx2,N,1);
  
  Mx.for_each( [](mat::elem_type& val) { val = R::digamma(val); } );
  Mx2.for_each( [](mat::elem_type& val) { val = R::digamma(val); } );
  Mx = Mx - Mx2;
  
  Mx = Mx - Rx;
  
  Mx = Mx.t();
  M = vectorise(Mx);
  
  return M;
  
}


#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

mat J(vec x, int m, mat Rxj) {
  
  int N = x.n_elem;
  mat jac(m,N,fill::zeros);
  vec hh(N,fill::zeros);
  int i;
  
  for (i=0;i!=N;i++)
    {
      hh[i] = 0.000006;
      jac.col(i) = (F(x + hh,Rxj.n_cols,Rxj) - F(x - hh,Rxj.n_cols,Rxj)) / 0.000012;
      hh[i] = 0;
    }
   
  return jac;
  
}


#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

vec newt(vec x0, int Nn, mat Rxn) {
  
  int niter = 0;
  double err = 0.000000015 + 1;
  int mm = Nn * Nn;
  vec Ft(mm,fill::zeros);
  vec delta(mm,fill::zeros);
  mat jacob(mm,Nn,fill::zeros);
  
  while (err >= 0.000000015 && niter < 10) {
    
    niter = niter + 1;
    Ft = F(x0,Nn,Rxn);
    jacob = J(x0,mm,Rxn);
    delta = -1 * solve(jacob,Ft);
    x0 = x0 + delta;
    err = sqrt(accu(delta % delta));
  
  }

  return x0;
  
}


#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

mat mKIapprox(mat w, mat vold) {

  mat R(vold.n_rows,vold.n_cols,fill::zeros);
  int Nr = vold.n_rows; 
  int Nc = vold.n_cols;
  int i;
  int j;
  mat vR(vold.n_rows,vold.n_cols,fill::zeros);
  double v;
  mat M(vold.n_rows,vold.n_cols,fill::eye);
  mat vRr(vold.n_rows,vold.n_cols,fill::zeros);
  mat vnew(vold.n_rows,vold.n_cols,fill::zeros);
  vec vnew2(Nr * Nc,fill::zeros);
  
  for (i=0;i!=Nr;i++)
    {
      for (j=0;j!=Nc;j++)
        {
          v = vold(i,j);
          vR.fill(v);
          vR = vR + M;
        
          vRr = sum(vR,0);
          vRr = repmat(vRr,Nc,1);
          vR.for_each( [](mat::elem_type& val) { val = R::digamma(val); } );
          vRr.for_each( [](mat::elem_type& val) { val = R::digamma(val); } );
          vR = vR - vRr;
          
          v = accu(w % vR);
          R(i,j) = v;
        }
    }
  
  vold = vold.t();
  vnew2 = newt(vectorise(vold),Nr,R);
  vnew = reshape(vnew2,Nr,Nc);
  vnew = vnew.t();

  return vnew;

}

