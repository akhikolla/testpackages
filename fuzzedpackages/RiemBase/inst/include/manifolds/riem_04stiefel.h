#ifndef RIEM_04STIEFEL_H
#define RIEM_04STIEFEL_H

#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

// 01. dim(x)
inline int stiefel_dim(arma::mat x){
  int n = x.n_rows;
  int p = x.n_cols;
  int output = ((n*p)-(p*(p+1)/2));
  return(output);
}
// 02. inner(x,d1,d2)
inline double stiefel_inner(arma::mat x, arma::mat d1, arma::mat d2){
  return(arma::as_scalar(arma::dot(arma::vectorise(d1), arma::vectorise(d2))));
}
// 03. norm(x,d)
inline double stiefel_norm(arma::mat x, arma::mat d){
  return(arma::norm(d, "fro"));
}

// 05. proj(x,u)
inline arma::mat stiefel_proj(arma::mat x, arma::mat u){
  arma::mat A = (x.t()*u);
  return(u - x*((A+A.t())/2.0));
}
// 06. tangent(x,u)
inline arma::mat stiefel_tangent(arma::mat x, arma::mat u){
  return(stiefel_proj(x,u));
}
// 07. tangent2ambient
// 08. rand(x)
inline arma::mat stiefel_rand(arma::mat x){
  int n = x.n_rows;
  int p = x.n_cols;
  
  arma::mat A(n,p,fill::randn);
  arma::mat Q,R;
  arma::qr(Q,R,A);

  return(Q);
}
// 09. randvec(x)
inline arma::mat stiefel_randvec(arma::mat x){
  int n = x.n_rows;
  int p = x.n_cols;
  
  arma::mat A(n,p,fill::randn);
  arma::mat U;
  U = stiefel_proj(x,A);
  U /= arma::norm(U,"fro");
  return(U);
}
// 10. zerovec(x)
inline arma::mat stiefel_zerovec(arma::mat x){
  int n = x.n_rows;
  int p = x.n_cols;
  
  arma::mat out(n,p,fill::zeros);
  return(out);
}
// 11. vec(x,u_mat)
inline arma::mat stiefel_vec(arma::mat x, arma::mat u_mat){
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat out = arma::reshape(u_mat,(n*p),1);
  return(out);
}
// 12. mat(x,u_vec)
inline arma::mat stiefel_mat(arma::mat x, arma::mat u_vec){
  arma::mat out = arma::reshape(u_vec, x.n_rows, x.n_cols);
  return(out);
}
// 13. nearest(x) {heuristic}
inline arma::mat stiefel_nearest(arma::mat x){
  arma::mat Q,R;
  arma::qr_econ(Q,R,x);
  return(Q);
}
// 14. exp(x,d,t)
inline arma::mat stiefel_exp(arma::mat x, arma::mat u, double t){
  const int n = x.n_rows;
  const int p = x.n_cols;
  
  arma::mat Ip(p,p,fill::eye);
  arma::mat Zp(p,p,fill::zeros);
  
  arma::mat tu = t*u;
  arma::mat term1 = arma::join_horiz(x, tu);
  
  arma::mat term21 = arma::join_horiz((x.t()*tu), -((tu.t())*tu));
  arma::mat term22 = arma::join_horiz(Ip, (x.t()*tu));
  arma::mat term2  = arma::expmat(arma::join_vert(term21, term22));
  
  arma::mat term3  = arma::join_vert(arma::expmat(-(x.t()*tu)), Zp);
  
  arma::mat output = term1*term2*term3;
  return(output);
}
// 15. log(x,y) + // 04. dist(x,y)
inline arma::mat stiefel_log(arma::mat U0, arma::mat U1){
  const int n = U0.n_rows;
  const int p = U0.n_cols;
  const double tau = 1e-6;   // default convergence threshold
  const int maxiter = 12345; // default maximum number of iterations
  
  // 1.
  arma::mat M = U0.t()*U1;   
  // 2. thin QR of normal component of U1
  arma::mat Q,N;
  arma::qr_econ(Q,N,U1-(U0*M));
  // 3. orthogonal completion + procrustes preprocessing ------------- no QR_ECON ?
  arma::mat V, Vaway;
  arma::qr(V, Vaway, arma::join_vert(M,N)); 
  
  arma::mat D, R; 
  arma::vec vecS;
  arma::svd(D,vecS,R,V.submat(p,p,(2*p)-1,(2*p)-1));
  arma::mat S = arma::diagmat(vecS);
  V.cols(p,(2*p)-1) = V.cols(p,(2*p)-1)*(R*D.t());
  V = arma::join_horiz(arma::join_vert(M,N),V.cols(p,(2*p)-1)); 
  
  // 4. for looping
  arma::cx_mat LVcx;
  arma::mat LV, C, Phi;
  double normC;
  for (int k=0;k<maxiter;k++){
    LVcx = arma::logmat(V);
    LV   = arma::real(LVcx);
    
    // lower (pxp) diagonal block
    C = LV.submat(p,p,(2*p)-1,(2*p)-1);
    // convergence check
    normC = arma::norm(C, 2);
    if (normC < tau){
      break;
    }
    // matrix exponential
    Phi = arma::expmat(-C);
    // update last p columns
    V.cols(p,(2*p)-1) = V.cols(p,(2*p)-1)*Phi;
  }
  
  // 5. prepare output
  arma::mat Delta = (U0*LV.submat(0,0,(p-1),(p-1))) + (Q*LV.submat(p,0,(2*p)-1,p-1));
  return(Delta);
}
inline double stiefel_dist(arma::mat x, arma::mat y){
  arma::mat delta = stiefel_log(x,y);
  double output = stiefel_norm(x, delta);
  return(output);
}
// 16. retr(x,d,t)
inline arma::mat stiefel_retr(arma::mat x, arma::mat u, double t){
  arma::mat y = x + (t*u);
  arma::mat Q,R;
  arma::qr(Q,R,y);
  return(Q);
}
// 17. invretr(x,y)
inline arma::mat stiefel_invretr(arma::mat x, arma::mat y){
  int n = x.n_rows;
  int p = x.n_cols;
  
  arma::mat xty = (x.t()*y);
  arma::mat Ip(p,p,fill::eye);
  
  arma::mat A = xty;
  arma::mat B = xty.t();
  arma::mat C = -2.0*Ip;
  
  arma::mat M = arma::syl(A,B,C);
  arma::mat U = y*M-x;
  return(U);
}

// 18. equiv(x,m,n)
inline arma::vec stiefel_equiv(arma::mat x, int m, int n){
  arma::vec output = arma::vectorise(x,0);
  return(output);
}

// 19. invequiv(x,m,n)
inline arma::mat stiefel_invequiv(arma::vec x, int m, int n){
  arma::mat mu = arma::reshape(x,m,n);
  arma::mat rhs = arma::pinv(arma::real(arma::sqrtmat(mu.t()*mu)));
  arma::mat output = mu*rhs;
  return(output);
}

// 20. extdist(x,y)
inline double stiefel_extdist(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = stiefel_equiv(x, m, n);
  arma::vec yext = stiefel_equiv(y, m, n);
  
  return(arma::as_scalar(arma::norm(xext-yext,"fro")));
}

#endif
