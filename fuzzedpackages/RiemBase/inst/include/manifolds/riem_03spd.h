#ifndef RIEM_03SPD_H
#define RIEM_03SPD_H

#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

// auxiliary functions
inline arma::mat spdaux_symm(arma::mat x){
  return((x+x.t())/2.0);
}

// 01. dim(x)
inline int spd_dim(arma::mat x){
  const int n = x.n_rows;
  return((n*(n+1))/2);
}
// 02. inner(x,d1,d2)
inline double spd_inner(arma::mat x, arma::mat u, arma::mat v){
  arma::mat sol1 = arma::solve(x,u);
  arma::mat sol2 = arma::solve(x,v);
  // use of trinner
  return(arma::as_scalar(arma::trace(sol1.t()*sol2)));
}
// 03. norm(x,d)
inline double spd_norm(arma::mat x, arma::mat eta){
  arma::mat sol  = arma::solve(x,eta);
  arma::vec rvec = arma::diagvec(sol*sol);

  return(std::sqrt(arma::as_scalar(arma::sum(rvec))));
}
// 04. dist(x,y)
inline double spd_dist(arma::mat x, arma::mat y){
  arma::mat sol = arma::solve(x,y);
  arma::cx_mat cxXY = arma::logmat(sol);
  arma::mat XY = arma::real(cxXY);
  return(std::sqrt(arma::as_scalar(arma::trace(XY*XY))));
}
// 05. proj(x,u)
inline arma::mat spd_proj(arma::mat x, arma::mat u){
  return(spdaux_symm(u));
}
// 06. tangent(x,u)
inline arma::mat spd_tangent(arma::mat x, arma::mat u){
  return(spd_proj(x,u));
}
// 07. tangent2ambient
inline arma::mat spd_tangent2ambient(arma::mat x, arma::mat eta){
  return(eta);
}
// 08. rand(x)
inline arma::mat spd_rand(arma::mat x){
  int n = x.n_cols;
  arma::vec d(n,fill::randu);
  arma::mat D = arma::diagmat(1.0+d);

  arma::mat normQ(n,n,fill::randn);
  arma::mat Q,R;
  arma::qr(Q,R,normQ);

  return(Q*D*Q.t());
}
// 09. randvec(x)
inline arma::mat spd_randvec(arma::mat x){
  int n = x.n_cols;
  arma::mat tmpR(n,n,fill::randn);
  tmpR = spd_proj(x, tmpR); // Project
  tmpR /= spd_norm(x,tmpR);
  return(tmpR);
}
// 10. zerovec(x)
inline arma::mat spd_zerovec(arma::mat x){
  int n = x.n_rows;
  arma::mat out(n,n,fill::zeros);
  return(out);
}
// 11. vec(x,u_mat)
inline arma::mat spd_vec(arma::mat x, arma::mat u_mat){
  int m = x.n_rows;
  int n = x.n_cols;
  arma::mat out = arma::reshape(u_mat,m*n,1);
  return(out);
}
// 12. mat(x,u_vec)
inline arma::mat spd_mat(arma::mat x, arma::mat u_vec){
  arma::mat out = arma::reshape(u_vec, x.n_rows, x.n_cols);
  return(out);
}
// 13. nearest(x)
inline arma::mat spd_nearest(arma::mat x){
  const int n = x.n_rows;
  arma::mat y = (x+x.t())/2.0;

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval,eigvec,y);

  for (int i=0;i<n;i++){
    if (eigval(i) <= 0.0){
      eigval(i) = 1e-12;
    }
  }

  arma::mat output = eigvec*(arma::diagmat(eigval))*eigvec.t();
  return(output);
}
// 14. exp(x,d,t)
inline arma::mat spd_exp(arma::mat x, arma::mat eta, double t){
  arma::mat teta = t*eta;
  arma::mat tmp  = arma::expmat(arma::solve(x, teta));
  arma::mat yy   = x*tmp;
  return((yy+yy.t())/2.0);
}
// 15. log(x,y)
inline arma::mat spd_log(arma::mat x, arma::mat y){
  arma::mat tmp  = arma::real(arma::logmat(arma::solve(x,y)));
  arma::mat yy   = x*tmp;
  return((yy+yy.t())/2.0);
}
// 16. retr(x,d,t)
inline arma::mat spd_retr(arma::mat x, arma::mat eta, double t){
  arma::mat teta = t*eta;
  return(spdaux_symm(x + teta + ((teta/2.0)*(arma::solve(x, teta)))));
}
// // 17. invretr(x,y)



// 18. equiv(x,m,n)
inline arma::vec spd_equiv(arma::mat x, int m, int n){
  arma::cx_mat logx = arma::logmat(x);
  arma::mat realogx = arma::real(logx);
  arma::vec out = arma::vectorise(realogx,0);
  return(out);
}

// 19. invequiv(x,m,n)
inline arma::mat spd_invequiv(arma::vec x, int m, int n){
  arma::mat tmpx = arma::reshape(x,m,n);
  arma::mat output = arma::expmat(((tmpx + tmpx.t())/2.0));
  return(output);
}

// 20. extdist(x,y)
inline double spd_extdist(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = spd_equiv(x, m, n);
  arma::vec yext = spd_equiv(y, m, n);
  
  return(arma::as_scalar(arma::norm(xext-yext,"fro")));
}

#endif
