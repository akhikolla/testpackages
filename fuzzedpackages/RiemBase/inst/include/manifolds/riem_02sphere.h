#ifndef RIEM_02SPHERE_H
#define RIEM_02SPHERE_H

#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

// 01. dim(x)
inline int sphere_dim(arma::mat x){
  return((x.n_rows-1));
}
// 02. inner(x,d1,d2)
inline double sphere_inner(arma::mat x, arma::mat d1, arma::mat d2){
  return(arma::dot(d1,d2));
}
// 03. norm(x,d)
inline double sphere_norm(arma::mat x, arma::mat d){
  return(arma::norm(d,"fro"));
}
// 04. dist(x,y)
inline double sphere_dist(arma::mat x, arma::mat y){
  return(acos(arma::as_scalar(x.t()*y)));
}
// 05. proj(x,u)
inline arma::mat sphere_proj(arma::mat x, arma::mat u){
  return(u-x*(arma::dot(x,u)));
}
// 06. tangent(x,u)
inline arma::mat sphere_tangent(arma::mat x, arma::mat u){
  return(sphere_proj(x,u));
}
// 07. tangent2ambient
// 08. rand(x)
inline arma::mat sphere_rand(arma::mat x){
  int n = x.n_rows;
  arma::mat out = arma::randn(n,1);
  out /= arma::norm(out,"fro");
  return(out);
}
// 09. randvec(x)
inline arma::mat sphere_randvec(arma::mat x){
  int n = x.n_rows;
  arma::mat d   = arma::randn(n,1);
  arma::mat out = d-(x*(arma::dot(x,d)));
  out /= arma::norm(out,"fro");
  return(out);
}
// 10. zerovec(x)
inline arma::mat sphere_zerovec(arma::mat x){
  int n = x.n_rows;
  arma::mat out(n,1,fill::zeros);
  return(out);
}
// 11. vec(x,u_mat)
inline arma::mat sphere_vec(arma::mat x, arma::mat u_mat){
  return(u_mat);
}
// 12. mat(x,u_vec)
inline arma::mat sphere_mat(arma::mat x, arma::mat u_vec){
  return(u_vec);
}
// 13. nearest(x)
inline arma::mat sphere_nearest(arma::mat x){
  arma::mat out = x;
  out /= arma::norm(out,"fro");
  return(out);
}
// 14. exp(x,d,t=1.0)
inline arma::mat sphere_exp(arma::mat x, arma::mat d, double t){
  double nrm_td = arma::norm(t*d, "fro"); // theta
  arma::mat out;
  if (nrm_td < 1e-15){ // very close
    out = x;
  } else {
    out = cos(nrm_td)*x + ((sin(nrm_td))/nrm_td)*t*d;
  }
  return(out);
}
// 15. log(x,y)
inline arma::mat sphere_log(arma::mat x, arma::mat y){
  arma::mat v = sphere_proj(x,y-x);
  double di = sphere_dist(x,y);
  if (di > 1e-6){
    double nv = arma::norm(v, "fro");
    v = v*(di/nv);
  }
  return(v);
}
// 16. retr(x,d,t)
inline arma::mat sphere_retr(arma::mat x, arma::mat d, double t){
  arma::mat out = x + t*d;
  return(sphere_nearest(out));
}
// 17. invretr(x,y)
inline arma::mat sphere_invretr(arma::mat x, arma::mat y){
  arma::mat d = (y/arma::dot(x,y))-x;
  return(d);
}

// 18. equiv(x,m,n)
inline arma::vec sphere_equiv(arma::mat x, int m, int n){
  arma::vec out = arma::vectorise(x,0);
  return(out);
}

// 19. invequiv(x,m,n)
inline arma::mat sphere_invequiv(arma::vec x, int m, int n){
  arma::mat out  = arma::reshape(x,m,n);
  double outsize = arma::norm(out,"fro");
  return((out/outsize));
}

// 20. extdist(x,y)
inline double sphere_extdist(arma::mat x, arma::mat y){
  int m = x.n_rows;
  int n = x.n_cols;
  
  arma::vec xext = sphere_equiv(x, m, n);
  arma::vec yext = sphere_equiv(y, m, n);
  
  return(arma::as_scalar(arma::norm(xext-yext,"fro")));
}

#endif
