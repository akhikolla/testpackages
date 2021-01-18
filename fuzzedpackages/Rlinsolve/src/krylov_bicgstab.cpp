#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Krylov Iterative Solver 3-1. Biconjugate Gradient Stabilized method
// From NETLIB : http://www.netlib.org/templates/matlab//bicgstab.m
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.bicgstab.single)]]
Rcpp::List single_bicgstab(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                           const double reltol, const int maxiter, const arma::mat& M){
  // 1-1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt; -1 breakdown rho=0; -2 breakdown omega=0;

  // 1-2. Preiteration
  double bnrm2 = norm(b);
  if (bnrm2==0){
    bnrm2=1.0;
  }

  arma::colvec x = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    x = xtmp;
  }
  arma::colvec r = b-A*x;
  double error = norm(r)/bnrm2;

  double omega = 1.0;
  arma::colvec r_tld = r;

  // 1-3. main iteration
  // arma::mat Minv = pinv(M);
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;
  double resid = 0;
  
  arma::colvec p(n,fill::zeros);
  arma::colvec p_hat(n,fill::zeros);
  arma::colvec v(n,fill::zeros);
  arma::colvec s(n,fill::zeros);
  arma::colvec s_hat(n,fill::zeros);
  arma::colvec t(n,fill::zeros);

  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;

  while (iter<maxiter){
    rho = dot(r_tld,r); // direction vector
    if (rho==0.0){
      break;
    }
    if (iter>0){
      beta = (rho/rho_1)*(alpha/omega);
      p    = r+beta*(p-omega*v);
    } else {
      p    = r;
    }

    p_hat = solve(M,p);
    v     = A*p_hat;
    alpha = rho/(dot(r_tld,v));
    s     = r - alpha*v;
    if (norm(s)<reltol){ // early convergence check
      x = x+alpha*p_hat;
      resid = norm(s)/bnrm2;
      break;
    }

    s_hat = solve(M,s); // stabilizer
    t     = A*s_hat;
    omega = dot(t,s)/dot(t,t);

    x = x+alpha*p_hat+omega*s_hat;
    r = s-omega*t;
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }

    if (error<reltol){
      break;
    }
    if (omega==0.0){
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  // 1-4. postprocess
  if ((error<=reltol)||all(s<=reltol)){
    if (all(s<=reltol)){
      flag = 0;
    }
  } else if (omega==0.0){
    flag = -2;
  } else if (rho==0.0){
    flag = -1;
  } else {
    flag = 1;
  }

  // 1-4. return outputs
  List res;
  res["x"] = x;
  res["iter"] = iter;
  if (iter>=maxiter){
    res["errors"] = errors;
  } else {
    arma::vec newerrors = errors.subvec(0,iter);
    res["errors"] = newerrors;
  }
  res["flag"] = flag;
  return res;
}

//------------------------------------------------------------------------------
// Krylov Iterative Solver 3-2. Biconjugate Gradient Stabilized method : sprase way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.bicgstab.single.sparse)]]
Rcpp::List single_bicgstab_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                                  const double reltol, const int maxiter, const arma::sp_mat M){
  // 1-1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt; -1 breakdown rho=0; -2 breakdown omega=0;
  // 1-2. Preiteration
  double bnrm2 = norm(b);
  if (bnrm2==0){
    bnrm2=1.0;
  }

  arma::colvec x = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    x = xtmp;
  }
  arma::colvec r = b-A*x;
  double error = norm(r)/bnrm2;

  double omega = 1.0;
  arma::colvec r_tld = r;

  // 1-3. main iteration
  // arma::mat Minv = pinv(M);
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;
  double resid = 0;
  
  arma::colvec p(n,fill::zeros);
  arma::colvec p_hat(n,fill::zeros);
  arma::colvec v(n,fill::zeros);
  arma::colvec s(n,fill::zeros);
  arma::colvec s_hat(n,fill::zeros);
  arma::colvec t(n,fill::zeros);

  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;

  while (iter<maxiter){
    rho = dot(r_tld,r); // direction vector
    if (rho==0.0){
      break;
    }
    if (iter>0){
      beta = (rho/rho_1)*(alpha/omega);
      p    = r+beta*(p-omega*v);
    } else {
      p    = r;
    }

    p_hat = spsolve(M,p,"lapack");
    v     = A*p_hat;
    alpha = rho/(dot(r_tld,v));
    s     = r - alpha*v;
    if (norm(s)<reltol){ // early convergence check
      x = x+alpha*p_hat;
      resid = norm(s)/bnrm2;
      break;
    }

    s_hat = spsolve(M,s,"lapack"); // stabilizer
    t     = A*s_hat;
    omega = dot(t,s)/dot(t,t);

    x = x+alpha*p_hat+omega*s_hat;
    r = s-omega*t;
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }

    if (error<reltol){
      break;
    }
    if (omega==0.0){
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  // 1-4. postprocess
  if ((error<=reltol)||all(s<=reltol)){
    if (all(s<=reltol)){
      flag = 0;
    }
  } else if (omega==0.0){
    flag = -2;
  } else if (rho==0.0){
    flag = -1;
  } else {
    flag = 1;
  }

  // 1-4. return outputs
  List res;
  res["x"] = x;
  res["iter"] = iter;
  if (iter>=maxiter){
    res["errors"] = errors;
  } else {
    arma::vec newerrors = errors.subvec(0,iter);
    res["errors"] = newerrors;
  }
  res["flag"] = flag;
  return res;
}

