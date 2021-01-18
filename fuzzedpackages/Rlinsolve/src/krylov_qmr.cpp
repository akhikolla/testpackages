#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Krylov Iterative Solver 5. QMR
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.qmr.single)]]
Rcpp::List single_qmr(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                      const double reltol, const int maxiter, const arma::mat& M){
  // 1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt; -1 breakdown

  // 2. Preiteration
  double bnrm2 = norm(b);
  if (bnrm2==0){
    bnrm2=1.0;
  }

  // 3. error preparation
  arma::colvec x = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    x = xtmp;
  }
  arma::colvec r = b-A*xinit;
  arma::vec errors(maxiter,fill::zeros);
  double error = norm(r)/bnrm2;
  errors(0) = error;

  // 4. LU decomposition
  arma::mat M1(n,n,fill::zeros);
  arma::mat M2(n,n,fill::zeros);
  lu(M1,M2,M);

  // 5. others I do not get related yet.
  arma::colvec v_tld = r;
  arma::colvec y = solve(M1,v_tld);
  double rho = norm(y);

  arma::colvec w_tld = r;
  arma::colvec z = solve(M2.t(),w_tld);
  double xi = norm(z);

  double gamma = 1.0;
  double eta   = -1.0;
  double theta = 0.0;

  // 6. Main Iteration
  arma::colvec p(n,fill::zeros);
  arma::colvec q(n,fill::zeros);
  arma::colvec v(n,fill::zeros);
  arma::colvec w(n,fill::zeros);
  arma::colvec y_tld(n,fill::zeros);
  arma::colvec z_tld(n,fill::zeros);
  arma::colvec p_tld(n,fill::zeros);
  arma::colvec d(n,fill::zeros);
  arma::colvec s(n,fill::zeros);
  double delta = 0.0;
  double ep    = 0.0;
  double beta  = 0.0;
  double rho_1 = 0.0;
  double gamma_1=0.0;
  double theta_1=0.0;
  double beta_1 =0.0;
  while (iter<maxiter){
    if ((rho==0.0)|(xi==0.0)){
      break;
    }

    v = v_tld/rho;
    y = y/rho;

    w = w_tld/xi;
    z = z/xi;

    delta = dot(z,y);
    if (delta==0.0){
      break;
    }

    y_tld = solve(M2,y);
    z_tld = solve(M1.t(),z);

    if (iter>0){ // direction vectors
      p = y_tld - (xi*delta/ep)*p;
      q = z_tld - (rho*delta/ep)*q;
    } else {
      p = y_tld;
      q = z_tld;
    }

    p_tld = A*p;
    ep = dot(q,p_tld);
    if (ep==0.0){
      break;
    }

    beta = ep/delta;
    if (beta==0.0){
      break;
    }

    v_tld = p_tld - beta*v;
    y = solve(M1, v_tld);

    rho_1 = rho;
    rho = norm(y);
    w_tld = (A.t()*q) - (beta*w);
    z = solve(M2.t(), w_tld);

    xi = norm(z);

    gamma_1 = gamma;
    theta_1 = theta;

    theta = rho/(gamma_1*beta);
    gamma = 1.0/sqrt(1.0+(theta*theta));
    if (gamma==0.0){
      break;
    }

    eta = -eta*rho_1*(gamma*gamma)/(beta*(gamma_1*gamma_1));
    if (iter>0){                      // compute adjustment
      d = eta*p + ((pow(theta_1*gamma,2))*d);
      s = eta*p_tld + ((pow(theta_1*gamma,2))*s);
    } else {
      d = eta*p;
      s = eta*p_tld;
    }

    x = x + d; // update approximation
    r = r - s; // update residual
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error<reltol){
      break;
    }
    iter += 1;
  }

  // 7. flag information
  if (error <= reltol){
    flag = 0;
  } else if (rho==0.0){
    flag = -1;
  } else if (beta==0.0){
    flag = -2;
  } else if (gamma==0.0){
    flag = -3;
  } else if (delta==0.0){
    flag = -4;
  } else if (ep==0.0){
    flag = -5;
  } else if (xi==0.0){
    flag = -6;
  } else {
    flag = 1;
  }

  // 8. return outputs
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
// Krylov Iterative Solver 5. QMR : Spasre Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.qmr.single.sparse)]]
Rcpp::List single_qmr_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                             const double reltol, const int maxiter, const arma::sp_mat M,
                             const arma::sp_mat M1, const arma::sp_mat M2){
  // 1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt; -1 breakdown

  // 2. Preiteration
  double bnrm2 = norm(b);
  if (bnrm2==0){
    bnrm2=1.0;
  }

  // 3. error preparation
  arma::colvec x = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    x = xtmp;
  }
  arma::colvec r = b-A*xinit;
  arma::vec errors(maxiter,fill::zeros);
  double error = norm(r)/bnrm2;
  errors(0) = error;

  // 4. LU decomposition
  // M1 and M2 from sparse decompositions
  //arma::mat M1(n,n,fill::zeros);
  //arma::mat M2(n,n,fill::zeros);

  // 5. others I do not get related yet.
  arma::colvec v_tld = r;
  arma::colvec y = spsolve(M1,v_tld,"lapack");
  double rho = norm(y);

  arma::colvec w_tld = r;
  arma::colvec z = spsolve(M2.t(),w_tld,"lapack");
  double xi = norm(z);

  double gamma = 1.0;
  double eta   = -1.0;
  double theta = 0.0;

  // 6. Main Iteration
  arma::colvec p(n,fill::zeros);
  arma::colvec q(n,fill::zeros);
  arma::colvec v(n,fill::zeros);
  arma::colvec w(n,fill::zeros);
  arma::colvec y_tld(n,fill::zeros);
  arma::colvec z_tld(n,fill::zeros);
  arma::colvec p_tld(n,fill::zeros);
  arma::colvec d(n,fill::zeros);
  arma::colvec s(n,fill::zeros);
  double delta = 0.0;
  double ep    = 0.0;
  double beta  = 0.0;
  double rho_1 = 0.0;
  double gamma_1=0.0;
  double theta_1=0.0;
  double beta_1 =0.0;
  while (iter<maxiter){
    if ((rho==0.0)|(xi==0.0)){
      break;
    }

    v = v_tld/rho;
    y = y/rho;

    w = w_tld/xi;
    z = z/xi;

    delta = dot(z,y);
    if (delta==0.0){
      break;
    }

    y_tld = spsolve(M2,y,"lapack");
    z_tld = spsolve(M1.t(),z,"lapack");

    if (iter>0){ // direction vectors
      p = y_tld - (xi*delta/ep)*p;
      q = z_tld - (rho*delta/ep)*q;
    } else {
      p = y_tld;
      q = z_tld;
    }

    p_tld = A*p;
    ep = dot(q,p_tld);
    if (ep==0.0){
      break;
    }

    beta = ep/delta;
    if (beta==0.0){
      break;
    }

    v_tld = p_tld - beta*v;
    y = spsolve(M1, v_tld,"lapack");

    rho_1 = rho;
    rho = norm(y);
    w_tld = (A.t()*q) - (beta*w);
    z = spsolve(M2.t(), w_tld,"lapack");

    xi = norm(z);

    gamma_1 = gamma;
    theta_1 = theta;

    theta = rho/(gamma_1*beta);
    gamma = 1.0/sqrt(1.0+(theta*theta));
    if (gamma==0.0){
      break;
    }

    eta = -eta*rho_1*(gamma*gamma)/(beta*(gamma_1*gamma_1));
    if (iter>0){                      // compute adjustment
      d = eta*p + ((pow(theta_1*gamma,2))*d);
      s = eta*p_tld + ((pow(theta_1*gamma,2))*s);
    } else {
      d = eta*p;
      s = eta*p_tld;
    }

    x = x + d; // update approximation
    r = r - s; // update residual
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error<reltol){
      break;
    }
    iter += 1;
  }

  // 7. flag information
  if (error <= reltol){
    flag = 0;
  } else if (rho==0.0){
    flag = -1;
  } else if (beta==0.0){
    flag = -2;
  } else if (gamma==0.0){
    flag = -3;
  } else if (delta==0.0){
    flag = -4;
  } else if (ep==0.0){
    flag = -5;
  } else if (xi==0.0){
    flag = -6;
  } else {
    flag = 1;
  }

  // 8. return outputs
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
