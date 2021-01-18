#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Auxliary : Givens rotation matrix parameters for a and b
arma::vec rotmat(double a, double b){
  double c, s, temp;
  if (b==0.0){
    c = 1.0;
    s = 0.0;
  } else if (std::pow(b,2.0)>std::pow(a,2.0)){
    temp = a/b;
    s = 1.0/std::sqrt(static_cast<float>(1.0+std::pow(temp,2)));
    c = temp*s;
  } else {
    temp = b/a;
    c = 1.0/std::sqrt(static_cast<float>(1.0+std::pow(temp,2.0)));
    s = temp*c;
  }
  arma::vec res(2,fill::zeros);
  res(0) = c;
  res(1) = s;
  return(res);
}

//------------------------------------------------------------------------------
// Krylov Iterative Solver 7-1. GMRES
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.gmres.single)]]
Rcpp::List single_gmres(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                        const double reltol, const int maxiter, arma::mat& M, const int restrt){
  // 1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt;

  // 2. other basics
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

  // 3. initialize workspace
  int m = restrt;
  double temp = 0;

  arma::mat V(n,m+1,fill::zeros);
  arma::mat H(m+1,m,fill::zeros);
  arma::vec cs(m,fill::zeros);
  arma::vec sn(m,fill::zeros);
  arma::vec e1(n,fill::zeros);
  e1(0) = 1.0;
  arma::vec s(n,fill::zeros);
  arma::colvec w(n,fill::zeros);
  arma::vec tmpcssn(2,fill::zeros);
  arma::colvec y(m,fill::zeros);

  // 4. main iteration
  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){

    r = solve(M,b-A*x);
    V.col(0) = r/norm(r);
    s = norm(r)*e1;
    // construct orthonormal basis using Gram-Schmidt
    int i;
    for (i=0;i<m;i++){
      w = solve(M,A*V.col(i));
      for (int k=0;k<i;k++){
        H(k,i) = dot(w,V.col(k));
        w = w - H(k,i)*V.col(k);
      }
      H(i+1,i) = norm(w);
      V.col(i+1) = w/H(i+1,i);
      for (int k=0;k<(i-1);k++){ // apply givens rotation
        temp = cs(k)*H(k,i) + sn(k)*H(k+1,i);
        H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
        H(k,i)   = temp;
      }
      tmpcssn = rotmat(H(i,i),H(i+1,i));
      cs(i) = tmpcssn(0);
      sn(i) = tmpcssn(1);

      temp = cs(i)*s(i);
      s(i+1) = -sn(i)*s(i);
      s(i)   = temp;
      H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
      H(i+1,i) = 0.0;
      error  = sqrt(pow(s(i+1),2)) / bnrm2;
      if (error <= reltol){
        y = solve(H.submat(0,0,i-1,i-1), s.subvec(0,i-1));
        x = x + V.cols(0,(i-1))*y;
        break;
      }
    }


    if (error <= reltol){
      break;
    }


    y = solve(H(span(0,m-1),span(0,m-1)), s(span(0,m-1)));

    x = x + V.cols(0,(m-1))*y;

    r = solve(M, b-A*x);

    error  = norm(r) / bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error <= reltol){
      break;
    }
    iter += 1;
  }

  // 5. finishing up
  if (error<=reltol){
    flag = 0;
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
// Krylov Iterative Solver 7-2. GMRES : a sparse way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.gmres.single.sparse)]]
Rcpp::List single_gmres_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                               const double reltol, const int maxiter, arma::sp_mat M, const int restrt){
  // 1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt;

  // 2. other basics
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

  // 3. initialize workspace
  int m = restrt;
  double temp = 0;

  arma::mat V(n,m+1,fill::zeros);
  arma::mat H(m+1,m,fill::zeros);
  arma::vec cs(m,fill::zeros);
  arma::vec sn(m,fill::zeros);
  arma::vec e1(n,fill::zeros);
  e1(0) = 1.0;
  arma::vec s(n,fill::zeros);
  arma::colvec w(n,fill::zeros);
  arma::vec tmpcssn(2,fill::zeros);
  arma::colvec y(m,fill::zeros);

  // 4. main iteration

  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){

    r = spsolve(M,b-A*x,"lapack");
    V.col(0) = r/norm(r);
    s = norm(r)*e1;
    // construct orthonormal basis using Gram-Schmidt
    int i;
    for (i=0;i<m;i++){
      w = spsolve(M,A*V.col(i),"lapack");
      for (int k=0;k<i;k++){
        H(k,i) = dot(w,V.col(k));
        w = w - H(k,i)*V.col(k);
      }
      H(i+1,i) = norm(w);
      V.col(i+1) = w/H(i+1,i);
      for (int k=0;k<(i-1);k++){ // apply givens rotation
        temp = cs(k)*H(k,i) + sn(k)*H(k+1,i);
        H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
        H(k,i)   = temp;
      }
      tmpcssn = rotmat(H(i,i),H(i+1,i));
      cs(i) = tmpcssn(0);
      sn(i) = tmpcssn(1);

      temp = cs(i)*s(i);
      s(i+1) = -sn(i)*s(i);
      s(i)   = temp;
      H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
      H(i+1,i) = 0.0;
      error  = sqrt(pow(s(i+1),2)) / bnrm2;
      if (error <= reltol){
        y = solve(H.submat(0,0,i-1,i-1), s.subvec(0,i-1));
        x = x + V.cols(0,(i-1))*y;
        break;
      }
    }

    if (error <= reltol){
      break;
    }
    y = solve(H(span(0,m-1),span(0,m-1)), s(span(0,m-1)));
    x = x + V.cols(0,(m-1))*y;
    r = spsolve(M, b-A*x,"lapack");
    error  = norm(r) / bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error <= reltol){
      break;
    }
    iter += 1;
  }

  // 5. finishing up
  if (error<=reltol){
    flag = 0;
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
