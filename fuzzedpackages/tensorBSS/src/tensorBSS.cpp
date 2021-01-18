#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP rnormCube(SEXP varp1, SEXP varp2, SEXP varp3) {
  
  int p1 = as<int>(varp1);
  int p2 = as<int>(varp2);
  int p3 = as<int>(varp3);
  
  cube ncube(p1, p2, p3, fill::randn);

  
  return Rcpp::List::create(Rcpp::Named("ncube") = ncube//,
                            //Rcpp::Named("varp1") = p1,
                            //Rcpp::Named("varp2") = p2,
                            //Rcpp::Named("varp3") = p3
                            );
  
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP eigenVectors(SEXP varx) {
  
  mat x = as<mat>(varx);
  
  mat eig_vec;
  vec eig_val;
  
  eig_sym(eig_val, eig_vec, x);
  
  return Rcpp::wrap(fliplr(eig_vec));
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP symmetricPower(SEXP varx, SEXP varr) {
  
  mat x = as<mat>(varx);
  float r = as<float>(varr);
  
  mat eig_vec;
  vec eig_val;
  
  eig_sym(eig_val, eig_vec, x);
  
  mat pow_val = diagmat(pow(sort(eig_val, "descend"), r));
  
  return Rcpp::wrap(fliplr(eig_vec)*pow_val*fliplr(eig_vec).t());
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP mFOBIMatrix(SEXP varx) {
  
  cube xcube = as<cube>(varx);
  int rows = xcube.n_rows;
  int cols = xcube.n_cols;
  int slices = xcube.n_slices;
  
  mat matFOBI(rows, rows, fill::zeros);
  
  for (int i = 0; i < slices; i++)
  {
    matFOBI = matFOBI + xcube.slice(i)*xcube.slice(i).t()*xcube.slice(i)*xcube.slice(i).t();
  }
  
  matFOBI = matFOBI/(cols*slices);
  
  return Rcpp::wrap(matFOBI);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP mFOBIMatrixNorm(SEXP varx) {
  
  cube xcube = as<cube>(varx);
  int rows = xcube.n_rows;
  int cols = xcube.n_cols;
  int slices = xcube.n_slices;
  
  mat matFOBI(rows, rows, fill::zeros);
  
  for (int i = 0; i < slices; i++)
  {
    matFOBI = matFOBI + pow(norm(xcube.slice(i), "fro"), 2)*xcube.slice(i)*xcube.slice(i).t();
  }
  
  matFOBI = matFOBI/(cols*slices);
  
  return Rcpp::wrap(matFOBI);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP mJADEMatrix(SEXP varx, SEXP vari, SEXP varj, SEXP varcov) {
  
  cube xcube = as<cube>(varx);
  int rows = xcube.n_rows;
  int cols = xcube.n_cols;
  int slices = xcube.n_slices;
  
  float i = as<float>(vari) - 1;
  float j = as<float>(varj) - 1;
  
  int delta = (i == j);
  
  mat cov = as<mat>(varcov);
  
  mat matJADE(rows, rows, fill::zeros);
  
  mat matEij(rows, rows, fill::zeros);
  matEij(i, j) = 1;
  
  arma::mat matEji(rows, rows, fill::zeros);
  matEji(j, i) = 1;
  
  arma::mat matEye(rows, rows);
  matEye.eye();
  
  for (int t = 0; t < slices; t++)
  {
    matJADE = matJADE + dot(xcube.slice(t).row(i), xcube.slice(t).row(j).t())*xcube.slice(t)*xcube.slice(t).t();
  }
  matJADE = matJADE/(cols*slices);
  
  matJADE = matJADE - cov*(matEij + matEji + delta*cols*matEye)*cov;
  
  return Rcpp::wrap(matJADE);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP matrixCovariance(SEXP varx) {
  
  cube xcube = as<cube>(varx);
  int rows = xcube.n_rows;
  int cols = xcube.n_cols;
  int slices = xcube.n_slices;
  
  mat matCov(rows, rows, fill::zeros);
  
  for (int i = 0; i < slices; i++)
  {
    matCov = matCov + xcube.slice(i)*xcube.slice(i).t();
  }
  
  matCov = matCov/(cols*slices);
  
  return Rcpp::wrap(matCov);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP mAutoCovMatrix(SEXP varx, SEXP varlag) {
  
  cube xcube = as<cube>(varx);
  int rows = xcube.n_rows;
  int cols = xcube.n_cols;
  int slices = xcube.n_slices;
    
  int lag = as<int>(varlag);
	
  mat matAutoCov(rows, rows, fill::zeros);
  
  for (int t = 0; t < (slices - lag); t++)
  {
    matAutoCov = matAutoCov + xcube.slice(t)*xcube.slice(t + lag).t();
  }
  matAutoCov = matAutoCov/(cols*(slices - lag));
   
  return Rcpp::wrap(matAutoCov);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP mTGFOBIMatrix(SEXP varx, SEXP varlag) {
  
  cube xcube = as<cube>(varx);
  int rows = xcube.n_rows;
  int cols = xcube.n_cols;
  int slices = xcube.n_slices;
    
  int lag = as<int>(varlag);
	
  mat matTGFOBI(rows, rows, fill::zeros);
  
  for (int t = 0; t < (slices - lag); t++)
  {
    matTGFOBI = matTGFOBI + xcube.slice(t)*xcube.slice(t + lag).t()*xcube.slice(t + lag)*xcube.slice(t).t();
  }
  
  matTGFOBI = matTGFOBI/(cols*(slices - lag));
  
  return Rcpp::wrap(matTGFOBI);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP mTGJADEMatrix(SEXP varx, SEXP vari, SEXP varj, SEXP varlags) {
  
  cube xcube = as<cube>(varx);
  int rows = xcube.n_rows;
  int cols = xcube.n_cols;
  int slices = xcube.n_slices;
  
  float i = as<float>(vari) - 1;
  float j = as<float>(varj) - 1;
 
  vec lags = as<vec>(varlags);
  
  mat matJADE(rows, rows, fill::zeros);
  
  for (int t = 0; t < (slices - lags(1)); t++)
  {
    matJADE = matJADE + dot(xcube.slice(t + lags(0)).row(i), xcube.slice(t + lags(1)).row(j).t())*xcube.slice(t + lags(2))*xcube.slice(t + lags(3)).t();
  }
  matJADE = matJADE/(cols*(slices - lags(1)));
  
  return Rcpp::wrap(matJADE);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP computeh(SEXP varuk, SEXP varxm, SEXP varnl) {
  
  cube xm = as<cube>(varxm);
  int n = xm.n_slices;
  
  int nl = as<int>(varnl);
  
  vec uk = as<vec>(varuk);
  
  double h = 0;
  
  // Find X^t u and add its squared norm to the sum
  if(nl == 1){
    for (int i = 0; i < n; i++)
    {
      h = h + pow(norm((xm.slice(i).t())*uk, 2), 4);
    }
  }
  if(nl == 2){
    for (int i = 0; i < n; i++)
    {
      h = h + pow(norm((xm.slice(i).t())*uk, 2), 3);
    }
  }
  if(nl == 3){
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
      temp = norm((xm.slice(i).t())*uk, 2);
      h = h + log(cosh(temp));
    }
  }
  
  h = h/n;
  
  return Rcpp::wrap(h);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP computeT(SEXP varuk, SEXP varxm, SEXP varnl) {
  
  cube xm = as<cube>(varxm);
  int p = xm.n_rows;
  int q = xm.n_cols;
  int n = xm.n_slices;
  
  int nl = as<int>(varnl);
  
  vec uk = as<vec>(varuk);
  
  vec Tcol(p, fill::zeros);
  
  vec temp(q, fill::zeros);
  
  // Find X^t u and add its squared norm to the sum
  if(nl == 1){
    for (int i = 0; i < n; i++)
    {
      temp = (xm.slice(i).t())*uk;
      Tcol = Tcol + 2*pow(norm(temp, 2), 2)*xm.slice(i)*temp;
    }
  }
  if(nl == 2){
    for (int i = 0; i < n; i++)
    {
      temp = (xm.slice(i).t())*uk;
      Tcol = Tcol + 1.5*norm(temp, 2)*xm.slice(i)*temp;
    }
  }
  if(nl == 3){
    double temp2 = 0;
    for (int i = 0; i < n; i++)
    {
      temp = (xm.slice(i).t())*uk;
      temp2 = norm(temp, 2);
      Tcol = Tcol + 0.5*(tanh(temp2)/temp2)*xm.slice(i)*temp;
    }
  }
  
  Tcol = Tcol/n;
  
  return Rcpp::wrap(Tcol);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP computed(SEXP varuk, SEXP varxm, SEXP varnl) {
  
  cube xm = as<cube>(varxm);
  int n = xm.n_slices;
  
  int nl = as<int>(varnl);
  
  vec uk = as<vec>(varuk);
  
  double d = 0;
  
  // Find X^t u and add its squared norm to the sum
  if(nl == 1){
    for (int i = 0; i < n; i++)
    {
      d = d + 2*pow(norm((xm.slice(i).t())*uk, 2), 2);
    }
  }
  if(nl == 2){
    for (int i = 0; i < n; i++)
    {
      d = d + 1.5*norm((xm.slice(i).t())*uk, 2);
    }
  }
  if(nl == 3){
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
      temp = norm((xm.slice(i).t())*uk, 2);
      d = d + 0.5*(tanh(temp)/temp);
    }
  }
  
  d = d/n;
  
  return Rcpp::wrap(d);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP computeb(SEXP varuk, SEXP varxm, SEXP varnl) {
  
  cube xm = as<cube>(varxm);
  int n = xm.n_slices;
  
  int nl = as<int>(varnl);
  
  vec uk = as<vec>(varuk);
  
  double b = 0;
  
  // Find X^t u and add its squared norm to the sum
  if(nl == 1){
    for (int i = 0; i < n; i++)
    {
      b = b + 2*pow(norm((xm.slice(i).t())*uk, 2), 2);
    }
  }
  if(nl == 2){
    for (int i = 0; i < n; i++)
    {
      b = b + 0.75*norm((xm.slice(i).t())*uk, 2);
    }
  }
  if(nl == 3){
    double temp = 0;
    for (int i = 0; i < n; i++)
    {
      temp = norm((xm.slice(i).t())*uk, 2);
      b = b + 0.25*((1 - pow(tanh(temp), 2))*temp - tanh(temp))/pow(temp, 3);
    }
  }
  
  b = b/n;
  
  return Rcpp::wrap(b);
}
