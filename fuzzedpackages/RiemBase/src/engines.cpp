#include "RcppArmadillo.h"
#include "riemfactory.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//////////////////////////////////////////////////////////////
// 1. pdist : pairwise distance
//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat engine_pdist(arma::cube data, std::string name){
  // XPtr<distPtr> xpfun = SetDistPtr(name);
  // distPtr fun = *xpfun;
  
  const int N = data.n_slices;
  arma::mat output(N,N,fill::zeros);
  arma::mat x,y;
  double distval;
  for (int i=0;i<(N-1);i++){
    x = data.slice(i);
    for (int j=(i+1);j<N;j++){
      y = data.slice(j);
      distval = riemfunc_dist(x,y,name);
      output(i,j) = distval;
      output(j,i) = distval;
    }
  }
  return(output);
}

// [[Rcpp::export]]
arma::mat engine_pdist_openmp(arma::cube data, std::string name, int nCores){
  // XPtr<distPtr> xpfun = SetDistPtr(name);
  // distPtr fun = *xpfun;
  
  const int N = data.n_slices;
  arma::mat output(N,N,fill::zeros);
  arma::mat x,y;
  double distval;
  
#ifdef _OPENMP
  #pragma omp parallel for num_threads(nCores) collapse(2) shared(output) private(x,y,distval)
  for (int i=0;i<(N-1);i++){
    for (int j=0;j<N;j++){
      if (i<j){
        x = data.slice(i);
        y = data.slice(j);
        distval = riemfunc_dist(x,y,name);
        output(i,j) = distval;
        output(j,i) = distval; 
      }
    }
  }
#else
  for (int i=0;i<(N-1);i++){
    for (int j=0;j<N;j++){
      if (i<j){
        x = data.slice(i);
        y = data.slice(j);
        distval = riemfunc_dist(x,y,name);
        output(i,j) = distval;
        output(j,i) = distval; 
      }
    }
  }
#endif
  return(output);
}

//////////////////////////////////////////////////////////////
// 2. pdist2 : pairwise data between two sets of data
//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat engine_pdist2(arma::cube data1, arma::cube data2, std::string name){
  // XPtr<distPtr> xpfun = SetDistPtr(name);
  // distPtr fun = *xpfun;
  
  const int M = data1.n_slices;
  const int N = data2.n_slices;
  arma::mat output(M,N,fill::zeros);
  arma::mat x,y;
  double distval;
  
  for (int i=0;i<M;i++){
    x = data1.slice(i);
    for (int j=0;j<N;j++){
      y = data2.slice(j);
      if (arma::norm(x-y,"fro")>1e-16){
        distval = riemfunc_dist(x,y,name);
        output(i,j) = distval;  
      }
    }
  }
  return(output);
}

// [[Rcpp::export]]
arma::mat engine_pdist2_openmp(arma::cube data1, arma::cube data2, std::string name, int nCores){
  // XPtr<distPtr> xpfun = SetDistPtr(name);
  // distPtr fun = *xpfun;
  
  const int M = data1.n_slices;
  const int N = data2.n_slices;
  arma::mat output(M,N,fill::zeros);
  arma::mat x,y;
  double distval;
  
#ifdef _OPENMP
  #pragma omp parallel for num_threads(nCores) collapse(2) shared(output) private(x,y,distval)
  for (int i=0;i<M;i++){
    for (int j=0;j<N;j++){
      x = data1.slice(i);
      y = data2.slice(j);
      if (arma::norm(x-y,"fro")>1e-16){
        distval = riemfunc_dist(x,y,name);
        output(i,j) = distval;  
      }
    }
  }
#else
  for (int i=0;i<M;i++){
    for (int j=0;j<N;j++){
      x = data1.slice(i);
      y = data2.slice(j);
      if (arma::norm(x-y,"fro")>1e-16){
        distval = riemfunc_dist(x,y,name);
        output(i,j) = distval;  
      }
    }
  }
#endif
  return(output);
}


//////////////////////////////////////////////////////////////
// 3. median : geometric median using Weiszfeld algorithm
//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List engine_median(arma::cube data, std::string name, int maxiter, double eps, arma::mat init){
  // get parameters
  const int N = data.n_slices;
  int iter = 0;
  
  // initialize
  arma::mat mold = init;
  arma::mat mnew;   mnew.copy_size(mold);  mnew.fill(0); 
  arma::mat dtmp;   dtmp.copy_size(mold);  dtmp.fill(0); // on TpM
  arma::cube tvecs; tvecs.copy_size(data); tvecs.fill(0);
  arma::vec normvec(N,fill::zeros);
  
  // let's iterate !
  arma::uvec nonsingular;
  double increment = 10000.00;
  while (increment > eps){
    // 1. compute log-pulled vectors and norm
    for (int i=0;i<N;i++){
      tvecs.slice(i) = riemfunc_log(mold, data.slice(i), name);
      normvec(i) = riemfunc_norm(mold, tvecs.slice(i), name);
    }
    // 2. find the one with non-singular distance
    nonsingular = arma::find(normvec>1e-10);  
    int M = nonsingular.n_elem;
    if (M==0){
      break;
    }
    // 3. update numerator
    dtmp.reset();
    dtmp.copy_size(mold); dtmp.fill(0);
    for (int j=0;j<M;j++){
      dtmp += tvecs.slice(nonsingular(j))/normvec(nonsingular(j));
    }
    // 4. update denominator
    double denom = 0.0;
    for (int j=0;j<M;j++){
      denom += 1/normvec(nonsingular(j));
    }
    dtmp /= denom;
    // 5. update using exponential map
    mnew = riemfunc_exp(mold, dtmp, 1.0, name);
    // 6. iterate condition
    increment = riemfunc_dist(mold, mnew, name);
    iter += 1;
    // 8. update
    mold = mnew;
    if (iter >= maxiter){
      break;
    }
  }
  
  arma::mat moutput = riemfunc_nearest(mold, name);
  return(Rcpp::List::create(Rcpp::Named("x")=moutput,
                            Rcpp::Named("iteration")=iter));
}


// [[Rcpp::export]]
Rcpp::List engine_median_openmp(arma::cube data, std::string name, int maxiter, double eps, int nCores, arma::mat init){
  // get parameters
  const int N = data.n_slices;
  int iter = 0;
  
  // initialize
  arma::mat mold = init;
  arma::mat mnew;   mnew.copy_size(mold);  mnew.fill(0); 
  arma::mat dtmp;   dtmp.copy_size(mold);  dtmp.fill(0); // on TpM
  arma::cube tvecs; tvecs.copy_size(data); tvecs.fill(0);
  arma::vec normvec(N,fill::zeros);
  
  // let's iterate !
  arma::uvec nonsingular;
  double increment = 10000.00;
  while (increment > eps){
    // 1. compute log-pulled vectors and norm
#ifdef _OPENMP
    #pragma omp parallel for num_threads(nCores) shared(tvecs, normvec, name, mold, data)
    for (int i=0;i<N;i++){
      tvecs.slice(i) = riemfunc_log(mold, data.slice(i), name);
      normvec(i) = riemfunc_norm(mold, tvecs.slice(i), name);
    }
#else
    for (int i=0;i<N;i++){
      tvecs.slice(i) = riemfunc_log(mold, data.slice(i), name);
      normvec(i) = riemfunc_norm(mold, tvecs.slice(i), name);
    }
#endif
    // 2. find the one with non-singular distance
    nonsingular = arma::find(normvec>1e-10);  
    int M = nonsingular.n_elem;
    if (M==0){
      break;
    }
    // 3. update numerator
    dtmp.reset();
    dtmp.copy_size(mold); dtmp.fill(0);
    for (int j=0;j<M;j++){
      dtmp += tvecs.slice(nonsingular(j))/normvec(nonsingular(j));
    }
    // 4. update denominator
    double denom = 0.0;
    for (int j=0;j<M;j++){
      denom += 1/normvec(nonsingular(j));
    }
    dtmp /= denom;
    // 5. update using exponential map
    mnew = riemfunc_exp(mold, dtmp, 1.0, name);
    // 6. iterate condition
    increment = riemfunc_dist(mold, mnew, name);
    iter += 1;
    // 8. update
    mold = mnew;
    if (iter >= maxiter){
      break;
    }
  }
  
  arma::mat moutput = riemfunc_nearest(mold, name);
  return(Rcpp::List::create(Rcpp::Named("x")=moutput,
                            Rcpp::Named("iteration")=iter));
}



//////////////////////////////////////////////////////////////
// 4. mean : Karcher Mean
//////////////////////////////////////////////////////////////
double engine_mean_eval(arma::mat tgt, arma::cube data, std::string name){
  const int N = data.n_slices;
  
  double output = 0.0;
  double tmpout = 0.0;
  for (int i=0;i<N;i++){
    tmpout = riemfunc_dist(tgt, data.slice(i), name);
    output += tmpout*tmpout;
  }
  return(output);
}
double engine_mean_stepsize(arma::mat mold, arma::mat grad, arma::cube data, std::string name, double evalmold){
  double stepsize = 1.0;
  arma::mat initexp = riemfunc_exp(mold, grad, -1.0*stepsize, name);
  double evalmnew = engine_mean_eval(initexp, data, name);
  if (evalmnew < evalmold){
    return(stepsize);
  } else {
    int iter = 0;
    while (evalmnew >= evalmold){
      stepsize = 0.8*stepsize;
      initexp  = riemfunc_exp(mold, grad, -1.0*stepsize, name);
      evalmnew = engine_mean_eval(initexp, data, name);
      if (iter >= 10){
        break;
      }
      iter += 1;
    }
    return(stepsize);
  }
}

arma::mat engine_extrinsicmean(arma::cube data, std::string name){
  int dnrow = data.n_rows;
  int dncol = data.n_cols;
  int nslice = data.n_slices;
  
  arma::vec testdata = riemfunc_equiv(data.slice(0),dnrow,dncol,name);
  int L = testdata.n_elem; 
  
  arma::mat Xext(L,nslice,fill::zeros);
  for (int i=0;i<nslice;i++){
    Xext.col(i) = riemfunc_equiv(data.slice(i),dnrow,dncol,name);
  }
  
  arma::vec extmean = arma::mean(Xext, 1); // find mean for each row.
  arma::mat inveqmu = riemfunc_invequiv(extmean,dnrow,dncol,name); // projected mean
  return(inveqmu);
}

// [[Rcpp::export]]
Rcpp::List engine_mean(arma::cube data, std::string name, int maxiter, double eps){
  // get parameters
  int N = data.n_slices;
  int iter = 0;
  
  // initialize
  arma::mat mold = engine_extrinsicmean(data, name); // extrinsic mean as an initializer
  arma::mat mnew;   mnew.copy_size(mold);  mnew.fill(0); 
  arma::cube tvecs; tvecs.copy_size(data); tvecs.fill(0);
  arma::mat dtmp; dtmp.copy_size(mold);  dtmp.fill(0); // on TpM

  // let's iterate !
  double sqnorm = 10000.00;
  while (sqnorm > eps){
    // 1. compute log-pulled vectors
    for (int i=0;i<N;i++){
      tvecs.slice(i) = riemfunc_log(mold, data.slice(i), name);
    }
    // 2. compute updating scheme
    dtmp = arma::mean(tvecs,2);
    // 3. update using exponential map and compute
    mnew = riemfunc_exp(mold, dtmp, 1.0, name);
    // 4. iteration : update sqnorm
    sqnorm = riemfunc_dist(mold,mnew,name);

    // 5. iteration : iter
    iter += 1;
    // 6. update others
    mold = mnew;
    if (iter >= maxiter){
      break;
    }
  }
  
  return(Rcpp::List::create(Rcpp::Named("x")=mold,
                            Rcpp::Named("iteration")=iter));
}

// [[Rcpp::export]]
Rcpp::List engine_mean_openmp(arma::cube data, std::string name, int maxiter, double eps, int nCores){
  // get parameters
  int N = data.n_slices;
  int iter = 0;
  
  // initialize
  arma::mat mold = engine_extrinsicmean(data, name); // extrinsic mean as an initializer
  arma::mat mnew;   mnew.copy_size(mold);  mnew.fill(0); 
  arma::cube tvecs; tvecs.copy_size(data); tvecs.fill(0);
  
  arma::mat dtmp; dtmp.copy_size(mold);  dtmp.fill(0); // on TpM
  
  // let's iterate !
  double sqnorm = 10000.00;
  while (sqnorm > eps){
    // 1. compute log-pulled vectors
#ifdef _OPENMP
    #pragma omp parallel for num_threads(nCores) shared(tvecs, mold, data, name)
    for (int i=0;i<N;i++){
      tvecs.slice(i) = riemfunc_log(mold, data.slice(i), name);
    }
#else
    for (int i=0;i<N;i++){
      tvecs.slice(i) = riemfunc_log(mold, data.slice(i), name);
    }
#endif
    
    // 2. compute updating scheme
    dtmp = arma::mean(tvecs,2);
    // 3. update using exponential map and compute
    mnew = riemfunc_exp(mold, dtmp, 1.0, name);
    // 4. iteration : update sqnorm
    sqnorm = riemfunc_dist(mold,mnew,name);
    
    // 5. iteration : iter
    iter += 1;
    // 6. update others
    mold = mnew;
    if (iter >= maxiter){
      break;
    }
  }
  
  return(Rcpp::List::create(Rcpp::Named("x")=mold,
                            Rcpp::Named("iteration")=iter));
}



//////////////////////////////////////////////////////////////
// 5. curvedist : L2 distance between two curves
//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double engine_curvedist(arma::cube data1, arma::cube data2, arma::vec vect, std::string name, double p){
  // get parameters
  int nrow = data1.n_rows;
  int ncol = data2.n_cols;
  int N = vect.n_elem;
  
  // 5-1. compute pairwise distance for each slice
  arma::mat tmp1(nrow,ncol,fill::zeros);
  arma::mat tmp2(nrow,ncol,fill::zeros);
  arma::vec dvec(N,fill::zeros);
  if (name=="intrinsic"){
    for (int i=0;i<N;i++){
      tmp1 = data1.slice(i);
      tmp2 = data2.slice(i);
      dvec(i) = std::pow(riemfunc_extdist(tmp1, tmp2, name), p);
    }
  } else {
    for (int i=0;i<N;i++){
      tmp1 = data1.slice(i);
      tmp2 = data2.slice(i);
      dvec(i) = std::pow(riemfunc_dist(tmp1, tmp2, name), p);
    }
  }
  
  // 5-2. summing up using trapezoidal rule
  double dt = 0.0;
  double tmpval = 0.0;
  for (int i=0;i<(N-1);i++){
    dt      = vect(i+1)-vect(i);
    tmpval += dt*(dvec(i+1)+dvec(i))/2.0;
  }
  double outval = std::pow(tmpval, 1.0/p);
  return(outval);
}
