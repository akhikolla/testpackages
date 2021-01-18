/* 
 Copyright (C) 2014 Matteo Fasiolo  matteo.fasiolo@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

#include "mvnfast.h"

/*
* Simulate random variables from a multivariate student-t distribution
*
* See ?rmixt() for a description of the arguments and output.
*/

RcppExport SEXP rmixtCpp(SEXP n_,  
                         SEXP mu_,  
                         SEXP sigma_,
                         SEXP df_,
                         SEXP indV_,
                         SEXP ncores_,
                         SEXP isChol_, 
                         SEXP retInd_,
                         SEXP A_)  
{ 
  using namespace Rcpp;
  
  try{
    
    uint32_t n = as<uint32_t>(n_);
    arma::mat mu = as<arma::mat>(mu_);  
    List sigma = List(sigma_);
    double df = as<double>(df_);
    IntegerVector indV = as<IntegerVector>(indV_);
    unsigned int  ncores = as<unsigned int>(ncores_); 
    bool isChol = as<bool>(isChol_); 
    bool retInd = as<bool>(retInd_);
    NumericMatrix A = NumericMatrix(A_);
    
    unsigned int d = mu.n_cols; // Dimension of random vectors
    unsigned int m = mu.n_rows; // Number of mixture components
    
    if(df < 0.0) stop("df must be positive.");
    if(n < 1) stop("n should be a positive integer");
    if(ncores == 0) stop("ncores has to be positive");
    if(d != A.ncol()) stop("mu.n_cols != A.ncol()");
    if(n != A.nrow()) stop("n != A.nrow()");
    if(mu.n_rows != m) stop("mu.n_rows != m");
    if(indV.length()!= n) stop("indV.length()!= n");
    
    // Get list of Cholesky decompositions of covariance matrices
    arma::mat tmpMat;
    std::vector< arma::mat > cholDec;
    for(int ii = 0; ii < m; ii++)
    {
      tmpMat = as<arma::mat>(wrap(sigma[ii]));
      
      if(d != tmpMat.n_cols) stop("mu.n_cols != sigma[ii].n_cols");
      if(d != tmpMat.n_rows) stop("mu.n_cols != sigma[ii].n_rows");
      
      // Calculate cholesky dec unless sigma is already a cholesky dec.
      if( isChol == false ) {
        cholDec.push_back( trimatu(arma::chol(tmpMat)) );
      }
      else{
        cholDec.push_back( trimatu(tmpMat) );
      }
    }
    
    // The A matrix that will be filled with firstly with standard normal rvs,
    // and finally with multivariate normal rvs.
    // We A wrap into a arma::mat "tmp" without making a copy.
    arma::mat tmp( A.begin(), A.nrow(), A.ncol(), false );
    
    RNGScope scope; // Declare RNGScope after the output in order to avoid a known Rcpp bug.
    
    // What I do to seed the sitmo::prng_engine is tricky. I produce "ncores" uniform numbers between 1 and the largest uint32_t,
    // which are the seeds. I put the first one in "coreSeed". If there is no support for OpenMP only this seed
    // will be used, as the computations will be sequential. If there is support for OpenMP, "coreSeed" will
    // be over-written, so that each core will get its own seed.
    NumericVector seeds = runif(ncores, 1.0, std::numeric_limits<uint32_t>::max());
    
#ifdef _OPENMP
#pragma omp parallel num_threads(ncores) if(ncores > 1)
{
#endif
  
  double acc;
  uint32_t irow, icol, ii, index;
  arma::rowvec work(d);
  
  uint32_t coreSeed = static_cast<uint32_t>(seeds[0]);
  
  // (Optionally) over-writing the seed here
#ifdef _OPENMP
  coreSeed = static_cast<uint32_t>( seeds[omp_get_thread_num()] );
#endif
  
  // Create the parallel random engine, normal distributions and chi-square distribution. 
  sitmo::prng_engine engine( coreSeed );
  boost::normal_distribution<> normal(0.0, 1.0);
  boost::random::chi_squared_distribution<> chisq( df );
  
  // Filling "out" with standard normal rvs
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (irow = 0; irow < n; irow++) 
    for (icol = 0; icol < d; icol++) 
      A(irow, icol) = normal(engine);
  
  // Multiplying "out"" by cholesky decomposition of covariance and adding the
  // mean to obtain the desired multivariate student-t data rvs.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for(irow = 0; irow < n; irow++)
  {
    index = indV[irow];
    
    for(icol = 0; icol < d; icol++)
    {
      acc = 0.0;
      
      for(ii = 0; ii <= icol; ii++) acc += tmp.at(irow, ii) * cholDec[index].at(ii, icol); 
      
      work.at(icol) = acc; 
      
    }
    
    tmp(arma::span(irow), arma::span::all) = work / sqrt(chisq(engine)/df) + mu.row(index);       
  }
  
#ifdef _OPENMP
}
#endif

// (Optionally) Add mixture indexes as attributes
if( retInd ){ A.attr("index") = indV+1; }

return R_NilValue;

  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return wrap(NA_REAL);
}
