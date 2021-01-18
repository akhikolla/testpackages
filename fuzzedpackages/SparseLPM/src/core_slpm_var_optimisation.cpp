#include "core_slpm_var.h"

void slpm_var::Optimisation()
{
  if (verbose) Rcpp::Rcout << "\nOptimisation has started ..." << std::endl;
  arma::wall_clock timer;
  timer.tic();
  unsigned int iter = 0;
  elbo_values_store.zeros(n_iter_max+1);
  elbo_values_store.at(iter) = elbo_value;
  ++iter;
  bool stop_condition = false;
  while (!stop_condition)
  {
    // ///// UPDATE LAMBDAS
    for (unsigned int i=0; i<M; ++i) for (unsigned int j=0; j<N; ++j) UpdateLambda(i,j);

    // ///// UPDATE As
    for (unsigned int k=0; k<K; ++k) UpdateA(k);

    // ///// UPDATE Bs
    for (unsigned int k=0; k<K; ++k) UpdateB(k);

    // ///// UPDATE DELTAS
    for (unsigned int k=0; k<K; ++k) UpdateDelta();

    // ///// UPDATE VARIATIONAL GAUSSIANS FOR SENDERS
    for (unsigned int i=0; i<M; ++i) for (unsigned int k=0; k<K; ++k) UpdateAlphaBetaU(i,k);
    
    // ///// UPDATE VARIATIONAL GAUSSIANS FOR RECEIVERS
    for (unsigned int j=0; j<N; ++j) for (unsigned int k=0; k<K; ++k) UpdateAlphaBetaV(j,k);
    
    ///// STORE ELBO VALUE AND PRINT OUT CURRENT ITERATION
    elbo_values_store.at(iter) = elbo_value;
    if (verbose) Rcpp::Rcout << "Elapsed Time " << floor(10*timer.toc())/10 << "\tEnd of iteration " << iter << "\t\tCurrent ELBO  =  " << elbo_value << std::endl;
    if (iter >= n_iter_max) 
    {
      Rcpp::Rcout << "WARNING: " << n_iter_max << " iterations reached" << std::endl;
      stop_condition = true;
    }
    if (elbo_value <= elbo_values_store.at(iter-1) + tol) stop_condition = true;
    ++iter;
  }
  elbo_values_store.resize(iter);
  if (verbose) Rcpp::Rcout << "... optimisation has terminated after " << floor(10*timer.toc())/10 << " seconds\n" << std::endl;
  if (debug_mode) CheckValues();////    DEBUG ONLY
}

