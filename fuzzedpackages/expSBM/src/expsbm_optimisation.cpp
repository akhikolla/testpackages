#include "expsbm.h"

void expsbm::Optimisation()
{
  if (verbose) Rcpp::Rcout << "\nOptimisation has started ..." << std::endl;
  arma::wall_clock timer;
  timer.tic();
  unsigned int iter = 0;
  elbo_values_store.zeros(n_iter_max+1);
  elbo_values_store.at(iter) = elbo_value;
  ++iter;
  bool stop_condition = false;
  if (n_iter_max > 0) while (!stop_condition)
  {
    ///// UPDATE PARAMETERS
    for (unsigned int g=0; g<K; ++g) for (unsigned int h=0; h<K; ++h) UpdateMu(g,h);
    for (unsigned int g=0; g<K; ++g) for (unsigned int h=0; h<K; ++h) UpdateNu(g,h);
    
    UpdateLambda();
    for (unsigned int i=0; i<N; ++i) UpdateZ(i);
    
    ///// UPDATE SUMMARIES, STATISTICS, AND ELBO
    EvaluateStatistics();
    EvaluateELBO();
    
    ///// STORE ELBO VALUE AND PRINT OUT CURRENT ITERATION
    elbo_values_store.at(iter) = elbo_value;
    if (verbose) Rcpp::Rcout << "Elapsed Time " << floor(10*timer.toc())/10 << "\tEnd of iteration " << iter << "\t\tCurrent ELBO  =  " << elbo_value << std::endl;
    if (iter >= n_iter_max) 
    {
      Rcpp::Rcout << "WARNING: " << n_iter_max << " iterations reached" << std::endl;
      stop_condition = true;
    }
    double relative_increment = (elbo_value - elbo_values_store.at(iter-1)) / elbo_value;
    if (relative_increment < 0) relative_increment = -relative_increment;
    if (relative_increment <= tol) stop_condition = true;
    ++iter;
  }
  elbo_values_store.resize(iter);
  if (verbose) Rcpp::Rcout << "... optimisation has terminated after " << floor(10*timer.toc())/10 << " seconds\n" << std::endl;
}


