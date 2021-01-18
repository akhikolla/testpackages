#include "dblpm.h"

void dblpm::MCMC(bool verbose)
{
  if (debug) Rcpp::Rcout << "dblpm::MCMC has been called" << std::endl;
  if (verbose) Rcpp::Rcout << "\nGibbs sampling has started ..." << std::endl;
  unsigned int i, j, t, d, index_inner, index_outer;
  arma::wall_clock timer;
  timer.tic();
  index_outer = 0;
  index_inner = 0;
  while (index_inner < niter)
  {
    
    for (i=0; i<N; ++i) for (d=0; d<D; ++d) UpdateX(i,x_var);
    for (t=0; t<T; ++t) for (j=0; j<M; ++j) for (d=0; d<D; ++d) UpdateW(t,j,w_var);
    for (t=0; t<T; ++t) UpdateGamma(t,gamma_var);
    for (t=0; t<T; ++t) UpdateBeta(t,beta_var);
    UpdateTauw();
    UpdateTauw0();
    UpdateTaugamma();
    UpdateTaugamma0();
    UpdateTaubeta();
    UpdateTaubeta0();
    
    if (index_outer > burnin) if (index_outer % thin == 0)
    {
      x_store.at(index_inner) = x;
      w_store.at(index_inner) = w;
      gamma_store.row(index_inner) = gamma.t();
      beta_store.row(index_inner) = beta.t();
      tauw_store.at(index_inner) = tauw;
      tauw0_store.at(index_inner) = tauw0;
      taugamma_store.at(index_inner) = taugamma;
      taugamma0_store.at(index_inner) = taugamma0;
      taubeta_store.at(index_inner) = taubeta;
      taubeta0_store.at(index_inner) = taubeta0;
      posterior_store.at(index_inner) = posterior_value;
      ++index_inner;
    }
    
    if (verbose) if (index_outer % 100 == 0) Rcpp::Rcout << "Elapsed Time " << floor(10*timer.toc())/10 << "\t\tEnd of iteration " << index_outer << " out of " << total_niter << std::endl;
    ++index_outer;
  }
  if (verbose) Rcpp::Rcout << "... Gibbs sampling has terminated after " << floor(10*timer.toc())/10 << " seconds\n" << std::endl;
  if (debug) Rcpp::Rcout << "dblpm::MCMC has terminated" << std::endl;
}

