#include "dblpm.h"

void dblpm::EvaluateSumOfSquares()
{
  if (debug) Rcpp::Rcout << "dblpm::EvaluateSumOfSquares has been called" << std::endl;
  unsigned int t,j,d;
  w0_ss = 0.0;
  w_innovation_ss = 0.0;
  gamma_innovation_ss = 0.0;
  beta_innovation_ss = 0.0;
  for (j=0; j<M; ++j) for (d=0; d<D; ++d) w0_ss += std::pow(w.at(j,d,0),2);
  for (t=1; t<T; ++t)
  {
    for (j=0; j<M; ++j) for (d=0; d<D; ++d) w_innovation_ss += std::pow(w.at(j,d,t)-w.at(j,d,t-1),2);
    gamma_innovation_ss += std::pow(gamma.at(t)-gamma.at(t-1),2);
    beta_innovation_ss += std::pow(beta.at(t)-beta.at(t-1),2);
  }
  if (debug) Rcpp::Rcout << "dblpm::EvaluateSumOfSquares has been terminated" << std::endl;
}

void dblpm::SetNoMissingData()
{
  if (debug) Rcpp::Rcout << "dblpm::SetNoMissingData has been called" << std::endl;
  unsigned int l, t, i, j;
  out_degrees.set_size(N,T);
  out_degrees.fill(0);
  out_tot_degrees.set_size(N);
  out_tot_degrees.fill(0);
  in_degrees.set_size(M,T);
  in_degrees.fill(0);
  in_tot_degrees.set_size(M);
  in_tot_degrees.fill(0);
  for (l=0; l<L; ++l)
  {
    t = edgelist.at(l,0);
    i = edgelist.at(l,1);
    j = edgelist.at(l,2);
    out_degrees.at(i,t) += 1;
    out_tot_degrees.at(i) += 1;
    in_degrees.at(j,t) += 1;
    in_tot_degrees.at(j) += 1;
  }
  j_activity_table.set_size(M,T);
  j_activity_table.fill(1);
  i_activity_table.set_size(N,T);
  i_activity_table.fill(1);
  j_first_activity.set_size(M);
  j_first_activity.fill(0);
  j_last_activity.set_size(M);
  j_last_activity.fill(T-1);

  i_activity_list.set_size(T);
  for (t=0; t<T; ++t) i_activity_list.at(t) = arma::linspace<arma::vec>(0, N-1, N);
  N_active = N;
  M_active = M;
  i_active = arma::linspace<arma::vec>(0, N-1, N);
  j_active = arma::linspace<arma::vec>(0, M-1, M);
  if (debug) Rcpp::Rcout << "dblpm::SetNoMissingData has terminated" << std::endl;
}

void dblpm::FillActivity()
{
  if (debug) Rcpp::Rcout << "dblpm::FillActivity has been called" << std::endl;
  unsigned int t, l, j, i, index;
  out_degrees.set_size(N,T);
  out_degrees.fill(0);
  out_tot_degrees.set_size(N);
  out_tot_degrees.fill(0);
  in_degrees.set_size(M,T);
  in_degrees.fill(0);
  in_tot_degrees.set_size(M);
  in_tot_degrees.fill(0);
  j_activity_table.set_size(M,T);
  j_activity_table.fill(0);
  i_activity_table.set_size(N,T);
  i_activity_table.fill(0);
  for (l=0; l<L; ++l)
  {
    t = edgelist.at(l,0);
    i = edgelist.at(l,1);
    j = edgelist.at(l,2);
    out_degrees.at(i,t) += 1;
    out_tot_degrees.at(i) += 1;
    in_degrees.at(j,t) += 1;
    in_tot_degrees.at(j) += 1;
    j_activity_table.at(j,t) = 1;
    i_activity_table.at(i,t) = 1;
  }
  j_first_activity.set_size(M);
  j_first_activity.fill(T);
  j_last_activity.set_size(M);
  j_last_activity.fill(0);
  for (j=0; j<M; ++j)
  {
    index = 0;
    while (index < T)
    {
      if (j_activity_table.at(j,index) > 0)
      {
        j_first_activity.at(j) = index;
        index = T;
      }
      ++index;
    }
    index = 0;
    while (index < T)
    {
      if (j_activity_table.at(j,T-1-index) > 0)
      {
        j_last_activity.at(j) = T-1-index;
        index = T;
      }
      ++index;
    }
  }

  i_activity_list.set_size(T);
  unsigned int n_active;
  for (t=0; t<T; ++t)
  {
    n_active = 0;
    for (i=0; i<N; ++i) n_active += i_activity_table.at(i,t);
    i_activity_list.at(t).set_size(n_active);
    i_activity_list.at(t).fill(0);
    index = 0;
    for (i=0; i<N; ++i) if (i_activity_table.at(i,t) > 0)
    {
      i_activity_list.at(t).at(index) = i;
      ++index;
    }
  }

  N_active = 0;
  for (i=0; i<N; ++i) if (out_tot_degrees.at(i) > 0) ++N_active;
  M_active = 0;
  for (j=0; j<M; ++j) if (in_tot_degrees.at(j) > 0) ++M_active;

  index = 0;
  i_active.set_size(N_active);
  for (i=0; i<N; ++i) if (out_tot_degrees.at(i) > 0)
  {
    i_active.at(index) = i;
    ++index;
  }
  index = 0;
  j_active.set_size(M_active);
  for (j=0; j<M; ++j) if (in_tot_degrees.at(j) > 0)
  {
    j_active.at(index) = j;
    ++index;
  }
  if (debug) Rcpp::Rcout << "dblpm::FillActivity has terminated" << std::endl;
}

void dblpm::FillY()
{
  if (debug) Rcpp::Rcout << "dblpm::FillY has been called" << std::endl;
  y.set_size(N,M,T);
  y.fill(0);
  unsigned int l;
  for (l=0; l<L; ++l)
  {
    y.at(edgelist.at(l,1), edgelist.at(l,2), edgelist.at(l,0)) ++;
  }
  if (debug) Rcpp::Rcout << "dblpm::FillY has terminated" << std::endl;
}

void dblpm::Likelihood()
{
  if (debug) Rcpp::Rcout << "dblpm::Likelihood has been called" << std::endl;
  double res = 0;
  unsigned int t, d;
  arma::vec::iterator i_itr, j_itr;
  double dist, eta, prob;
  double delta_temp;
  for (j_itr=j_active.begin(); j_itr<j_active.end(); ++j_itr)
  {
    for (t=j_first_activity.at(*j_itr); t<j_last_activity.at(*j_itr)+1; ++t)
    {
      for (i_itr=i_activity_list.at(t).begin(); i_itr<i_activity_list.at(t).end(); ++i_itr)
      {
        //                    std::cout << "Adding contribution given by \ti=" << *i_itr << "\tj=" << *j_itr << "\tt=" << t << std::endl;
        if (t == j_first_activity.at(*j_itr)) delta_temp = delta;
        else delta_temp = y.at(*i_itr,*j_itr,t-1);
        dist = 0;
        for (d=0; d<D; ++d) dist += std::pow(x.at(*i_itr,d)-w.at(*j_itr,d,t),2);
        dist = sqrt(dist);
        eta = gamma.at(t)*delta_temp + beta.at(t)*(1-delta_temp) - dist;
        prob = exp(eta) / (1+exp(eta));
        if (y.at(*i_itr,*j_itr,t) == 1) res += log(prob);
        else res += log(1-prob);
      }
    }
  }
  likelihood_value = res;
  if (debug) Rcpp::Rcout << "dblpm::Likelihood has terminated" << std::endl;
}

void dblpm::Posterior ()
{
  if (debug) Rcpp::Rcout << "dblpm::Posterior has been called" << std::endl;
  unsigned int t, i, j, d;
  Likelihood();
  double res = likelihood_value;

  for (d=0; d<D; ++d) for (i=0; i<N; ++i) res += R::dnorm(x.at(i,d),0,1/sqrt(taux),1);
  for (d=0; d<D; ++d) for (j=0; j<M; ++j) res += R::dnorm(w.at(j,d,0),0,1/sqrt(tauw0),1);
  for (t=1; t<T; ++t) for (d=0; d<D; ++d) for (j=0; j<M; ++j) res += R::dnorm(w.at(j,d,t)-w.at(j,d,t-1),0,1/sqrt(tauw),1);
  res += R::dnorm(gamma.at(0),0,1/sqrt(taugamma0),1);
  for (t=1; t<T; ++t) res += R::dnorm(gamma.at(t)-gamma.at(t-1),0,1/sqrt(taugamma),1);
  res += R::dnorm(beta.at(0),0,1/sqrt(taubeta0),1);
  for (t=1; t<T; ++t) res += R::dnorm(beta.at(t)-beta.at(t-1),0,1/sqrt(taubeta),1);

  res += R::dgamma(tauw,aw,1/bw,1) + R::dgamma(tauw0,aw,1/bw,1);
  res += R::dgamma(taugamma,agamma,1/bgamma,1) + R::dgamma(taugamma0,agamma,1/bgamma,1);
  res += R::dgamma(taubeta,abeta,1/bbeta,1) + R::dgamma(taubeta0,abeta,1/bbeta,1);

  posterior_value = res;
  if (debug) Rcpp::Rcout << "dblpm::Posterior has terminated" << std::endl;
}
