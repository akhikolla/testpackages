#include "dblpm.h"

void dblpm::UpdateX (unsigned int i, double variance)
{
  unsigned int t, d;
  arma::vec::iterator j_itr;
  
  arma::rowvec x_prop(D);
  for (d=0; d<D; ++d) x_prop.at(d) = R::rnorm(x.at(i,d),sqrt(variance));
  
  double hr = 0;
  double dist, dist_prop, eta, eta_prop, prob, prob_prop;
  double delta_temp;
  for (j_itr=j_active.begin(); j_itr<j_active.end(); ++j_itr)
  {
    t = j_first_activity.at(*j_itr);
    for (t=j_first_activity.at(*j_itr); t<j_last_activity.at(*j_itr)+1; ++t)
    {
      if (i_activity_table.at(i,t) > 0)
      {
        if (t == j_first_activity.at(*j_itr)) delta_temp = delta;
        else delta_temp = y.at(i,*j_itr,t-1);
        dist = 0;
        dist_prop = 0;
        for (d=0; d<D; ++d) dist += std::pow(x.at(i,d)-w.at(*j_itr,d,t),2);
        for (d=0; d<D; ++d) dist_prop += std::pow(x_prop.at(i,d)-w.at(*j_itr,d,t),2);
        dist = sqrt(dist);
        dist_prop = sqrt(dist_prop);
        eta = gamma.at(t)*delta_temp + beta.at(t)*(1-delta_temp) - dist;
        eta_prop = gamma.at(t)*delta_temp + beta.at(t)*(1-delta_temp) - dist_prop;
        prob = exp(eta) / (1+exp(eta));
        prob_prop = exp(eta_prop) / (1+exp(eta_prop));
        if (y.at(i,*j_itr,t) == 1) hr += log(prob_prop) - log(prob);
        else hr += log(1-prob_prop) - log(1-prob);
      }
    }
  }
  
  for (d=0; d<D; d++) hr += -0.5 * taux * (std::pow(x_prop.at(d),2) - std::pow(x.at(i,d),2));
  
  double u = R::runif(0,1);
  if (log(u) < hr) for (d=0; d<D; d++) x.at(i,d) = x_prop.at(d);
}


void dblpm::UpdateW (unsigned int t, unsigned int j, double variance)
{
  unsigned int d;
  arma::vec::iterator i_itr;
  
  arma::vec w_prop (D, arma::fill::zeros);
  for (d=0; d<D; ++d) w_prop.at(d) = R::rnorm(w.at(j,d,t),sqrt(variance));
  
  double hr = 0;
  double dist, dist_prop, eta, eta_prop, prob, prob_prop;
  double delta_temp;
  if (t >= j_first_activity.at(j) && t < j_last_activity.at(j)+1)
  {
    for (i_itr=i_activity_list.at(t).begin(); i_itr<i_activity_list.at(t).end(); ++i_itr)
    {
      if (t == j_first_activity.at(j)) delta_temp = delta;
      else delta_temp = y.at(*i_itr,j,t-1);
      dist = 0;
      dist_prop = 0;
      for (d=0; d<D; ++d) dist += std::pow(x.at(*i_itr,d)-w.at(j,d,t),2);
      for (d=0; d<D; ++d) dist_prop += std::pow(x.at(*i_itr,d)-w_prop.at(d),2);
      dist = sqrt(dist);
      dist_prop = sqrt(dist_prop);
      eta = gamma.at(t)*delta_temp + beta.at(t)*(1-delta_temp) - dist;
      eta_prop = gamma.at(t)*delta_temp + beta.at(t)*(1-delta_temp) - dist_prop;
      prob = exp(eta) / (1+exp(eta));
      prob_prop = exp(eta_prop) / (1+exp(eta_prop));
      if (y.at(*i_itr,j,t) == 1) hr += log(prob_prop) - log(prob);
      else hr += log(1-prob_prop) - log(1-prob);
    }
  }
  
  for (d=0; d<D; ++d)
  {
    if (t == 0) hr += -0.5 * tauw0 * (std::pow(w_prop.at(d),2)-std::pow(w.at(j,d,t),2)) - 0.5 * tauw * (std::pow(w.at(j,d,t+1)-w_prop.at(d),2)-std::pow(w.at(j,d,t+1)-w.at(j,d,t),2));
    if (t > 0) if (t < T-1) hr += - 0.5 * tauw * (std::pow(w.at(j,d,t+1)-w_prop.at(d),2) + std::pow(w_prop.at(d)-w.at(j,d,t-1),2)) + 0.5 * tauw * (std::pow(w.at(j,d,t+1)-w.at(j,d,t),2) + std::pow(w.at(j,d,t)-w.at(j,d,t-1),2));
    if (t == T-1) hr += - 0.5 * tauw * (std::pow(w_prop.at(d)-w.at(j,d,t-1),2)) + 0.5 * tauw * (std::pow(w.at(j,d,t)-w.at(j,d,t-1),2));
  }
  
  double u = R::runif(0,1);
  if (log(u) < hr)
  {
    if (t == 0) for (d=0; d<D; ++d) w0_ss += std::pow(w_prop.at(d),2) - std::pow(w.at(j,d,t),2);
    if (t > 0) for (d=0; d<D; ++d) w_innovation_ss += std::pow(w_prop.at(d)-w.at(j,d,t-1),2) - std::pow(w.at(j,d,t)-w.at(j,d,t-1),2);
    if (t < T-1) for (d=0; d<D; ++d) w_innovation_ss += std::pow(w.at(j,d,t+1)-w_prop.at(d),2) - std::pow(w.at(j,d,t+1)-w.at(j,d,t),2);
    for (d=0; d<D; ++d) w.at(j,d,t) = w_prop.at(d);
  }
}

void dblpm::UpdateGamma (unsigned int t, double variance)
{
  unsigned int d;
  arma::vec::iterator i_itr, j_itr;
  
  double gamma_prop = R::rnorm(gamma.at(t),sqrt(variance));
  
  double hr = 0;
  double dist, eta, eta_prop, prob, prob_prop;
  double delta_temp;
  for (j_itr=j_active.begin(); j_itr<j_active.end(); ++j_itr)
  {
    if (t >= j_first_activity.at(*j_itr) && t < j_last_activity.at(*j_itr)+1)
    {
      for (i_itr=i_activity_list.at(t).begin(); i_itr<i_activity_list.at(t).end(); ++i_itr)
      {
        if (t == j_first_activity.at(*j_itr)) delta_temp = delta;
        else delta_temp = y.at(*i_itr,*j_itr,t-1);
        dist = 0;
        for (d=0; d<D; ++d) dist += std::pow(x.at(*i_itr,d)-w.at(*j_itr,d,t),2);
        dist = sqrt(dist);
        eta = gamma.at(t)*delta_temp + beta.at(t)*(1-delta_temp) - dist;
        eta_prop = gamma_prop*delta_temp + beta.at(t)*(1-delta_temp) - dist;
        prob = exp(eta) / (1+exp(eta));
        prob_prop = exp(eta_prop) / (1+exp(eta_prop));
        if (y.at(*i_itr,*j_itr,t) == 1) hr += log(prob_prop) - log(prob);
        else hr += log(1-prob_prop) - log(1-prob);
      }
    }
  }
  
  if (t == 0) hr += -0.5 * taugamma0 * (std::pow(gamma_prop,2)-std::pow(gamma.at(t),2)) - 0.5 * taugamma * (std::pow(gamma.at(t+1)-gamma_prop,2)-std::pow(gamma.at(t+1)-gamma.at(t),2));
  if (t > 0) if (t < T-1) hr += - 0.5 * taugamma * (std::pow(gamma.at(t+1)-gamma_prop,2) + std::pow(gamma_prop-gamma.at(t-1),2)) + 0.5 * taugamma * (std::pow(gamma.at(t+1)-gamma.at(t),2) + std::pow(gamma.at(t)-gamma.at(t-1),2));
  if (t == T-1) hr += - 0.5 * taugamma * (std::pow(gamma_prop-gamma.at(t-1),2)) + 0.5 * taugamma * (std::pow(gamma.at(t)-gamma.at(t-1),2));
  
  double u = R::runif(0,1);
  if (log(u) < hr)
  {
    if (t > 0) gamma_innovation_ss += std::pow(gamma_prop-gamma.at(t-1),2) - std::pow(gamma.at(t)-gamma.at(t-1),2);
    if (t < T-1) gamma_innovation_ss += std::pow(gamma.at(t+1)-gamma_prop,2) - std::pow(gamma.at(t+1)-gamma.at(t),2);
    gamma.at(t) = gamma_prop;
  }
}

void dblpm::UpdateBeta (unsigned int t, double variance)
{
  unsigned int d;
  arma::vec::iterator i_itr, j_itr;
  
  double beta_prop = R::rnorm(beta.at(t),sqrt(variance));
  
  double hr = 0;
  double dist, eta, eta_prop, prob, prob_prop;
  double delta_temp;
  for (j_itr=j_active.begin(); j_itr<j_active.end(); ++j_itr)
  {
    if (t >= j_first_activity.at(*j_itr) && t < j_last_activity.at(*j_itr)+1)
    {
      for (i_itr=i_activity_list.at(t).begin(); i_itr<i_activity_list.at(t).end(); ++i_itr)
      {
        if (t == j_first_activity.at(*j_itr)) delta_temp = delta;
        else delta_temp = y.at(*i_itr,*j_itr,t-1);
        dist = 0;
        for (d=0; d<D; ++d) dist += std::pow(x.at(*i_itr,d)-w.at(*j_itr,d,t),2);
        dist = sqrt(dist);
        eta = gamma.at(t)*delta_temp + beta.at(t)*(1-delta_temp) - dist;
        eta_prop = gamma.at(t)*delta_temp + beta_prop*(1-delta_temp) - dist;
        prob = exp(eta) / (1+exp(eta));
        prob_prop = exp(eta_prop) / (1+exp(eta_prop));
        if (y.at(*i_itr,*j_itr,t) == 1) hr += log(prob_prop) - log(prob);
        else hr += log(1-prob_prop) - log(1-prob);
      }
    }
  }
  
  if (t == 0) hr += -0.5 * taubeta0 * (std::pow(beta_prop,2)-std::pow(beta.at(t),2)) - 0.5 * taubeta * (std::pow(beta.at(t+1)-beta_prop,2)-std::pow(beta.at(t+1)-beta.at(t),2));
  if (t > 0) if (t < T-1) hr += - 0.5 * taubeta * (std::pow(beta.at(t+1)-beta_prop,2) + std::pow(beta_prop-beta.at(t-1),2)) + 0.5 * taubeta * (std::pow(beta.at(t+1)-beta.at(t),2) + std::pow(beta.at(t)-beta.at(t-1),2));
  if (t == T-1) hr += - 0.5 * taubeta * (std::pow(beta_prop-beta.at(t-1),2)) + 0.5 * taubeta * (std::pow(beta.at(t)-beta.at(t-1),2));
  
  double u = R::runif(0,1);
  if (log(u) < hr)
  {
    if (t > 0) beta_innovation_ss += std::pow(beta_prop-beta.at(t-1),2) - std::pow(beta.at(t)-beta.at(t-1),2);
    if (t < T-1) beta_innovation_ss += std::pow(beta.at(t+1)-beta_prop,2) - std::pow(beta.at(t+1)-beta.at(t),2);
    beta.at(t) = beta_prop;
  }
  
}

void dblpm::UpdateTauw ()
{
  tauw = R::rgamma(aw+M*D*(T-1)/2, 1/(bw+w_innovation_ss/2));
}

void dblpm::UpdateTauw0 ()
{
  tauw0 = R::rgamma(aw+M*D/2, 1/(bw+w0_ss/2));
}

void dblpm::UpdateTaugamma ()
{
  taugamma = R::rgamma(agamma+(T-1)/2, 1/(bgamma+gamma_innovation_ss/2));
}

void dblpm::UpdateTaugamma0 ()
{
  taugamma0 = R::rgamma(agamma+0.5,1/(bgamma+0.5*std::pow(gamma.at(0),2)));
}

void dblpm::UpdateTaubeta ()
{
  taubeta = R::rgamma(abeta+(T-1)/2, 1/(bbeta+beta_innovation_ss/2));
}

void dblpm::UpdateTaubeta0 ()
{
  taubeta0 = R::rgamma(abeta+0.5,1/(bbeta+0.5*std::pow(beta.at(0),2)));
}
