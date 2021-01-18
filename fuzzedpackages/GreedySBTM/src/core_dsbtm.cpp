#include "core_dsbtm.h"

dsbtm::dsbtm()
{
  adj.load("data/network.bin");
  z.load("data/allocations.csv",arma::csv_ascii);
  max_n_iter = 100;
  verbose = false;
  T = z.n_rows;
  N = z.n_cols;
  Kup = z.max() + 1;
  eta0 = 0.5;
  zeta0 = 0.5;
  ap = 0.5;
  aq = 0.5;
  bp = 0.5;
  bq = 0.5;
  delta = 0.5;
  UpdateAllValues();
}

dsbtm::dsbtm(arma::cube adj_, arma::mat z_, unsigned int max_n_iter_, bool verbose_)
{
  adj = adj_;
  z = z_;
  max_n_iter = max_n_iter_;
  verbose = verbose_;
  T = z.n_rows;
  N = z.n_cols;
  Kup = z.max() + 1;
  eta0 = 0.5;
  zeta0 = 0.5;
  ap = 0.5;
  aq = 0.5;
  bp = 0.5;
  bq = 0.5;
  delta = 0.5;
  UpdateAllValues();
}

void dsbtm::Print()
{
  std::ostringstream strs;
  strs << "\nclass dsbtm\n";
  strs << "\nT\t=\t" << T;
  strs << "\nN\t=\t" << N;
  strs << "\n\nAdjacency cube:\n";
  adj.print(strs);
  strs << "\n\nAdjacency cube indicating observed edges:\n";
  active.print(strs);
  strs << "\n\nAllocations:\n";
  z.print(strs);
  strs << "\n\nGroup counts:\n";
  counts.print(strs);
  strs << "\n\nGroup total counts:\n";
  total_counts.t().print(strs);
  strs << "\n\nNon empty groups:\n";
  non_empty_groups.t().print(strs);
  strs << "\n\nEdges between groups counts (eta):\n";
  eta.print(strs);
  strs << "\n\nNon-edges between groups counts (zeta):\n";
  zeta.print(strs);
  strs << "\n\nFailed to create counts (u00):\n";
  u00.print(strs);
  strs << "\n\nCreate counts (u01):\n";
  u01.print(strs);
  strs << "\n\nDestroy counts (u10):\n";
  u10.print(strs);
  strs << "\n\nFailed to destroy counts (u11):\n";
  u11.print(strs);
  strs << "\n\nLog-prior value\t=\t" << prior_value << "\n";
  strs << "\n\nLog-likelihood value\t=\t" << likelihood_value << "\n";
  strs << "\n\nLog-posterior value\t=\t" << posterior_value << "\n";
  Rcpp::Rcout << strs.str() << std::endl << std::endl << std::endl;
}

void dsbtm::Summary()
{
  std::ostringstream strs;
  strs << "\nclass dsbtm\n";
  strs << "\n\nLog-prior value\t=\t" << prior_value << "\n";
  strs << "\n\nLog-likelihood value\t=\t" << likelihood_value << "\n";
  strs << "\n\nLog-posterior value\t=\t" << posterior_value << "\n";
  Rcpp::Rcout << strs.str() << std::endl << std::endl << std::endl;
}

void dsbtm::EvaluateActive()
{
  unsigned int t, i, j;
  active.zeros(N,N,T);
  for (t=0; t<T; ++t) for (i=0; i<N-1; ++i) for (j=i+1; j<N; ++j) if (z.at(t,i) > 0) if (z.at(t,j) > 0)
  {
    active.at(i,j,t) = 1;
    active.at(j,i,t) = 1;
  }
}

void dsbtm::EvaluateCountsN()
{
  unsigned int t, i;
  counts.zeros(T,Kup);
  total_counts.zeros(Kup);
  for (t=0; t<T; ++t) for (i=0; i<N; ++i) 
  {
    counts.at(t,z.at(t,i)) ++;
    total_counts.at(z.at(t,i)) ++;
  }
}

void dsbtm::EvaluateNonEmptyGroups()
{
  unsigned int g, index;
  K = 0;
  for (g=0; g<Kup; ++g) if (total_counts.at(g) > 0) K++;
  K_no_0 = K;
  if (total_counts.at(0) > 0) K_no_0--;
  non_empty_groups.zeros(K);
  non_empty_groups_no_0.zeros(K_no_0);
  index = 0;
  for (g=0; g<Kup; ++g) if (total_counts.at(g) > 0)
  {
    non_empty_groups.at(index) = g;
    ++index;
  }
  index = 0;
  for (g=1; g<Kup; ++g) if (total_counts.at(g) > 0)
  {
    non_empty_groups_no_0.at(index) = g;
    ++index;
  }
}

void dsbtm::EvaluateCountsPi()
{
  unsigned int t, i;
  transition_counts.zeros(Kup,Kup);
  for (i=0; i<N; ++i) for (t=1; t<T; ++t) transition_counts.at(z.at(t-1,i),z.at(t,i))++;
  transition_total_counts = sum(transition_counts,1);
}

void dsbtm::EvaluateCountsSBM()
{
  unsigned int t, i, j, g, h;
  eta.zeros(Kup,Kup);
  zeta.zeros(Kup,Kup);
  for (i=0; i<N-1; ++i) for (j=i+1; j<N; ++j) if (active.at(i,j,0) > 0)
  {
    g = z.at(0,i);
    h = z.at(0,j);
    eta.at(g,h) += adj.at(i,j,0);
    if (g != h) eta.at(h,g) += adj.at(i,j,0);
    zeta.at(g,h) += 1-adj.at(i,j,0);
    if (g != h) zeta.at(h,g) += 1-adj.at(i,j,0);
  }
  for (t=1; t<T; ++t) for (i=0; i<N-1; ++i) for (j=i+1; j<N; ++j) if (active.at(i,j,t-1) == 0) if (active.at(i,j,t) > 0)
  {
    g = z.at(t,i);
    h = z.at(t,j);
    eta.at(g,h) += adj.at(i,j,t);
    if (g != h) eta.at(h,g) += adj.at(i,j,t);
    zeta.at(g,h) += 1-adj.at(i,j,t);
    if (g != h) zeta.at(h,g) += 1-adj.at(i,j,t);
  }
}

void dsbtm::EvaluateCountsSBTM()
{
  unsigned int t, i, j, g, h;
  u00.zeros(Kup,Kup);
  u01.zeros(Kup,Kup);
  u10.zeros(Kup,Kup);
  u11.zeros(Kup,Kup);
  for (t=1; t<T; ++t) for (i=0; i<N-1; ++i) for (j=i+1; j<N; ++j) if (active.at(i,j,t-1) > 0) if (active.at(i,j,t) > 0)
  {
    g = z.at(t,i);
    h = z.at(t,j);
    if (adj.at(i,j,t-1) == 0) if (adj.at(i,j,t) == 0) 
    {
      u00.at(g,h) ++;
      if (g != h) u00.at(h,g) ++;
    }
    if (adj.at(i,j,t-1) == 0) if (adj.at(i,j,t) == 1) 
    {
      u01.at(g,h) ++;
      if (g != h) u01.at(h,g) ++;
    }
    if (adj.at(i,j,t-1) == 1) if (adj.at(i,j,t) == 0) 
    {
      u10.at(g,h) ++;
      if (g != h) u10.at(h,g) ++;
    }
    if (adj.at(i,j,t-1) == 1) if (adj.at(i,j,t) == 1) 
    {
      u11.at(g,h) ++;
      if (g != h) u11.at(h,g) ++;
    }
  }
}

void dsbtm::UpdateAllValues()
{
  EvaluateActive();
  EvaluateCountsN();
  EvaluateNonEmptyGroups();
  EvaluateCountsPi();
  EvaluateCountsSBM();
  EvaluateCountsSBTM();
  EvaluatePrior();
  EvaluateLikelihood();
  EvaluatePosterior();
}

void dsbtm::EvaluatePrior()
{
  arma::vec::iterator itr1, itr2;
  prior_value = 0;
  for (itr1=non_empty_groups.begin(); itr1<non_empty_groups.end(); ++itr1)
  {
    prior_value += counts.at(0,*itr1) * log(total_counts.at(*itr1)/(N*T));
    prior_value += lgamma(delta*K) - K*lgamma(delta) - lgamma(delta*K+transition_total_counts.at(*itr1));
    for (itr2=non_empty_groups.begin(); itr2<non_empty_groups.end(); ++itr2) prior_value += lgamma(delta+transition_counts.at(*itr1,*itr2));
  }
}

void dsbtm::EvaluateLikelihood()
{
  unsigned int g, h;
  likelihood_value = 0;
  for (g=1; g<Kup; ++g) for (h=g; h<Kup; ++h)
  {
    likelihood_value += lgamma(eta0+zeta0) - lgamma(eta0) - lgamma(zeta0);
    likelihood_value += lgamma(ap+bp) - lgamma(ap) - lgamma(bp);
    likelihood_value += lgamma(aq+bq) - lgamma(aq) - lgamma(bq);
    likelihood_value += lgamma(eta0 + eta.at(g,h)) + lgamma(zeta0 + zeta.at(g,h)) - lgamma(eta0 + zeta0 + eta.at(g,h) + zeta.at(g,h));
    likelihood_value += lgamma(ap + u01.at(g,h)) + lgamma(bp + u00.at(g,h)) - lgamma(ap + bp + u01.at(g,h) + u00.at(g,h));
    likelihood_value += lgamma(aq + u10.at(g,h)) + lgamma(bq + u11.at(g,h)) - lgamma(aq + bq + u10.at(g,h) + u11.at(g,h));
  }
}

void dsbtm::EvaluatePosterior()
{
  posterior_value = prior_value + likelihood_value;
}

void dsbtm::SetUpNodeInfoForUpdate(unsigned int t, unsigned int i)
{
  unsigned int j;
  n_edges_to_group.zeros(Kup);
  n_non_edges_to_group.zeros(Kup);
  if (t == 0) for (j=0; j<N; ++j) if (i != j) if (active.at(i,j,t) > 0) 
  {
    if (adj.at(i,j,t) > 0) n_edges_to_group.at(z.at(t,j))++;
    else n_non_edges_to_group.at(z.at(t,j))++;
  }
  if (t > 0) for (j=0; j<N; ++j) if (i != j)  if (active.at(i,j,t-1) == 0) if (active.at(i,j,t) > 0) 
  {
    if (adj.at(i,j,t) > 0) n_edges_to_group.at(z.at(t,j))++;
    else n_non_edges_to_group.at(z.at(t,j))++;
  }
  w00.zeros(Kup);
  w01.zeros(Kup);
  w10.zeros(Kup);
  w11.zeros(Kup);
  if (t > 0) for (j=0; j<N; ++j) if (i != j) if (active.at(i,j,t-1) > 0) if (active.at(i,j,t) > 0)
  {
    if (adj.at(i,j,t-1) == 0) if (adj.at(i,j,t) == 0) w00.at(z.at(t,j)) += 1;
    if (adj.at(i,j,t-1) == 0) if (adj.at(i,j,t) == 1) w01.at(z.at(t,j)) += 1;
    if (adj.at(i,j,t-1) == 1) if (adj.at(i,j,t) == 0) w10.at(z.at(t,j)) += 1;
    if (adj.at(i,j,t-1) == 1) if (adj.at(i,j,t) == 1) w11.at(z.at(t,j)) += 1;
  }
  prior_value_deltas.zeros(Kup);
  likelihood_value_deltas.zeros(Kup);
}

void dsbtm::EvaluatePriorDelta(unsigned int t, unsigned int i, unsigned int h)
{
  arma::vec::iterator itr;
  unsigned int g;
  g = z.at(t,i);
  double res = 0;
  if (g != h) if (g > 0) if (h > 0)
  {
    //  delta for the initial allocations
    if (t == 0)
    {
      if (total_counts.at(g) > 1) res += (counts.at(0,g)-1)*log(total_counts.at(g)-1) - counts.at(0,g)*log(total_counts.at(g));
      if (total_counts.at(h) > 0) res += (counts.at(0,h)+1)*log(total_counts.at(h)+1) - counts.at(0,h)*log(total_counts.at(h));
    } else {
      if (total_counts.at(g) > 1) res += counts.at(0,g)*log(total_counts.at(g)-1) - counts.at(0,g)*log(total_counts.at(g));
      if (total_counts.at(h) > 0) res += counts.at(0,h)*log(total_counts.at(h)+1) - counts.at(0,h)*log(total_counts.at(h));
    }
    
    //  delta for the part depending on transitions. Note that the complexity of this part is constant
    //  only if the number of group does not change, otherwise it is O(K)
    unsigned int K1;
    K1 = K;
    if (total_counts.at(g) == 1) K1--;
    if (total_counts.at(h) == 0) K1++;
    unsigned int previous = 0, next = 0;
    if (t > 0) previous = z.at(t-1,i);
    if (t < T-1) next = z.at(t+1,i);
    //    unsigned int case_number = 0;
    if (t > 0 && t < T-1)
    {
      if (g != previous && h != previous)
      {
        //            case_number = 3;
        res += lgamma(delta+transition_counts.at(previous,g)-1) - lgamma(delta+transition_counts.at(previous,g)) + lgamma(delta+transition_counts.at(previous,h)+1) - lgamma(delta+transition_counts.at(previous,h));
        res += lgamma(delta+transition_counts.at(g,next)-1) - lgamma(delta+transition_counts.at(g,next)) - lgamma(delta*K1+transition_total_counts.at(g)-1) + lgamma(delta*K+transition_total_counts.at(g));
        res += lgamma(delta+transition_counts.at(h,next)+1) - lgamma(delta+transition_counts.at(h,next)) - lgamma(delta*K1+transition_total_counts.at(h)+1) + lgamma(delta*K+transition_total_counts.at(h));
        if (K1 != K)
        {
          res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
          res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(previous)) + lgamma(delta*K+transition_total_counts(previous));
          for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != g) if (*itr != h) if (*itr != previous)
          {
            res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
          }
        }
      }//end of case 3
      else
      {
        if (g == previous)
        {
          if (h != next && g != next)
          {
            //                    case_number = 4;
            res += lgamma(delta+transition_counts.at(previous,g)-1) - lgamma(delta+transition_counts.at(previous,g)) +
              lgamma(delta+transition_counts.at(previous,next)-1) - lgamma(delta+transition_counts.at(previous,next)) +
              lgamma(delta+transition_counts.at(previous,h)+1) - lgamma(delta+transition_counts.at(previous,h)) -
              lgamma(delta*K1+transition_total_counts.at(previous)-1) + lgamma(delta*K+transition_total_counts.at(previous));
            res += lgamma(delta+transition_counts.at(h,next)+1) - lgamma(delta+transition_counts.at(h,next)) - lgamma(delta*K1+transition_total_counts.at(h)+1) + lgamma(delta*K+transition_total_counts.at(h));
            if (K1 != K)
            {
              res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
              for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != h) if (*itr != previous)
              {
                res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
              }
            }
          }//end of case 4
          else if (g == next)
          {
            //                    case_number = 6;
            res += lgamma(delta+transition_counts.at(previous,g)-2) - lgamma(delta+transition_counts.at(previous,g)) +
              lgamma(delta+transition_counts.at(previous,h)+1) - lgamma(delta+transition_counts.at(previous,h)) -
              lgamma(delta*K1+transition_total_counts.at(previous)-1) + lgamma(delta*K+transition_total_counts.at(previous));
            res += lgamma(delta+transition_counts.at(h,next)+1) - lgamma(delta+transition_counts.at(h,next)) -
              lgamma(delta*K1+transition_total_counts.at(h)+1) + lgamma(delta*K+transition_total_counts.at(h));
            if (K1 != K)
            {
              res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
              for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != h) if (*itr != previous)
              {
                res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
              }
            }
          }//end of case 6
          else if (h == next)
          {
            //                    case_number = 8;
            res += lgamma(delta+transition_counts.at(previous,previous)-1) - lgamma(delta+transition_counts.at(previous,previous)) -
              lgamma(delta*K1+transition_total_counts.at(previous)-1) + lgamma(delta*K+transition_total_counts.at(previous));
            res += lgamma(delta+transition_counts.at(next,next)+1) - lgamma(delta+transition_counts.at(next,next)) -
              lgamma(delta*K1+transition_total_counts.at(next)+1) + lgamma(delta*K+transition_total_counts.at(next));
            if (K1 != K)
            {
              res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
              for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != next) if (*itr != previous)
              {
                res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
              }
            }
          }//end of case 8
        }//end if g==previous
        else if (h == previous)
        {
          if (g != next && h != next)
          {
            //                    case_number = 5;
            res += lgamma(delta+transition_counts.at(previous,g)-1) - lgamma(delta+transition_counts.at(previous,g)) +
              lgamma(delta+transition_counts.at(previous,next)+1) - lgamma(delta+transition_counts.at(previous,next)) +
              lgamma(delta+transition_counts.at(previous,h)+1) - lgamma(delta+transition_counts.at(previous,h)) -
              lgamma(delta*K1+transition_total_counts.at(previous)+1) + lgamma(delta*K+transition_total_counts.at(previous));
            res += lgamma(delta+transition_counts.at(g,next)-1) - lgamma(delta+transition_counts.at(g,next)) -
              lgamma(delta*K1+transition_total_counts.at(g)-1) + lgamma(delta*K+transition_total_counts.at(g));
            if (K1 != K)
            {
              res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
              for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != g) if (*itr != previous)
              {
                res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
              }
            }
          }//end of case 5
          else if (h == next)
          {
            //                    case_number = 7;
            res += lgamma(delta+transition_counts.at(previous,g)-1) - lgamma(delta+transition_counts.at(previous,g)) +
              lgamma(delta+transition_counts.at(previous,h)+2) - lgamma(delta+transition_counts.at(previous,h)) -
              lgamma(delta*K1+transition_total_counts.at(previous)+1) + lgamma(delta*K+transition_total_counts.at(previous));
            res += lgamma(delta+transition_counts.at(g,next)-1) - lgamma(delta+transition_counts.at(g,next)) -
              lgamma(delta*K1+transition_total_counts.at(g)-1) + lgamma(delta*K+transition_total_counts.at(g));
            if (K1 != K)
            {
              res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
              for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != g) if (*itr != previous)
              {
                res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
              }
            }
          }//end of case 7
          else if (g == next)
          {
            //                    case_number = 9;
            res += lgamma(delta+transition_counts.at(previous,previous)+1) - lgamma(delta+transition_counts.at(previous,previous)) -
              lgamma(delta*K1+transition_total_counts.at(previous)+1) + lgamma(delta*K+transition_total_counts.at(previous));
            res += lgamma(delta+transition_counts.at(next,next)-1) - lgamma(delta+transition_counts.at(next,next)) -
              lgamma(delta*K1+transition_total_counts.at(next)-1) + lgamma(delta*K+transition_total_counts.at(next));
            if (K1 != K)
            {
              res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
              for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != next) if (*itr != previous)
              {
                res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
              }
            }
          }//end of case 9
        }//end of h==previous
      }//end of else for cases 4-6-8 and 5-7-9
    }//end of if (t > 0 && t < T-1)
    else
    {
      if (t == 0)
      {
        //            case_number = 1;
        res += lgamma(delta+transition_counts.at(g,next)-1) - lgamma(delta+transition_counts.at(g,next)) - lgamma(delta*K1+transition_total_counts.at(g)-1) + lgamma(delta*K+transition_total_counts.at(g));
        res += lgamma(delta+transition_counts.at(h,next)+1) - lgamma(delta+transition_counts.at(h,next)) - lgamma(delta*K1+transition_total_counts.at(h)+1) + lgamma(delta*K+transition_total_counts.at(h));
        if (K1 != K)
        {
          res += 2*( lgamma(delta*K1) - lgamma(delta*K) );
          for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != g) if (*itr != h)
          {
            res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
          }
        }
      }//end of case 1
      else if (t == T-1)
      {
        //            case_number = 2;
        res += lgamma(delta+transition_counts.at(previous,g)-1) - lgamma(delta+transition_counts.at(previous,g)) + lgamma(delta+transition_counts.at(previous,h)+1) - lgamma(delta+transition_counts.at(previous,h));
        if (K1 != K)
        {
          res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(previous)) + lgamma(delta*K+transition_total_counts(previous));
          for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr) if (*itr != previous)
          {
            res += lgamma(delta*K1) - lgamma(delta*K) - lgamma(delta*K1+transition_total_counts.at(*itr)) + lgamma(delta*K+transition_total_counts(*itr));
          }
        }
      }//end of case 2
    }//end of else for cases 1-2
  }//end of statement if (g != h) 
  prior_value_deltas.at(h) = res;
}

void dsbtm::EvaluateLikelihoodDelta(unsigned int t, unsigned int i, unsigned int h)
{
  arma::vec::iterator itr;
  unsigned int g;
  g = z.at(t,i);
  double res = 0;
  if (g != h) if (g > 0) if (h > 0)
  {
    for (itr=non_empty_groups_no_0.begin(); itr<non_empty_groups_no_0.end(); ++itr) if (*itr != g) if (*itr != h)
    {
      // SBM part
      res += BernoulliMarginalDelta(eta.at(g,*itr)+zeta.at(g,*itr)-n_edges_to_group.at(*itr)-n_non_edges_to_group.at(*itr), eta.at(g,*itr)-n_edges_to_group.at(*itr), eta.at(g,*itr)+zeta.at(g,*itr), eta.at(g,*itr), eta0, zeta0);
      res += BernoulliMarginalDelta(eta.at(h,*itr)+zeta.at(h,*itr)+n_edges_to_group.at(*itr)+n_non_edges_to_group.at(*itr), eta.at(h,*itr)+n_edges_to_group.at(*itr), eta.at(h,*itr)+zeta.at(h,*itr), eta.at(h,*itr), eta0, zeta0);
      // SBTM part
      res += BernoulliMarginalDelta(u00.at(g,*itr)-w00.at(*itr)+u01.at(g,*itr)-w01.at(*itr), u01.at(g,*itr)-w01.at(*itr), u00.at(g,*itr)+u01.at(g,*itr), u01.at(g,*itr), ap, bp);
      res += BernoulliMarginalDelta(u10.at(g,*itr)-w10.at(*itr)+u11.at(g,*itr)-w11.at(*itr), u10.at(g,*itr)-w10.at(*itr), u10.at(g,*itr)+u11.at(g,*itr), u10.at(g,*itr), aq, bq);
      res += BernoulliMarginalDelta(u00.at(h,*itr)+w00.at(*itr)+u01.at(h,*itr)+w01.at(*itr), u01.at(h,*itr)+w01.at(*itr), u00.at(h,*itr)+u01.at(h,*itr), u01.at(h,*itr), ap, bp);
      res += BernoulliMarginalDelta(u10.at(h,*itr)+w10.at(*itr)+u11.at(h,*itr)+w11.at(*itr), u10.at(h,*itr)+w10.at(*itr), u10.at(h,*itr)+u11.at(h,*itr), u10.at(h,*itr), aq, bq);
    }
    // SBM part
    res += BernoulliMarginalDelta(eta.at(g,g)+zeta.at(g,g)-n_edges_to_group.at(g)-n_non_edges_to_group.at(g), eta.at(g,g)-n_edges_to_group.at(g), eta.at(g,g)+zeta.at(g,g), eta.at(g,g), eta0, zeta0);
    res += BernoulliMarginalDelta(eta.at(g,h)+zeta.at(g,h)+n_edges_to_group.at(g)+n_non_edges_to_group.at(g)-n_edges_to_group.at(h)-n_non_edges_to_group.at(h), eta.at(g,h)-n_edges_to_group.at(h)+n_edges_to_group.at(g), eta.at(g,h)+zeta.at(g,h), eta.at(g,h), eta0, zeta0);
    res += BernoulliMarginalDelta(eta.at(h,h)+zeta.at(h,h)+n_edges_to_group.at(h)+n_non_edges_to_group.at(h), eta.at(h,h)+n_edges_to_group.at(h), eta.at(h,h)+zeta.at(h,h), eta.at(h,h), eta0, zeta0);
    // SBTM part
    res += BernoulliMarginalDelta(u00.at(g,g)-w00.at(g)+u01.at(g,g)-w01.at(g), u01.at(g,g)-w01.at(g), u00.at(g,g)+u01.at(g,g), u01.at(g,g), ap, bp);
    res += BernoulliMarginalDelta(u10.at(g,g)-w10.at(g)+u11.at(g,g)-w11.at(g), u10.at(g,g)-w10.at(g), u10.at(g,g)+u11.at(g,g), u10.at(g,g), aq, bq);
    res += BernoulliMarginalDelta(u00.at(g,h)+w00.at(g)-w00.at(h)+u01.at(g,h)-w01.at(h)+w01.at(g), u01.at(g,h)-w01.at(h)+w01.at(g), u00.at(g,h)+u01.at(g,h), u01.at(g,h), ap, bp);
    res += BernoulliMarginalDelta(u10.at(g,h)+w10.at(g)-w10.at(h)+u11.at(g,h)-w11.at(h)+w11.at(g), u10.at(g,h)-w10.at(h)+w10.at(g), u10.at(g,h)+u11.at(g,h), u10.at(g,h), aq, bq);
    res += BernoulliMarginalDelta(u00.at(h,h)+w00.at(h)+u01.at(h,h)+w01.at(h), u01.at(h,h)+w01.at(h), u00.at(h,h)+u01.at(h,h), u01.at(h,h), ap, bp);
    res += BernoulliMarginalDelta(u10.at(h,h)+w10.at(h)+u11.at(h,h)+w11.at(h), u10.at(h,h)+w10.at(h), u10.at(h,h)+u11.at(h,h), u10.at(h,h), aq, bq);
  }
  // Make sure the vector likelihood_value_deltas is initialised first!!
  likelihood_value_deltas.at(h) = res;
}

void dsbtm::Move(unsigned int t, unsigned int i, unsigned int h)
{
  unsigned int g;
  arma::vec::iterator itr;
  prior_value += prior_value_deltas.at(h);
  likelihood_value += likelihood_value_deltas.at(h);
  EvaluatePosterior();  
  // SetUpNodeInfoForUpdate(t,i);// remember to call this line somewhere before doing the move, this line is repeated here in case you want to DEBUG
  g = z.at(t,i);
  if (g == h) throw std::runtime_error("Attempting to move a node to the same group it is currently allocated to.");
  
  // Updating edge counts for the SBM part
  for (itr=non_empty_groups_no_0.begin(); itr<non_empty_groups_no_0.end(); ++itr)  if (*itr != g) if (*itr != h)
  {
    eta.at(g,*itr) -= n_edges_to_group.at(*itr);
    eta.at(*itr,g) -= n_edges_to_group.at(*itr);
    zeta.at(g,*itr) -= n_non_edges_to_group.at(*itr);
    zeta.at(*itr,g) -= n_non_edges_to_group.at(*itr);
    eta.at(h,*itr) += n_edges_to_group.at(*itr);
    eta.at(*itr,h) += n_edges_to_group.at(*itr);
    zeta.at(h,*itr) += n_non_edges_to_group.at(*itr);
    zeta.at(*itr,h) += n_non_edges_to_group.at(*itr);
  }
  eta.at(g,g) += - n_edges_to_group.at(g);
  eta.at(g,h) += n_edges_to_group.at(g) - n_edges_to_group.at(h);
  eta.at(h,g) += n_edges_to_group.at(g) - n_edges_to_group.at(h);
  eta.at(h,h) += n_edges_to_group.at(h);
  zeta.at(g,g) += - n_non_edges_to_group.at(g);
  zeta.at(g,h) += n_non_edges_to_group.at(g) - n_non_edges_to_group.at(h);
  zeta.at(h,g) += n_non_edges_to_group.at(g) - n_non_edges_to_group.at(h);
  zeta.at(h,h) += n_non_edges_to_group.at(h);
  
  // Updating creation and destruction counts
  for (itr=non_empty_groups.begin(); itr<non_empty_groups.end(); ++itr)  if (*itr != g) if (*itr != h)
  {
    u00.at(g,*itr) -= w00.at(*itr);
    u01.at(g,*itr) -= w01.at(*itr);
    u10.at(g,*itr) -= w10.at(*itr);
    u11.at(g,*itr) -= w11.at(*itr);
    u00.at(*itr,g) -= w00.at(*itr);
    u01.at(*itr,g) -= w01.at(*itr);
    u10.at(*itr,g) -= w10.at(*itr);
    u11.at(*itr,g) -= w11.at(*itr);
    u00.at(h,*itr) += w00.at(*itr);
    u01.at(h,*itr) += w01.at(*itr);
    u10.at(h,*itr) += w10.at(*itr);
    u11.at(h,*itr) += w11.at(*itr);
    u00.at(*itr,h) += w00.at(*itr);
    u01.at(*itr,h) += w01.at(*itr);
    u10.at(*itr,h) += w10.at(*itr);
    u11.at(*itr,h) += w11.at(*itr);
  }
  u00.at(g,g) -= w00.at(g);
  u01.at(g,g) -= w01.at(g);
  u10.at(g,g) -= w10.at(g);
  u11.at(g,g) -= w11.at(g);
  u00.at(g,h) += w00.at(g) - w00.at(h);
  u01.at(g,h) += w01.at(g) - w01.at(h);
  u10.at(g,h) += w10.at(g) - w10.at(h);
  u11.at(g,h) += w11.at(g) - w11.at(h);
  u00.at(h,g) += w00.at(g) - w00.at(h);
  u01.at(h,g) += w01.at(g) - w01.at(h);
  u10.at(h,g) += w10.at(g) - w10.at(h);
  u11.at(h,g) += w11.at(g) - w11.at(h);
  u00.at(h,h) += w00.at(h);
  u01.at(h,h) += w01.at(h);
  u10.at(h,h) += w10.at(h);
  u11.at(h,h) += w11.at(h);
  
  // Moving the node to the new group
  z.at(t,i) = h;

  // Allocation counts and non-empty groups
  counts.at(t,g) -= 1;
  counts.at(t,h) += 1;
  total_counts.at(g) -= 1;
  total_counts.at(h) += 1;
  EvaluateNonEmptyGroups();
  unsigned int previous_allocation, next_allocation;
  if (t > 0)
  {
    previous_allocation = z.at(t-1,i);
    transition_counts(previous_allocation,g) -= 1;
    transition_counts(previous_allocation,h) += 1;
  }
  if (t < T-1) {
    next_allocation = z.at(t+1,i);
    transition_counts(g,next_allocation) -= 1;
    transition_counts(h,next_allocation) += 1;
    transition_total_counts.at(g) -= 1;
    transition_total_counts.at(h) += 1;
  }
}

void dsbtm::GreedyMove(unsigned int t, unsigned int i)
{
  SetUpNodeInfoForUpdate(t,i);// This will prepare and initialise all of the quantities needed for the update

  arma::vec::iterator itr;
  unsigned int g = z.at(t,i);
  unsigned int h;
  arma::vec deltas;
  
  // Evaluate the deltas for the moves to each of the non-empty groups
  deltas.zeros(Kup);
  for (itr=non_empty_groups_no_0.begin(); itr<non_empty_groups_no_0.end(); ++itr) if (*itr != g)
  {
    h = *itr;
    EvaluatePriorDelta(t,i,h);
    EvaluateLikelihoodDelta(t,i,h);
    deltas.at(h) += prior_value_deltas.at(h) + likelihood_value_deltas.at(h);
  }
  
  // The following makes sure that the delta is evaluated also for an empty group
  unsigned int first_empty_group = Kup;
  bool stop_condition = false;
  h = 1;
  while (!stop_condition)
  {
    if (counts.at(h) == 0) 
    {
      first_empty_group = h;
      stop_condition = true;
    }
    ++h;
    if (h >= Kup) stop_condition = true;
  }
  h = first_empty_group;
  if (h < Kup)
  {
    EvaluatePriorDelta(t,i,h);
    EvaluateLikelihoodDelta(t,i,h);
    deltas.at(h) += prior_value_deltas.at(h) + likelihood_value_deltas.at(h);
  }
  
  // Move the node to the group that gives the best ICL increase
  unsigned int h_best = g;
  double best_delta = 0;
  for (h=1; h<Kup; ++h) if (deltas.at(h) > best_delta)
  {
    h_best = h;
    best_delta = deltas.at(h);
  }
  
  // Perform the best move
  if (g != h_best) 
  {
    if (verbose) Rcpp::Rcout << "Moving node (" << t << "," << i << ") from group " << g << " to group " << h_best << std::endl;
    Move(t,i,h_best);
  }
  else if (verbose) Rcpp::Rcout << "Node (" << t << "," << i << ") not moved " << std::endl;
}

void dsbtm::GreedyOptimisation()
{
  if (verbose) Rcpp::Rcout << "\n\nGreedy optimisation of ICLex started\n" << std::endl;
  unsigned int iter, index, t, i;
  arma::wall_clock timer;
  bool stop_condition;
  //  Define the set of labels for the nodes to randomise the update order
  unsigned int pool_length = T*N;
  arma::vec pool = arma::linspace<arma::vec>(0, pool_length-1, pool_length);
  arma::mat index_legend(pool_length,2,arma::fill::zeros);
  index = 0;
  for (t=0; t<T; ++t) for (i=0; i<N; ++i)
  {
    index_legend(index,0) = t;
    index_legend(index,1) = i;
    ++index;
  }// each row of index_legend identifies an allocation that may be changed/updated
  greedy_icl_store.zeros(max_n_iter*T*N+1);
  greedy_icl_store.at(0) = posterior_value;
  stop_condition = false;
  iter = 0;
  timer.tic();
  while (!stop_condition && iter < max_n_iter)
  {
    pool = RandomShuffle(pool);
    for (index=0; index<T*N; ++index)
    {
      t = index_legend.at(pool.at(index),0);
      i = index_legend.at(pool.at(index),1);
      if (z.at(t,i) != 0) GreedyMove(t,i);
      greedy_icl_store.at(1+T*N*iter+index) = posterior_value;
    }// end of iteration
    if (verbose) Rcpp::Rcout << "\nIteration " << iter << " ended after " << floor(10*timer.toc())/10 << " seconds " << std::endl;
    index = 0;
    ++iter;
    if (iter > 0) if (posterior_value <= greedy_icl_store.at(1+T*N*(iter-1))) stop_condition = true;// stop if no allocations were changed in the last iteration
  }// end of optimisation
  greedy_icl_store = greedy_icl_store.subvec(0,T*N*iter+index);
  if (verbose) Rcpp::Rcout << "\nGreedy optimisation finished after " << floor(10*timer.toc())/10 << " seconds " << std::endl;
}

void dsbtm::MergeUpdates()
{
  if (verbose)
  {
    Rcpp::Rcout << "\nStarting hierarchical clustering greedy steps" << std::endl;
    Rcpp::Rcout << "Current non-empty groups are:" << std::endl;
    std::ostringstream strs;
    non_empty_groups.t().print(strs);
    Rcpp::Rcout << strs.str() << std::endl;
  }
  unsigned int t, i, index, g, h;
  arma::vec::iterator itr1, itr2;
  double current_posterior;
  arma::mat current_z;
  unsigned int pool_length;
  arma::vec pool;
  arma::mat index_legend;
  double stop_condition;
  stop_condition = false;
  while(!stop_condition && K_no_0 > 1)
  {
    pool_length = K_no_0*(K_no_0-1)/2;
    pool = arma::linspace<arma::vec>(0, pool_length-1, pool_length);
    index_legend.zeros(pool_length,2);
    index = 0;
    for (itr1=non_empty_groups_no_0.begin(); itr1<non_empty_groups_no_0.end(); ++itr1)
    {
      for (itr2=non_empty_groups_no_0.begin(); itr2<non_empty_groups_no_0.end(); ++itr2)
      {
        if (*itr1 < *itr2)
        {
          index_legend.at(index,0) = *itr1;
          index_legend.at(index,1) = *itr2;
          index++;
        }
      }
    }
    pool = RandomShuffle(pool);
    stop_condition = true;
    for (itr1=pool.begin(); itr1<pool.end(); ++itr1)
    {
      g = index_legend.at(*itr1,0);
      h = index_legend.at(*itr1,1);
      if (g > 0) if (h > g)
      {
        current_posterior = posterior_value;
        current_z = z;
        for (t=0; t<T; ++t) for (i=0; i<N; ++i) if (z.at(t,i) == g) z.at(t,i) = h;
        UpdateAllValues();
        if (verbose) Rcpp::Rcout << "Testing merging of groups " << g << " and " << h << ": delta is equal to " << posterior_value - current_posterior << std::endl;
        if (current_posterior > posterior_value)
        {
          z = current_z;
          UpdateAllValues();
        }
        else
        {
          stop_condition = false;// keep better z and start over again!
          break;
        }
      }// end of condition if (g > 0) if (h > g)
    }// end of for (itr1=pool.begin(); itr1<pool.end(); ++itr1)
  }// end of while(!stop_condition && K_no_0 > 1)
  if (verbose) Rcpp::Rcout << "\nHierarchical clustering procedure ended.\n" << std::endl;
}

void dsbtm::DebugCheckAllValues()
{
  arma::mat counts_check = counts;
  arma::mat transition_counts_check = transition_counts;
  arma::mat eta_check = eta;
  arma::mat zeta_check = zeta;
  arma::mat u00_check = u00;
  arma::mat u01_check = u01;
  arma::mat u10_check = u10;
  arma::mat u11_check = u11;
  double prior_value_check = prior_value;
  double likelihood_value_check = likelihood_value;
  UpdateAllValues();
  if (verbose) 
  {
    Rcpp::Rcout << "\nDEBUG\n" << std::endl;
    Rcpp::Rcout << "Error on counts\t=\t" << accu(abs(counts_check - counts)) << std::endl;
    Rcpp::Rcout << "Error on transition counts\t=\t" << accu(abs(transition_counts_check - transition_counts)) << std::endl;
    Rcpp::Rcout << "Error on eta\t=\t" << accu(abs(eta_check - eta)) << std::endl;
    Rcpp::Rcout << "Error on zeta\t=\t" << accu(abs(zeta_check - zeta)) << std::endl;
    Rcpp::Rcout << "Error on u00\t=\t" << accu(abs(u00_check - u00)) << std::endl;
    Rcpp::Rcout << "Error on u01\t=\t" << accu(abs(u01_check - u01)) << std::endl;
    Rcpp::Rcout << "Error on u10\t=\t" << accu(abs(u10_check - u10)) << std::endl;
    Rcpp::Rcout << "Error on u11\t=\t" << accu(abs(u11_check - u11)) << std::endl;
    Rcpp::Rcout << "Error on prior value\t\t=\t" << (prior_value - prior_value_check) << std::endl;
    Rcpp::Rcout << "Error on likelihood value\t=\t" << std::abs(likelihood_value - likelihood_value_check) << std::endl;
    Rcpp::Rcout << std::endl << std::endl;
  }
}



