#include "core_binder.h"

binder::binder(arma::mat sample_, arma::vec weights_, arma::vec decision_) {
  loss_function_name = "Binder's";
  sample = sample_;
  weights = weights_;
  decision = decision_;
  niter = sample.n_rows;
  sum_of_weights = sum(weights);
  N = sample.n_cols;
  Kup = sample.max() + 1;
  if (decision.max() + 1 > Kup) Kup = decision.max() + 1;
  EvaluateCounts();
  EvaluateLosses();
}

void binder::EvaluateLosses()
{
  unsigned int iter;
  arma::vec::iterator g, h;
  losses.set_size(niter);
  losses.fill(0);
  for (iter=0; iter<niter; ++iter)
  {
    for (g=non_empty_groups_decision.begin(); g<non_empty_groups_decision.end(); ++g)
    {
      losses.at(iter) += std::pow(decision_counts.at(*g),2) / 2;
    }
    for (h=non_empty_groups_sample.at(iter).begin(); h<non_empty_groups_sample.at(iter).end(); ++h)
    {
      losses.at(iter) += std::pow(sample_counts.at(iter,*h),2) / 2;
    }
    for (g=non_empty_groups_decision.begin(); g<non_empty_groups_decision.end(); ++g)
    {
      for (h=non_empty_groups_sample.at(iter).begin(); h<non_empty_groups_sample.at(iter).end(); ++h)
      {
        if (contingency_tables.at(*g,*h,iter) > 0) losses.at(iter) -= std::pow(contingency_tables.at(*g,*h,iter),2);
      }
    }
  }
  epl_value = arma::as_scalar(losses.t() * weights);
  epl_value /= sum_of_weights;
}

double binder::EvaluateDelta(unsigned int i, unsigned int h)
{
  unsigned int iter, g;
  g = decision.at(i);
  double res;
  res = 0;
  if (g != h)
  {
    res += decision_counts.at(h) - decision_counts.at(g) - 1;
    for (iter=0; iter<niter; ++iter)
    {
      res += 2 * (contingency_tables.at(g,sample.at(iter,i),iter) - contingency_tables.at(h,sample.at(iter,i),iter)) * weights.at(iter) / sum_of_weights;
    }
  }
  return (res);
}

void binder::EvaluateDeltas(unsigned int i)
{
  unsigned int h;
  deltas.set_size(Kup);
  deltas.fill(0);
  for (h=0; h<Kup; ++h) if (decision_counts.at(h) > 0) deltas.at(h) += EvaluateDelta(i,h);
  unsigned int first_empty_group;
  first_empty_group = Kup;
  for (h=0; h<Kup; ++h) if (decision_counts.at(h) == 0)
  {
    first_empty_group = h;
    break;
  }
  if (first_empty_group < Kup) deltas.at(h) += EvaluateDelta(i,first_empty_group);
}

void binder::Move(unsigned int i, unsigned int h)
{
  unsigned int iter, g, r, s, K;
  g = decision.at(i);
  if (g != h)
  {
    decision.at(i) = h;
    decision_counts.at(g)--;
    decision_counts.at(h)++;
    if (decision.at(g) == 0 || decision_counts.at(h) == 1)
    {
      K = 0;
      for (r=0; r<Kup; ++r) if (decision_counts.at(r) > 0) K++;
      non_empty_groups_decision.set_size(K);
      s = 0;
      for (r=0; r<Kup; ++r) if (decision_counts.at(r) > 0)
      {
        non_empty_groups_decision.at(s) = r;
        ++s;
      }
    }
    for (iter=0; iter<niter; ++iter)
    {
      contingency_tables.at(g,sample.at(iter,i),iter)--;
      contingency_tables.at(h,sample.at(iter,i),iter)++;
    }
    epl_value += deltas.at(h);
  }
}

