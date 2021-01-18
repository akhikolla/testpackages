#include "core_variation_of_information.h"

variation_of_information::variation_of_information(arma::mat sample_, arma::vec weights_, arma::vec decision_) {
  loss_function_name = "Variation of Information";
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

void variation_of_information::EvaluateLosses()
{
  unsigned int iter;
  arma::vec::iterator g, h;
  losses.set_size(niter);
  losses.fill(0);
  for (iter=0; iter<niter; ++iter)
  {
    for (g=non_empty_groups_decision.begin(); g<non_empty_groups_decision.end(); ++g)
    {
      losses.at(iter) += Entropy(decision_counts.at(*g)) / N;
    }
    for (h=non_empty_groups_sample.at(iter).begin(); h<non_empty_groups_sample.at(iter).end(); ++h)
    {
      losses.at(iter) += Entropy(sample_counts.at(iter,*h)) / N;
    }
    for (g=non_empty_groups_decision.begin(); g<non_empty_groups_decision.end(); ++g)
    {
      for (h=non_empty_groups_sample.at(iter).begin(); h<non_empty_groups_sample.at(iter).end(); ++h)
      {
        losses.at(iter) -= 2 * Entropy(contingency_tables.at(*g,*h,iter)) / N;
      }
    }
  }
  epl_value = arma::as_scalar(losses.t() * weights);
  epl_value /= sum_of_weights;
}

double variation_of_information::EvaluateDelta(unsigned int i, unsigned int h)
{
  unsigned int iter, g;
  g = decision.at(i);
  double res;
  res = 0;
  if (g != h)
  {
    res += Entropy(decision_counts.at(g)-1) - Entropy(decision_counts.at(g));
    res += Entropy(decision_counts.at(h)+1) - Entropy(decision_counts.at(h));
    res /= N;
    for (iter=0; iter<niter; ++iter)
    {
      res += (   -2*Entropy(contingency_tables.at(g,sample.at(iter,i),iter)-1) +
        2*Entropy(contingency_tables.at(g,sample.at(iter,i),iter)) -
        2*Entropy(contingency_tables.at(h,sample.at(iter,i),iter)+1) +
        2*Entropy(contingency_tables.at(h,sample.at(iter,i),iter))      )*weights.at(iter)/sum_of_weights/N;
    }
  }
  return (res);
}

void variation_of_information::EvaluateDeltas(unsigned int i)
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

void variation_of_information::Move(unsigned int i, unsigned int h)
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
