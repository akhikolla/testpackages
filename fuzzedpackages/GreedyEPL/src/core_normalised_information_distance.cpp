#include "core_normalised_information_distance.h"

normalised_information_distance::normalised_information_distance(arma::mat sample_, arma::vec weights_, arma::vec decision_) {
  loss_function_name = "Normalised Information Distance";
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

void normalised_information_distance::EvaluateLosses()
{
  unsigned int iter;
  arma::vec::iterator g, h;
  entropy_decision = 0;
  for (g=non_empty_groups_decision.begin(); g<non_empty_groups_decision.end(); ++g) entropy_decision -= Entropy(decision_counts.at(*g)/ N);
  entropies_sample.zeros(niter);
  for (iter=0; iter<niter; ++iter) for (h=non_empty_groups_sample.at(iter).begin(); h<non_empty_groups_sample.at(iter).end(); ++h)
  {
    entropies_sample.at(iter) -= Entropy(sample_counts.at(iter,*h)/N);
  }
  joint_entropies.zeros(niter);
  for (iter=0; iter<niter; ++iter) for (g=non_empty_groups_decision.begin(); g<non_empty_groups_decision.end(); ++g) for (h=non_empty_groups_sample.at(iter).begin(); h<non_empty_groups_sample.at(iter).end(); ++h)
  {
    joint_entropies.at(iter) -= Entropy(contingency_tables.at(*g,*h,iter)/N);
  }
  losses.zeros(niter);
  for (iter=0; iter<niter; ++iter) losses.at(iter) += 1 - ( entropy_decision + entropies_sample.at(iter) - joint_entropies.at(iter) ) / std::max(entropy_decision, entropies_sample.at(iter));
  epl_value = arma::as_scalar(losses.t() * weights);
  epl_value /= sum_of_weights;
}

double normalised_information_distance::EvaluateDelta(unsigned int i, unsigned int h)
{
  unsigned int iter, g;
  g = decision.at(i);
  double res;
  res = 0;
  double new_entropy_decision, new_joint_entropy;
  if (g != h)
  {
    new_entropy_decision = entropy_decision - Entropy((decision_counts.at(g)-1)/N) + Entropy(decision_counts.at(g)/N) - Entropy((decision_counts.at(h)+1)/N) + Entropy((decision_counts.at(h))/N);
    for (iter=0; iter<niter; ++iter)
    {
      new_joint_entropy = joint_entropies.at(iter) - Entropy((contingency_tables.at(g,sample.at(iter,i),iter)-1)/N) + Entropy((contingency_tables.at(g,sample.at(iter,i),iter))/N) - Entropy((contingency_tables.at(h,sample.at(iter,i),iter)+1)/N) + Entropy((contingency_tables.at(h,sample.at(iter,i),iter))/N);
      res += (1 - ( new_entropy_decision + entropies_sample.at(iter) - new_joint_entropy ) / std::max(new_entropy_decision, entropies_sample.at(iter)) - losses.at(iter)) * weights.at(iter);
    }
    res /= sum_of_weights;
  }
  return (res);
}

void normalised_information_distance::EvaluateDeltas(unsigned int i)
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

void normalised_information_distance::Move(unsigned int i, unsigned int h)
{
  unsigned int iter, g, r, s, K;
  g = decision.at(i);
  if (g != h)
  {
    entropy_decision -= Entropy((decision_counts.at(g)-1)/N) - Entropy(decision_counts.at(g)/N) + Entropy((decision_counts.at(h)+1)/N) - Entropy(decision_counts.at(h)/N);
    for (iter=0; iter<niter; ++iter)
    {
      joint_entropies.at(iter) -= Entropy((contingency_tables.at(g,sample.at(iter,i),iter)-1)/N) - Entropy(contingency_tables.at(g,sample.at(iter,i),iter)/N) + Entropy((contingency_tables.at(h,sample.at(iter,i),iter)+1)/N) - Entropy(contingency_tables.at(h,sample.at(iter,i),iter)/N);
      losses.at(iter) = 1 - ( entropy_decision + entropies_sample.at(iter) - joint_entropies.at(iter) ) / std::max(entropy_decision, entropies_sample.at(iter));
    }
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

