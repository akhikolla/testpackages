# include "core_sample_of_partitions.h"

sample_of_partitions::sample_of_partitions(arma::mat sample_, arma::vec weights_, arma::vec decision_) {
  loss_function_name = "None";
  sample = sample_;
  weights = weights_;
  decision = decision_;
  niter = sample.n_rows;
  sum_of_weights = accu(weights);
  N = sample.n_cols;
  Kup = sample.max() + 1;
  if (decision.max() + 1 > Kup) Kup = decision.max() + 1;
  EvaluateCounts();
  EvaluateLosses();
}

/*
void sample_of_partitions::Print() {
  unsigned int iter, i, g;
  std::ostringstream strs;
  strs << "Class: sample_of_partitions\n";
  strs << "\nDimensions:\n" << "niter\t=\t" << niter << "\nN\t=\t" << N << "\nKup\t=\t" << Kup << "\n";
  strs << "\nPosterior sample:\n";
  sample.print(strs);
  strs << "\n\nSample weights:\n";
  weights.t().print(strs);
  strs << "\n\nCounts:\n";
  sample_counts.print(strs);
  strs << "\n\nNon-empty groups for the sample:\n";
  for (iter=0; iter<niter; ++iter) non_empty_groups_sample.at(iter).t().print(strs);
  strs << "\n\n";
  strs << "\nLoss function:\t" << loss_function_name << "\n";
  strs << "\nDecision:\n";
  for (i=0; i<N; ++i) strs << decision.at(i) << " ";
  strs << "\n";
  strs << "\nCounts:\n";
  for (g=0; g<Kup; ++g) strs << decision_counts.at(g) << " ";
  strs << "\n";
  strs << "\nNon-empty groups for the decision:\n";
  for (g=0; g<non_empty_groups_decision.size(); ++g) strs << non_empty_groups_decision.at(g) << " ";
  strs << "\n";
  strs << "\nLosses:\n";
  losses.t().print(strs);
  strs << "\nExpected Posterior Loss\t=\t" << epl_value;
  strs << "\n\n\n";
  std::cout << strs.str() << std::endl;
}

void sample_of_partitions::Summary()
{
  unsigned int i;
  std::ostringstream strs;
  strs << "Class: SampleOfPartitions\n";
  strs << "\nDimensions:\n" << "niter\t=\t" << niter << "\nN\t=\t" << N << "\nKup\t=\t" << Kup << "\n";
  strs << "\nDecision:\n";
  for (i=0; i<N; ++i) strs << decision.at(i) << " ";
  strs << "\n";
  strs << "\n\nNumber of non-empty groups\t=\t" << non_empty_groups_decision.size();
  strs << "\n\nExpected Posterior Loss\t=\t" << epl_value;
  strs << "\n\n\n";
  std::cout << strs.str() << std::endl;
}
*/

void sample_of_partitions::EvaluateCounts()
{
  unsigned int iter, i, g, h, K;
  sample_counts.set_size(niter,Kup);
  sample_counts.fill(0);
  decision_counts.set_size(Kup);
  decision_counts.fill(0);
  contingency_tables.set_size(Kup,Kup,niter);
  contingency_tables.fill(0);
  for (i=0; i<N; ++i)
  {
    decision_counts.at(decision.at(i))++;
    for (iter=0; iter<niter; ++iter)
    {
      sample_counts.at(iter,sample.at(iter,i))++;
      contingency_tables.at(decision.at(i),sample.at(iter,i),iter)++;
    }
  }
  non_empty_groups_sample.set_size(niter);
  for (iter=0; iter<niter; ++iter)
  {
    K = 0;
    for (g=0; g<Kup; ++g) if (sample_counts.at(iter,g) > 0) K++;
    non_empty_groups_sample.at(iter).set_size(K);
    h = 0;
    for (g=0; g<Kup; ++g)
    {
      if (sample_counts.at(iter,g) > 0)
      {
        non_empty_groups_sample.at(iter).at(h) = g;
        ++h;
      }
    }
  }
  K = 0;
  for (g=0; g<Kup; ++g) if (decision_counts.at(g) > 0) K++;
  non_empty_groups_decision.set_size(K);
  h = 0;
  for (g=0; g<Kup; ++g) if (decision_counts.at(g) > 0)
  {
    non_empty_groups_decision.at(h) = g;
    ++h;
  }
}

void sample_of_partitions::EvaluateLosses() {
  arma::vec::iterator g, h;
  losses.set_size(niter);
  losses.fill(0);
  epl_value = arma::as_scalar(losses.t() * weights);
  epl_value /= sum_of_weights;
}

double sample_of_partitions::EvaluateDelta(unsigned int i, unsigned int h)
{
  double res;
  res = 0;
  return (res);
}

void sample_of_partitions::EvaluateDeltas(unsigned int i)
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

void sample_of_partitions::Move(unsigned int i, unsigned int h)
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
