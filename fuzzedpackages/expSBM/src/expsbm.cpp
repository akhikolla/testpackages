#include "expsbm.h"

expsbm::expsbm(unsigned int N_, arma::mat edgelist_, arma::mat Z_, arma::vec lambda_, arma::mat mu_, arma::mat nu_, bool directed_, bool trunc_, double tol_, unsigned int n_iter_max_, bool verbose_)
{
  N = N_;
  edgelist = edgelist_;
  lambda = lambda_;
  Z = Z_;
  K = Z.n_cols;
  mu = mu_;
  nu = nu_;
  directed = directed_;
  trunc = trunc_;
  verbose = verbose_;
  tol = tol_;
  n_iter_max = n_iter_max_;
  ConstructAdjacency();
  EvaluateDataSummaries();
  EvaluateStatistics();
  EvaluateELBO();
}

void expsbm::ConstructAdjacency()
{
  W.zeros(N,N);
  for (unsigned int l=0; l<edgelist.n_rows; ++l) 
  {
    W.at(edgelist.at(l,0),edgelist.at(l,1))++;
    if (!directed) W.at(edgelist.at(l,1),edgelist.at(l,0))++;
  }
  A.set_size(N,N);
  X.set_size(N,N);
  for (unsigned int i=0; i<N; ++i) for (unsigned int j=0; j<N; ++j) 
  {
    A.at(i,j).zeros(W.at(i,j));
    X.at(i,j).zeros(W.at(i,j));
  }
  arma::mat W_temp;
  W_temp.zeros(N,N);
  for (unsigned int l=0; l<edgelist.n_rows; ++l) 
  {
    unsigned int i = edgelist.at(l,0);
    unsigned int j = edgelist.at(l,1);
    A.at(i,j).at(W_temp.at(i,j)) = edgelist.at(l,2);
    X.at(i,j).at(W_temp.at(i,j)) = edgelist.at(l,3);
    W_temp.at(i,j)++;
    if (!directed)
    {
      A.at(j,i).at(W_temp.at(j,i)) = edgelist.at(l,2);
      X.at(j,i).at(W_temp.at(j,i)) = edgelist.at(l,3);
      W_temp.at(j,i)++;
    }
  }
}

void expsbm::EvaluateDataSummaries()
{
  A1.zeros(N,N);
  A0.zeros(N,N);
  for (unsigned int i=0; i<N; ++i) for (unsigned int j=0; j<N; ++j) if (i != j)
  {
    if (W.at(i,j)>2) for (unsigned int w=1; w<W.at(i,j)-1; ++w)// cycle through embedded segments
    {
      A1.at(i,j) += A.at(i,j).at(w);
      A0.at(i,j) += 1 - A.at(i,j).at(w);
    }
    if (!trunc)// if the first and last segments are not truncated, they contribute to the sum
    {
      A1.at(i,j) += A.at(i,j).at(0);
      A0.at(i,j) += 1 - A.at(i,j).at(0);
      if (W.at(i,j) > 1)// if an edge takes only one value, it should not contribute twice to the sum, hence this if statement
      {
        A1.at(i,j) += A.at(i,j).at(W.at(i,j)-1);
        A0.at(i,j) += 1 - A.at(i,j).at(W.at(i,j)-1);
      }
    }
  }
  X1.zeros(N,N);
  X0.zeros(N,N);
  for (unsigned int i=0; i<N; ++i) for (unsigned int j=0; j<N; ++j) if (i != j) for (unsigned int w=0; w<W.at(i,j); ++w)
  {
    X1.at(i,j) += X.at(i,j).at(w) * A.at(i,j).at(w);
    X0.at(i,j) += X.at(i,j).at(w) * (1-A.at(i,j).at(w));
  }
}

void expsbm::EvaluateStatistics()
{
  L_mu.zeros(K,K);
  L_nu.zeros(K,K);
  for (unsigned int g=0; g<K; ++g) for (unsigned int h=0; h<K; ++h) for (unsigned int i=0; i<N; ++i) for (unsigned int j=0; j<N; ++j) if (i != j)
  {
    L_mu.at(g,h) += A1.at(i,j) * Z.at(i,g) * Z.at(j,h);
    L_nu.at(g,h) += A0.at(i,j) * Z.at(i,g) * Z.at(j,h);
  }
  eta.zeros(K,K);
  zeta.zeros(K,K);
  for (unsigned int g=0; g<K; ++g) for (unsigned int h=0; h<K; ++h) for (unsigned int i=0; i<N; ++i) for (unsigned int j=0; j<N; ++j) if (i != j)
  {
    eta.at(g,h) += X1.at(i,j) * Z.at(i,g) * Z.at(j,h);
    zeta.at(g,h) += X0.at(i,j) * Z.at(i,g) * Z.at(j,h);
  }
  if (!directed)
  {
    L_mu /= 2;
    L_nu /= 2;
    eta /= 2;
    zeta /= 2;
  }
}

void expsbm::EvaluateELBO()
{
  elbo_value = 0;
  for (unsigned int g=0; g<K; ++g) for (unsigned int h=0; h<K; ++h) 
  {
    double contribution = L_mu.at(g,h) * log(mu.at(g,h)) + L_nu.at(g,h) * log(nu.at(g,h)) - mu.at(g,h) * eta.at(g,h) - nu.at(g,h) * zeta.at(g,h);
    elbo_value += contribution;
    if (!directed) if (g > h) elbo_value -= contribution;
  }
}


