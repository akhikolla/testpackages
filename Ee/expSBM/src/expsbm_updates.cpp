#include "expsbm.h"

void expsbm::UpdateZ(unsigned int l)
{
  arma::vec derivative_term = log(lambda);
  for (unsigned int k=0; k<K; ++k)
  {
    for (unsigned int h=0; h<K; ++h) for (unsigned int j=0; j<N; ++j) if (j != l)
    {
      derivative_term.at(k) += log(mu.at(k,h)) * Z.at(j,h) * A1.at(l,j);
      derivative_term.at(k) += log(nu.at(k,h)) * Z.at(j,h) * A0.at(l,j);
      derivative_term.at(k) += - mu.at(k,h) * Z.at(j,h) * X1.at(l,j);
      derivative_term.at(k) += - nu.at(k,h) * Z.at(j,h) * X0.at(l,j);
    }
    if (directed) for (unsigned int g=0; g<K; ++g) for (unsigned int i=0; i<N; ++i) if (i != l)
    {
      derivative_term.at(k) += log(mu.at(g,k)) * Z.at(i,g) * A1.at(i,l);
      derivative_term.at(k) += log(nu.at(g,k)) * Z.at(i,g) * A0.at(i,l);
      derivative_term.at(k) += - mu.at(g,k) * Z.at(i,g) * X1.at(i,l);
      derivative_term.at(k) += - nu.at(g,k) * Z.at(i,g) * X0.at(i,l);
    }
  }
  double prop_const = derivative_term.max();
  arma::vec transformed_derivative_term = derivative_term;
  for (unsigned int k=0; k<K; ++k) transformed_derivative_term.at(k) = exp(derivative_term.at(k) - prop_const);
  double transformed_derivative_term_sum = sum(transformed_derivative_term);
  for (unsigned int k=0; k<K; ++k) Z.at(l,k) = transformed_derivative_term.at(k) / transformed_derivative_term_sum;
}

void expsbm::UpdateLambda()
{
  arma::rowvec colsums = sum(Z,0);
  lambda = colsums.t() / accu(colsums);
}

void expsbm::UpdateMu(unsigned int g, unsigned int h)
{
  if (L_mu.at(g,h) > 0) 
  {
    mu.at(g,h) = L_mu.at(g,h) / eta.at(g,h);
    if (!directed) mu.at(h,g) = mu.at(g,h);
  }
}

void expsbm::UpdateNu(unsigned int g, unsigned int h)
{
  if (L_nu.at(g,h) > 0) 
  {
    nu.at(g,h) = L_nu.at(g,h) / zeta.at(g,h);
    if (!directed) nu.at(h,g) = nu.at(g,h);
  }
}
