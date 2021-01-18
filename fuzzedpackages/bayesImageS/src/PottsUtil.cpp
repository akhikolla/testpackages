// ----------------------------------------------------------------------
// This file is part of the R package bayesImageS. It contains utility
// functions common to both Metropolis-Hastings and sequential Monte
// Carlo algorithms for the hidden Potts model.
// Copyright (C) 2013-2015  Matthew Moores
//
// bayesImageS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bayesImageS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------
#include "PottsUtil.h"

arma::uvec unsign(const Rcpp::IntegerVector & x)
{
  arma::uvec result(x.size());
  for (unsigned i=0; i<result.size(); i++)
  {
    result[i] = (unsigned) x[i];
  }
  return result;
}

arma::umat unsignMx(const Rcpp::IntegerMatrix & m)
{
  arma::umat result(m.nrow(), m.ncol());
  for (unsigned i=0; i<result.n_rows; i++)
  {
    for (unsigned j=0; j<result.n_cols; j++)
    {
      result(i,j) = (unsigned) m(i,j);
    }
  }
  return result;
}

arma::rowvec rgamma(const arma::rowvec & shape, const arma::rowvec & rate)
{
  arma::rowvec result(shape.n_elem);
  for (unsigned i=0; i<shape.n_elem; i++)
  {
    result[i] = ::Rf_rgamma(shape[i], 1 / rate[i]);
  }
  return result;
}

arma::rowvec rnorm(const arma::rowvec & mean, const arma::rowvec & stddev)
{
  Rcpp::NumericVector xR = Rcpp::rnorm(mean.n_elem);
  arma::rowvec x(xR.begin(), xR.size(), false);
  return x % stddev + mean; // element-wise multiplication and addition
}

arma::mat dnorm(const Rcpp::NumericVector & yunique, const arma::uvec & ymatch,
                const arma::rowvec & mean, const arma::rowvec & stddev)
{
  arma::vec prob(yunique.size());
  arma::mat probMx(ymatch.n_elem,mean.n_elem);
  for (unsigned i=0; i<mean.n_elem; i++)
  {
    for (int j=0; j<yunique.size(); j++)
    {
      prob[j] = ::Rf_dnorm4(yunique[j],mean[i],stddev[i],1);
    }
    probMx.col(i) = prob.elem(ymatch);
  }
  return probMx;
}

arma::umat randomIndices(const unsigned n, int k)
{
  Rcpp::NumericVector xR = Rcpp::runif(n, 0, k);
  arma::umat indices = arma::zeros<arma::umat>(n+1,k);
#pragma omp parallel for shared(indices)
  for (unsigned i=0; i<n; i++)
  {
    unsigned j = (unsigned)xR[i];
    indices(i,j) = 1;
  }
  return indices;
}

// the sufficient statistic of the Potts model: the number of identical pairs of neighbours
unsigned sum_ident(const arma::umat & z, const arma::umat & neigh, const std::vector<arma::uvec> & blocks)
{
  unsigned total = 0;
  const arma::uvec block = blocks[0];
#pragma omp parallel for reduction(+:total)
  for (unsigned i=0; i < block.n_elem; i++)
  {    
    for (unsigned j=0; j < z.n_cols; j++)
    {
      if (z(block(i),j) == 1)
      {
        unsigned sum_neigh = 0;
        for (unsigned k=0; k < neigh.n_cols; k++)
        {
          sum_neigh += z(neigh(block(i),k),j);
        }
        total += sum_neigh;
      }
    }
  }
  return total;
}

// the log sum of a vector of logs
// http://jblevins.org/log/log-sum-exp
double sum_logs(arma::vec log_prob)
{
  double suml = 0.0;
  double maxl = log_prob.max();
  for (unsigned i=0; i < log_prob.n_elem; i++)
  {
    if (arma::is_finite(log_prob(i)))
      suml += exp(log_prob(i) - maxl);
  }
  return log(suml) + maxl;
}

arma::rowvec gibbsMeans(const arma::rowvec & nZ, const arma::rowvec & sumY,
                        const arma::rowvec & pr_mu, const arma::rowvec & pr_mu_tau,
                        const arma::rowvec & sigma)
{
  arma::rowvec oldTau = arma::pow(sigma, -2);
  arma::rowvec newTau = pr_mu_tau + oldTau % nZ;
  arma::rowvec mean = (pr_mu_tau%pr_mu + oldTau%sumY) / newTau; // element-wise
  return rnorm(mean, arma::pow(newTau, -0.5));
}

arma::rowvec gibbsStdDev(const arma::rowvec & nZ, const arma::rowvec & sumY,
                         const arma::rowvec & sqDiff, const arma::rowvec & pr_sd_nu,
                         const arma::rowvec & pr_sd_SS, const arma::rowvec & mean)
{
  // avoid dividing by zero if one of the mixture components is empty
  arma::rowvec Ybar(sumY.n_elem);
  for (unsigned j=0; j < sumY.n_elem; j++)
  {
    if (nZ[j] == 0) Ybar[j] = 0;
    else Ybar[j] = sumY[j] / nZ[j];
  }
  arma::rowvec shape = (pr_sd_nu + nZ)/2;
  arma::rowvec rate = (pr_sd_SS + sqDiff + nZ%arma::square(Ybar - mean))/2;
  return arma::pow(rgamma(shape, rate), -0.5);
}

arma::rowvec gibbsDirichlet(const arma::rowvec & nZ, const arma::rowvec & pr_lambda)
{
  arma::rowvec shape = nZ + pr_lambda;
  arma::rowvec sim = rgamma(shape, arma::ones<arma::rowvec>(shape.n_elem));
  return sim/arma::sum(sim);
}

// updates labels Z and count of allocations alloc
void gibbsLabels(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                 arma::umat & z, arma::umat & alloc, const double beta,
                 const arma::mat & log_xfield)
{
  const Rcpp::NumericVector randU = Rcpp::runif(neigh.n_rows);

  // the blocks are conditionally independent
  for (unsigned b=0; b < blocks.size(); b++)
  {
    const arma::uvec block = blocks[b];
    // for each pixel in the block
#pragma omp parallel for
    for (unsigned i=0; i < block.size(); i++)
    {
      // compute posterior probability for each label j
      arma::vec log_prob(z.n_cols);
      for (unsigned j=0; j < z.n_cols; j++)
      {
        unsigned sum_neigh = 0;
        for (unsigned k=0; k < neigh.n_cols; k++)
        {
          sum_neigh += z(neigh(block[i],k),j);
        }
        log_prob[j] = log_xfield(block[i],j) + beta*sum_neigh;
      }
      double total_llike = sum_logs(log_prob);

      // update labels Z
      double cumProb = 0.0;
      z.row(block[i]).zeros();
      for (unsigned j=0; j < log_prob.n_elem; j++)
      {
        cumProb += exp(log_prob[j] - total_llike);
        if (randU[block[i]] < cumProb)
        {
          z(block[i],j) = 1;
          alloc(block[i],j) += 1;
          break;
        }
      }
    }
  }
}

void gibbsLabelsNoData(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                 arma::umat & z, arma::umat & alloc, const double beta)
{
  const Rcpp::NumericVector randU = Rcpp::runif(neigh.n_rows);

  // the blocks are conditionally independent
  for (unsigned b=0; b < blocks.size(); b++)
  {
    const arma::uvec block = blocks[b];
    // for each pixel in the block
#pragma omp parallel for
    for (unsigned i=0; i < block.n_elem; i++)
    {
      // compute posterior probability for each label j
      arma::vec log_prob(z.n_cols);
      for (unsigned j=0; j < z.n_cols; j++)
      {
        unsigned sum_neigh = 0;
        for (unsigned k=0; k < neigh.n_cols; k++)
        {
          sum_neigh += z(neigh(block[i],k),j);
        }
        log_prob[j] = beta*sum_neigh;
      }
      double total_llike = sum_logs(log_prob);

      // update labels Z
      double cumProb = 0.0;
      z.row(block[i]).zeros();
      for (unsigned j=0; j < log_prob.n_elem; j++)
      {
        cumProb += exp(log_prob[j] - total_llike);
        if (randU[block[i]] < cumProb)
        {
          z(block[i],j) = 1;
          alloc(block[i],j) += 1;
          break;
        }
      }
    }
  }
}

void swLabelsNoData(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                    const double beta, const double k, arma::umat & z, arma::umat & alloc)
{
  // the neighbourhood relation is symmetrical, so we only need one block
  const arma::uvec block = blocks[0];
  const Rcpp::NumericVector randU = Rcpp::runif(block.n_elem * neigh.n_cols);
  arma::umat bonds(block.n_elem * neigh.n_cols, 2);
  unsigned nbonds = 0;
  for (unsigned i=0; i < block.n_elem; i++)
  {
    for (unsigned n=0; n < neigh.n_cols; n++)
    {
      unsigned j = neigh(block(i),n);
      if (j <= neigh.n_rows)
      {
        double eq = arma::as_scalar(z.row(block(i)) * z.row(j).t()); // either 0 or 1
        double p = 1 - exp(-beta * eq);
        if (randU(i * neigh.n_cols + n) < p)
        {
          // add a bond
          bonds(nbonds,0) = block(i);
          bonds(nbonds,1) = j;
          nbonds++;
        }
      }
    }
  }
  
  // coalesce bonds into clusters
  arma::uvec clust(alloc.n_rows);
  for (unsigned c=0; c<alloc.n_rows; c++) clust(c) = c;
  for (unsigned b=0; b<nbonds; b++)
  {
    unsigned p0 = bonds(b,0);
    unsigned q0 = bonds(b,1);
    unsigned p1 = clust(p0);
    unsigned q1 = clust(q0);
    while(p1 != q1)
    {
      if(q1 < p1)
      {
        clust(p0) = q1;
        p0 = p1;
      	p1 = clust(p1);
      }
      else
      {
      	clust(q0) = p1;
      	q0 = q1;
      	q1 = clust(q1);
      }
    }
  }
  for(unsigned i = 0; i < alloc.n_rows; i++)
  {
    clust(i) = clust(clust(i));
  }

  arma::uvec uniq_clust = arma::unique(clust);
  Rcpp::NumericVector randJ = Rcpp::runif(uniq_clust.n_elem);
  arma::uvec newJ(alloc.n_rows);
  for (unsigned c=0; c<uniq_clust.n_elem; c++)
  {
    for (double j=0; j<k; j++)
    {
      if (randJ[c] < (j+1.0)/k)
      {
        newJ(uniq_clust(c)) = (unsigned)j;
        break;
      }
    }
  }

  // flip all cluster pixels to the same value
  for (unsigned i=0; i<alloc.n_rows; i++)
  {
    z.row(i).zeros();
    z(i,newJ(clust(i))) = 1;
    alloc(i,newJ(clust(i))) += 1;
  }
}

// Computes the number of neighbouring pixels allocated to component j, for pixel i. 
void neighbj(arma::mat & ne, arma::uvec & e, const arma::umat & z, const arma::umat & neigh)
{
#pragma omp parallel for
  for (unsigned i=0; i < z.n_rows-1; i++)
  {
    for (unsigned j=0; j < z.n_cols; j++)
    {
      unsigned sum_neigh = 0;
      for (unsigned k=0; k < neigh.n_cols; k++)
      {
        sum_neigh += z(neigh(i,k),j);
      }
      ne(j,i) = (double)sum_neigh;
      if (z(i,j) == 1)
      {
        e[i] = j;
      }
    }
  }
}

// pseudo log-likelihood of the Potts model
double pseudolike(const arma::mat & ne, const arma::uvec & e, const double b, const unsigned n, const unsigned k)
{
  double num = 0.0;
  double denom = 0.0;
#pragma omp parallel for reduction(+:num,denom)
  for (unsigned i=0; i < n; i++)
  {
    num=num+ne(e[i],i);
    double tdenom=0.0;
    for (unsigned j=0; j < k; j++)
    {
      tdenom=tdenom+exp(b*ne(j,i));
    }
    denom=denom+log(tdenom);
  }
  return b*num-denom;
}

// updates labels Z and count of allocations alloc
void classify(arma::umat & z, arma::umat & alloc, const arma::rowvec & lambda, const arma::mat & log_xfield)
{
  const Rcpp::NumericVector randU = Rcpp::runif(log_xfield.n_rows);

#pragma omp parallel for
  for (unsigned i=0; i < log_xfield.n_rows; i++)
  {
    // compute posterior probability for each label j
    arma::vec log_prob(z.n_cols);
    for (unsigned j=0; j < z.n_cols; j++)
    {
      log_prob(j) = log_xfield(i,j) + lambda(j);
    }
    double total_llike = sum_logs(log_prob);

    // update labels Z
    double cumProb = 0.0;
    z.row(i).zeros();
    for (unsigned j=0; j < log_prob.n_elem; j++)
    {
      cumProb += exp(log_prob[j] - total_llike);
      if (randU[i] < cumProb)
      {
        z(i,j) = 1;
        alloc(i,j) += 1;
        break;
      }
    }
  }
}

void updateStats(const arma::colvec & y, const arma::umat & z,
                 arma::rowvec & nZ, arma::rowvec & sumY, arma::rowvec & sqDiff)
{
  nZ.zeros();
  sumY.zeros();
  sqDiff.zeros();
  for (unsigned i=0; i < y.n_elem; i++)
  {
    for (unsigned j=0; j < z.n_cols; j++)
    {
      if (z(i,j)==1)
      {
        nZ[j]++;
        sumY[j] += y[i];
      }
    }
  }
  arma::rowvec ybar = sumY/nZ;
  for (unsigned i=0; i < y.n_elem; i++)
  {
    for (unsigned j=0; j < z.n_cols; j++)
    {
      if (z(i,j)==1)
      {
        sqDiff[j] += pow(y[i] - ybar[j],2);
      }
    }
  }
}

double rwmh(const double mean, const double stddev, const double prior[2])
{
  double proposal = ::Rf_rnorm(mean, stddev);
  if (proposal < prior[0])
  {
    proposal = std::min(prior[1], prior[0] + prior[0] - proposal);
  }
  if (proposal > prior[1])
  {
    proposal = std::max(prior[0], prior[1] + prior[1] - proposal);
  }
  return proposal;
}

// linear interpolation
double interp(double val, unsigned idx, const arma::mat & path)
{
  return(path(1,idx) + (val - path(0,idx))*(path(1,idx+1) - path(1,idx))/(path(0,idx+1)-path(0,idx)));
}
