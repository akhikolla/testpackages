// ----------------------------------------------------------------------
// This file is part of the R package bayesImageS. It contains sequential
// Monte Carlo algorithms for image segmentation using a hidden Potts
// model with approximate Bayesian computation (SMC-ABC).
// Copyright (C) 2013-2017  Matthew Moores
//
// bayesImageS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
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
#include "smcPotts.h"
#include "PottsUtil.h"

Rcpp::NumericVector updatePseudoGibbs(const double particle, const unsigned nValues,
               const unsigned k, const arma::umat & neigh, const std::vector<arma::uvec> & blocks)
{
  const unsigned burn = 500;
  const unsigned niter = burn + nValues;
  Rcpp::NumericVector values(nValues);
  arma::umat z = randomIndices(neigh.n_rows, k);
  arma::umat alloc = arma::zeros<arma::umat>(neigh.n_rows, k);
  for (unsigned it=0; it<niter; it++){
    // update labels
    gibbsLabelsNoData(neigh, blocks, z, alloc, particle);
    if (it >= burn)
    {
      values[it-burn] = sum_ident(z, neigh, blocks);
    }
  }
  return(values);
}

Rcpp::NumericVector updatePseudoSW(const double particle, const unsigned nValues,
               const unsigned k, const arma::umat & neigh, const std::vector<arma::uvec> & blocks)
{
  const unsigned burn = 100;
  const unsigned niter = burn + nValues;
  Rcpp::NumericVector values(nValues);
  arma::umat z = randomIndices(neigh.n_rows, k);
  arma::umat alloc = arma::zeros<arma::umat>(neigh.n_rows, k);
  for (unsigned it=0; it<niter; it++){
    // update labels
    swLabelsNoData(neigh, blocks, particle, k, z, alloc);
    if (it >= burn)
    {
      values[it-burn] = sum_ident(z, neigh, blocks);
    }
  }
  return(values);
}

Rcpp::NumericVector updatePseudoPath(const double particle, const unsigned nValues,
               const unsigned k, const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
               arma::mat path, arma::mat sdMx)
{
  unsigned ix = 0;
  while(path(0,ix) <= particle) ix++;
  double szMu = interp(particle, ix-1, path);
  double szSd = interp(particle, ix-1, sdMx);
  return(Rcpp::rnorm(nValues, szMu, szSd));
}  

Rcpp::NumericVector updatePseudo(const double particle, const unsigned nValues,
               const unsigned k, const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
               arma::mat pathMx, arma::mat sdMx, bool aux_sw)
{
  if (pathMx.n_rows > 0)
  {
    return updatePseudoPath(particle, nValues, k, neigh, blocks, pathMx, sdMx);
  }
  else if (aux_sw)
  {
    return updatePseudoSW(particle, nValues, k, neigh, blocks);
  }
  else
  {
    return updatePseudoGibbs(particle, nValues, k, neigh, blocks);
  }
}

Rcpp::NumericMatrix updatePseudo(const Rcpp::NumericVector &particles, const unsigned nValues,
               const unsigned k, const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
               arma::mat pathMx, arma::mat sdMx, bool aux_sw)
{
  Rcpp::NumericMatrix values(particles.size(), nValues);
  for (int p=0; p<particles.size(); p++)
  {
    values(p, Rcpp::_) = updatePseudo(particles[p], nValues, k, neigh, blocks, pathMx, sdMx, aux_sw);
  }
  return values;
}

unsigned surv(const Rcpp::NumericVector &pseudo, const unsigned stat, const double epsilon)
{
  unsigned x=0;
  for (int m=0; m<pseudo.size(); m++)
  {
    if (fabs(pseudo[m] - stat) < epsilon)
    {
      x++;
    }
  }
  return(x);
}

Rcpp::NumericVector survivors(Rcpp::NumericMatrix &pseudo, const unsigned stat, const double epsilon)
{
  Rcpp::NumericVector prop(pseudo.nrow());
  for (int p=0; p<pseudo.nrow(); p++)
  {
    const Rcpp::NumericVector ps = pseudo(p, Rcpp::_);
    prop[p] = surv(ps, stat, epsilon);
  }
  return(prop);
}

/**
 * Computes the effective sample size (ESS) according to Liu (2001)
 */
double effectiveSampleSize(const arma::vec &log_wt)
{
  double sum_wt = sum_logs(log_wt);
  double sum_sq = sum_logs(log_wt + log_wt);
  double res = exp(sum_wt + sum_wt - sum_sq);
  if (arma::is_finite(res)) return res;
  else return 0;
}

arma::vec calcWeights(const arma::vec &old_wt, const Rcpp::NumericVector &new_surv, const Rcpp::NumericVector &old_surv)
{
  arma::vec wt(old_wt.size());
  for (unsigned p=0; p<wt.n_elem; p++)
  {
    wt[p] = old_wt[p] + log(new_surv[p]) - log(old_surv[p]);
  }
  // normalize to sum to 1
  double sum_wt = sum_logs(wt);
  return(wt - sum_wt);
}

/**
 * Choose a new ABC tolerance such that the effective sample size (ESS)
 * is reduced by a factor of approximately ess_delta, using bisection.
 */
arma::vec updateImportanceWeights(const arma::vec &old_wt,
                       std::vector<double> &epsilon, std::vector<double> &ess,
                       Rcpp::NumericMatrix &pseudo, const unsigned stat, const double alpha)
{
  double a=0.0, b=epsilon.back(), c;
  const double tol = 0.01 * alpha * (double)pseudo.nrow();
  unsigned max_it = 1000; // sensible upper limit on computational cost vs. precision of epsilon
  Rcpp::NumericVector surv_old = survivors(pseudo,stat,b);
  arma::vec wt_c(old_wt.size());
  double ess_b = ess.back(), ess_c;
  const double ess_target = alpha * ess_b;
  Rcpp::Rcout << "previous epsilon " << b << " and ESS " << ess_b << " (target: " << ess_target << ")\n";
  unsigned i=0;
  do
  {
    c = (a + b)/2.0;
    Rcpp::NumericVector surv_c = survivors(pseudo,stat,c);
    wt_c = calcWeights(old_wt, surv_c, surv_old);
    ess_c = effectiveSampleSize(wt_c);
    if (ess_c < ess_target) a = c; // lower bound
    else b = c;
  } while (++i<=max_it && fabs(ess_c - ess_target) > tol);
  epsilon.push_back(c);
  ess.push_back(ess_c);
  return(wt_c);
}

/**
 * Simple multinomial sampling with replacement
 */
Rcpp::NumericVector subsample(Rcpp::NumericVector &particles, arma::vec &log_wt, unsigned n)
{
  const Rcpp::NumericVector randU = Rcpp::runif(n);
  Rcpp::NumericVector result(n);
  for (unsigned i=0; i<n; i++)
  {
    // select a particle at random, according to the weights
    double total = 0.0;
    for (int j=0; j < particles.size() && total <= randU[i]; j++)
    {
      if (arma::is_finite(log_wt(j)))
      {
        total += exp(log_wt(j));
      }
      result[i] = particles[j];
    }
  }
  return(result);
}

/**
 * Residual resampling to decrease the variance (Douc, Cappe & Moulines 2005)
 */
Rcpp::IntegerVector resample_resid(Rcpp::NumericVector &particles, arma::vec &log_wt,
                                   Rcpp::NumericMatrix &pseudo)
{
  const unsigned n = particles.size();
  Rcpp::IntegerVector idx(n);

  // first loop is deterministic: only accept particles with n*weight > 1
  unsigned r=0;
  for (unsigned i=0; i<n; i++)
  {
    if (arma::is_finite(log_wt(i)))
    {
      int tW = (int)trunc(exp(log_wt(i) + log((double)n)));
      for (int j=0; j < tW; j++)
      {
        idx[r+j] = i;
      }
      r += tW;
      log_wt(i) = log(exp(log_wt(i) + log((double)n)) - tW);
    }
  }
  // renormalize the weights
  log_wt = log_wt - log((double)(n-r));

  // second loop uses multinomial resampling for the remaining n-r particles
  const Rcpp::NumericVector randU = Rcpp::runif(n-r);
  for (unsigned i=r; i<n; i++)
  {
    // select a particle at random, according to the weights
    double total = 0.0;
    for (unsigned j=0; j < n && total <= randU[i-r]; j++)
    {
      if (arma::is_finite(log_wt(j)))
      {
        total += exp(log_wt(j));
      }
      idx[i] = j;
    }
  }
  return(idx.sort());
}

/**
 * weighted arithmetic mean of the particles
 */
double weightedMean(const Rcpp::NumericVector &particles, const arma::vec &log_weights)
{
  double suml = 0.0;
  double maxl = log_weights.max();
  for (unsigned i=0; i < log_weights.n_elem; i++)
  {
    if (arma::is_finite(log_weights[i]))
      suml += exp(log_weights(i) - maxl + log(particles[i]));
  }
  return suml * exp(maxl);
}

/**
 * weighted sample variance of the particles
 */
double weightedVariance(const Rcpp::NumericVector &particles, const arma::vec &log_weights, double mean)
{
  double suml = 0.0;
  double maxl = log_weights.max();
  for (unsigned i=0; i < log_weights.n_elem; i++)
  {
    if (arma::is_finite(log_weights[i]))
      suml += exp(log_weights(i) - maxl) * pow(particles[i] - mean, 2.0);
  }
  return suml * exp(maxl);
}

/**
 * update particles using a random walk with the specified bandwidth
 */
unsigned rwmhParticles(Rcpp::NumericVector &particles, const arma::vec &log_wt,
              Rcpp::NumericMatrix &pseudo, double bw, const double prior[2],
              const unsigned stat, const double epsilon,
              const unsigned k, const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
              arma::mat pathMx, arma::mat sdMx, bool aux_sw)
{
  unsigned acc=0;
  for (int i=0; i<particles.size(); i++)
  {
    if (arma::is_finite(log_wt(i)))
    {
      double prop = rwmh(particles[i], bw, prior);
      Rcpp::NumericVector pS = updatePseudo(prop, pseudo.ncol(), k, neigh, blocks, pathMx, sdMx, aux_sw);
      Rcpp::NumericVector oS = pseudo(i, Rcpp::_);
      double log_ratio = log((double)surv(pS, stat, epsilon)) - log((double)surv(oS, stat, epsilon));
      // accept/reject
      if (unif_rand() < exp(log_ratio))
      {
        acc++;
        particles[i] = prop;
        for (int m=0; m<pseudo.ncol(); m++)
        {
          pseudo(i,m) = pS[m];
        }
      }
    }
  }
  return(acc);
}

SEXP smcPotts(SEXP yS, SEXP nS, SEXP bS, SEXP parS, SEXP prS)
{
  BEGIN_RCPP
  time_t timer,timer2;
  time(&timer);
  Rcpp::NumericVector yR(yS);       // creates Rcpp vector from SEXP
  Rcpp::IntegerMatrix nR(nS);       // creates Rcpp matrix from SEXP
  Rcpp::List bR(bS), paR(parS), prR(prS);
  
  bool aux_sw = true;
  unsigned nP = Rcpp::as<unsigned>(paR["npart"]);
  unsigned nM = Rcpp::as<unsigned>(paR["nstat"]);
  arma::mat pathMx,sdMx;
  if (paR.containsElementNamed("path"))
  {
    Rcpp::NumericMatrix pR = paR["path"];
    pathMx = arma::mat(pR.begin(), pR.nrow(), pR.ncol(), false);
    Rcpp::NumericMatrix sdR = paR["sd"];
    sdMx = arma::mat(sdR.begin(), sdR.nrow(), sdR.ncol(), false);
  }
  else if (paR.containsElementNamed("aux_alg"))
  {
    std::string aux_alg = Rcpp::as<std::string>(paR["aux_alg"]);
    aux_sw = (aux_alg.compare(0, aux_alg.length(), "Swendsen-Wang", 0, aux_alg.length()) == 0);
  }
  double alpha = 0.95;
  if (paR.containsElementNamed("alpha"))
  {
    alpha = Rcpp::as<double>(paR["alpha"]);
  }
  // minimum ABC tolerance
  double minE = DBL_EPSILON;
  if (paR.containsElementNamed("epsilon"))
  {
    minE = Rcpp::as<double>(paR["epsilon"]);
  }
  // minimum weighted empirical variance
  double minVar = DBL_EPSILON;
  if (paR.containsElementNamed("minVar"))
  {
    minVar = Rcpp::as<double>(paR["minVar"]);
  }
  // maximum number of SMC iterations
  unsigned max_iter = 100;
  if (paR.containsElementNamed("iter"))
  {
    max_iter = Rcpp::as<unsigned>(paR["iter"]);
  }
  // number of MCMC iterations to update mu, sigma and z
  unsigned aux_iter = nP;
  if (paR.containsElementNamed("aux"))
  {
    aux_iter = Rcpp::as<unsigned>(paR["aux"]);
  }

  if (prR.length() == 0)
  {
    throw std::invalid_argument("prior is empty");
  }
  int nvert = nR.nrow();
  int k = Rcpp::as<int>(prR["k"]);
  Rcpp::NumericVector prior_mu = prR["mu"];
  Rcpp::NumericVector prior_mu_sd = prR["mu.sd"];
  Rcpp::NumericVector prior_sd = prR["sigma"];
  Rcpp::NumericVector prior_sd_nu = prR["sigma.nu"];
  Rcpp::NumericVector prior_beta = prR["beta"];
  
  Rcpp::NumericVector yunique = Rcpp::unique(yR);
  Rcpp::IntegerVector ymatchR = Rcpp::match(yR, yunique);
  // no easy conversion from IntegerVector to uvec
  arma::uvec ymatch = unsign(ymatchR) - 1;
  arma::umat neigh = unsignMx(nR) - 1;
  
  // block index vectors are not symmetric
  std::vector<arma::uvec> blocks;
  blocks.reserve(bR.length());
  for (int b=0; b<bR.length(); b++)
  {
    Rcpp::IntegerVector block = bR[b];
    arma::uvec ublock = unsign(block - 1);
    blocks.push_back(ublock);
  }
  
  arma::colvec y(yR.begin(), yR.size(), false); // reuses memory and avoids extra copy
  arma::rowvec pr_mu(prior_mu.begin(), prior_mu.size(), false);
  arma::rowvec pr_mu_sd(prior_mu_sd.begin(), prior_mu_sd.size(), false);
  arma::rowvec pr_mu_tau = arma::pow(pr_mu_sd, -2);
  arma::rowvec pr_sd(prior_sd.begin(), prior_sd.size(), false);
  arma::rowvec pr_sd_nu(prior_sd_nu.begin(), prior_sd_nu.size(), false);
  arma::rowvec pr_sd_SS = pr_sd_nu % arma::square(pr_sd); // Schur product
  double pr_beta[2];
  pr_beta[0] = prior_beta[0];
  pr_beta[1] = prior_beta[prior_beta.size()-1];
  
  Rcpp::RNGScope scope;               // initialize random number generator
  arma::rowvec mu = rnorm(pr_mu,pr_mu_sd);
  arma::rowvec sd = arma::pow(rgamma(pr_sd_nu/2, pr_sd_SS/2), -0.5);
  arma::umat z = randomIndices(nvert, k);
  arma::umat alloc = arma::zeros<arma::umat>(nR.nrow(), k);
  arma::rowvec nZ(k), sumY(k), sqDiff(k);
  // warm up z (at beta=0) to avoid local minima
  gibbsLabelsNoData(neigh, blocks, z, alloc, 0.0);
  updateStats(y, z, nZ, sumY, sqDiff);
  for (unsigned p=0; p < aux_iter; p++)
  {
    mu = gibbsMeans(nZ, sumY, pr_mu, pr_mu_tau, sd);
    sd = gibbsStdDev(nZ, sumY, sqDiff, pr_sd_nu, pr_sd_SS, mu);
    arma::mat alpha = dnorm(yunique, ymatch, mu, sd);
    gibbsLabels(neigh, blocks, z, alloc, 0.0, alpha);
    updateStats(y, z, nZ, sumY, sqDiff);
  }
  unsigned s_z = sum_ident(z, neigh, blocks);
  Rcpp::NumericVector particles = Rcpp::runif(nP, pr_beta[0], pr_beta[1]);
  arma::vec logWeights(nP);
  logWeights.fill(-log((double)nP));
  Rcpp::NumericMatrix pseudo = updatePseudo(particles, nM, k, neigh, blocks, pathMx, sdMx, aux_sw);
  
  arma::mat mu_save = arma::zeros(aux_iter, k); // history of simulated values of mu
  arma::mat sd_save = arma::zeros(aux_iter, k); // history of simulated values of sigma
  
  unsigned it = 0;
  std::vector<double> accept, ess, epsilon, sum, var;
  accept.reserve(20);
  accept.push_back(nP);
  ess.reserve(20);
  ess.push_back(nP);
  epsilon.reserve(20);
  sum.reserve(20);
  sum.push_back(s_z);
  var.reserve(20);
  double wMean = weightedMean(particles, logWeights);
  double wVar = weightedVariance(particles, logWeights, wMean);
  var.push_back(wVar);
  time(&timer2);
  Rcpp::Rcout << "Initialization took " << difftime(timer2,timer) << "sec\n";
  Rcpp::IntegerVector idx;
  
  do
  {
    it++;
    Rcpp::Rcout << "Iteration " << it << "\n";
    Rcpp::checkUserInterrupt();
    // adaptively select epsilon and update importance weights    
    time(&timer);
    if (epsilon.size() == 0)
    {
      double e = 0.0;
      for (int i=0; i<pseudo.nrow(); i++)
      {
        for (int m=0; m<pseudo.ncol(); m++)
        {
          double d = fabs(pseudo(i,m) - s_z);
          if (d > e) e = d;
        }
      }
      epsilon.push_back(e);
    }
    logWeights = updateImportanceWeights(logWeights,epsilon,ess,pseudo,s_z,alpha);
    time(&timer2);
    Rcpp::Rcout << "Took " << difftime(timer2,timer) << "sec to update epsilon=" << epsilon.back() << " (ESS=" << ess.back() << ")\n";

    // resample if necessary
    if (ess.back() < 0.5*(double)nP)
    {
      time(&timer);
      idx = resample_resid(particles, logWeights, pseudo);
      Rcpp::NumericVector newP(nP);
      Rcpp::NumericMatrix newM(nP, nM);
      for (unsigned p=0; p<nP; p++)
      {
        int j=idx[p];
        newP[p] = particles[j];
        for (unsigned m=0; m<nM; m++)
        {
          newM(p,m) = pseudo(j,m);
        }
      }
      particles = newP;
      pseudo = newM;
      logWeights.fill(-log((double)nP));
      ess.push_back(nP);
      time(&timer2);
      Rcpp::Rcout << "Took " << difftime(timer2,timer) << "sec to resample " << nP << " particles\n";
    }
    // update particles (& pseudo-data) using random walk proposals
    wMean = weightedMean(particles, logWeights);
    wVar = weightedVariance(particles, logWeights, wMean);
    var.push_back(wVar);
    unsigned acc = rwmhParticles(particles, logWeights, pseudo, sqrt(2*wVar), pr_beta, s_z,
                          epsilon.back(), k, neigh, blocks, pathMx, sdMx, aux_sw);
    accept.push_back(acc);
    time(&timer);
    Rcpp::Rcout << "Took " << difftime(timer,timer2) << "sec for " << acc << " RWMH updates (bw=" << sqrt(2*wVar) << ")\n";

    // update sufficient statistic using a random sample of particles
    Rcpp::NumericVector aux_beta = subsample(particles, logWeights, aux_iter);
    for (unsigned p=0; p < aux_iter; p++)
    {
      // update labels
      arma::mat alpha = dnorm(yunique, ymatch, mu, sd);
      gibbsLabels(neigh, blocks, z, alloc, aux_beta[p], alpha);
      updateStats(y, z, nZ, sumY, sqDiff);

      // update means
      mu = gibbsMeans(nZ, sumY, pr_mu, pr_mu_tau, sd);
      mu_save.row(p) = mu;

      // update standard deviations
      sd = gibbsStdDev(nZ, sumY, sqDiff, pr_sd_nu, pr_sd_SS, mu);
      sd_save.row(p) = sd;
    }
    s_z = sum_ident(z, neigh, blocks);
    sum.push_back(s_z);
    time(&timer2);
    Rcpp::Rcout << "Took " << difftime(timer2,timer) << "sec for " << aux_iter << " iterations to calculate S(z)=" << s_z << "\n";

  } while (accept[it]/((double)nP) > 0.015 && epsilon[it] > minE && it < max_iter && wVar > minVar);
  
  return Rcpp::List::create(
    Rcpp::Named("alloc")   = alloc,     // count of allocations to each component
    Rcpp::Named("mu")      = mu_save,   // means of each mixture component
    Rcpp::Named("sigma")   = sd_save,   // standard deviations
    Rcpp::Named("beta")    = particles, // inverse temperature
    Rcpp::Named("wt")      = arma::exp(logWeights), // importance weights
    Rcpp::Named("accept")  = accept,    // M-H acceptance
    Rcpp::Named("epsilon") = epsilon,   // ABC tolerance
    Rcpp::Named("ess")     = ess,       // effective sample size
    Rcpp::Named("sum")     = sum,       // sufficient statistic
    Rcpp::Named("variance")= var,       // empirical variance of the particles
    Rcpp::Named("pseudo")  = pseudo,   // pseudo-data
    Rcpp::Named("idx")     = idx
  );
  END_RCPP
}

// compare two particles based on their normalized distance from the sufficient statistic of the data
class ParticleComparitor
{
  double mSz;
  arma::mat mPath;
  arma::mat mSD;
public:

  // ctor
  ParticleComparitor(double sum_z, arma::mat path, arma::mat sdMx)
  {
    mSz = sum_z;
    mPath = path;
    mSD = sdMx;
  }
  
  void setSz(double sum_z)
  {
    mSz = sum_z;
  }
  
  double calcDist(const double beta) const
  {
    unsigned ix = 0;
    while(mPath(0,ix) <= beta) ix++;
    double mu_Sz = interp(beta, ix-1, mPath);
    double sd_Sz = interp(beta, ix-1, mSD);
    return(fabs((mSz - mu_Sz)/sd_Sz));
  }
  
  bool operator() (double i, double j)
  {
    return (calcDist(i) < calcDist(j));
  }
};

arma::vec calcWeights(const arma::vec &particles, const ParticleComparitor &pcomp)
{
  arma::vec weights(particles.size());
  for (unsigned i=0; i<particles.size(); i++)
  {
    weights[i] = pcomp.calcDist(particles[i]);
  }
  return(weights);
}

// this recursive function shifts the position of the label in row i.
// when it reaches the end of a row, it propagates changes to the row below.
void increment_labels(arma::umat &z, unsigned i)
{
  // find the current position
  unsigned j=0;
  while (j < z.n_cols && z(i,j) != 1) j++;

  z(i,j) = 0;
  if (j == z.n_cols-1)
  {
    z(i,0) = 1;
    increment_labels(z, i+1);
  }
  else
  {
    z(i,j+1) = 1;
  }
}

SEXP initSedki(SEXP yS, SEXP nS, SEXP bS, SEXP parS, SEXP prS)
{
BEGIN_RCPP
  time_t timer,timer2;
  time(&timer);
  Rcpp::NumericVector yR(yS);       // creates Rcpp vector from SEXP
  Rcpp::IntegerMatrix nR(nS);       // creates Rcpp matrix from SEXP
  Rcpp::List bR(bS), paR(parS), prR(prS);
  
  unsigned nP = Rcpp::as<unsigned>(paR["npart"]);
  arma::mat pathMx,sdMx;
  if (paR.containsElementNamed("path"))
  {
    Rcpp::NumericMatrix pR = paR["path"];
    pathMx = arma::mat(pR.begin(), pR.nrow(), pR.ncol(), false);
    Rcpp::NumericMatrix sdR = paR["sd"];
    sdMx = arma::mat(sdR.begin(), sdR.nrow(), sdR.ncol(), false);
  }
  
  if (prR.length() == 0)
  {
    throw std::invalid_argument("prior is empty");
  }
  int nvert = nR.nrow();
  int k = Rcpp::as<int>(prR["k"]);
  Rcpp::NumericVector prior_mu = prR["mu"];
  Rcpp::NumericVector prior_mu_sd = prR["mu.sd"];
  Rcpp::NumericVector prior_sd = prR["sigma"];
  Rcpp::NumericVector prior_sd_nu = prR["sigma.nu"];
  Rcpp::NumericVector prior_beta = prR["beta"];
  
  Rcpp::NumericVector yunique = Rcpp::unique(yR);
  Rcpp::IntegerVector ymatchR = Rcpp::match(yR, yunique);
  // no easy conversion from IntegerVector to uvec
  arma::uvec ymatch = unsign(ymatchR) - 1;
  arma::umat neigh = unsignMx(nR) - 1;
  
  // block index vectors are not symmetric
  std::vector<arma::uvec> blocks;
  blocks.reserve(bR.length());
  for (int b=0; b<bR.length(); b++)
  {
    Rcpp::IntegerVector block = bR[b];
    arma::uvec ublock = unsign(block - 1);
    blocks.push_back(ublock);
  }
  
  arma::colvec y(yR.begin(), yR.size(), false); // reuses memory and avoids extra copy
  arma::rowvec pr_mu(prior_mu.begin(), prior_mu.size(), false);
  arma::rowvec pr_mu_sd(prior_mu_sd.begin(), prior_mu_sd.size(), false);
  arma::rowvec pr_mu_tau = arma::pow(pr_mu_sd, -2);
  arma::rowvec pr_sd(prior_sd.begin(), prior_sd.size(), false);
  arma::rowvec pr_sd_nu(prior_sd_nu.begin(), prior_sd_nu.size(), false);
  arma::rowvec pr_sd_SS = pr_sd_nu % arma::square(pr_sd); // Schur product
  double pr_beta[2];
  pr_beta[0] = prior_beta[0];
  pr_beta[1] = prior_beta[prior_beta.size()-1];
  
  Rcpp::RNGScope scope;               // initialize random number generator
  arma::rowvec mu = rnorm(pr_mu,pr_mu_sd);
  arma::rowvec sd = arma::pow(rgamma(pr_sd_nu/2, pr_sd_SS/2), -0.5);
  arma::umat z = randomIndices(nvert, k);
  arma::mat mu_save = arma::zeros(nP, k); // history of simulated values of mu
  arma::mat sd_save = arma::zeros(nP, k); // history of simulated values of sigma
  arma::umat alloc = arma::zeros<arma::umat>(nR.nrow(), k);
  arma::rowvec nZ(k), sumY(k), sqDiff(k);

  Rcpp::NumericVector newBeta = Rcpp::runif(nP, pr_beta[0], pr_beta[1]);
  arma::vec particles(newBeta.begin(), newBeta.size());
  unsigned it = 0, s_z;
  std::vector<double> epsilon, sum, var;
  epsilon.reserve(20);
  sum.reserve(20);
  var.reserve(20);
  double varBeta = arma::var(particles);
  var.push_back(varBeta);
  time(&timer2);
  Rcpp::Rcout << "Initialization took " << difftime(timer2,timer) << "sec\n";
  Rcpp::IntegerVector idx;
  
  // update sufficient statistic
  for (unsigned p=0; p<nP; p++)
  {
    // update labels
    arma::mat alpha = dnorm(yunique, ymatch, mu, sd);
    gibbsLabels(neigh, blocks, z, alloc, particles[p], alpha);
    updateStats(y, z, nZ, sumY, sqDiff);
    
    // update means
    mu = gibbsMeans(nZ, sumY, pr_mu, pr_mu_tau, sd);
    mu_save.row(p) = mu;
    
    // update standard deviations
    sd = gibbsStdDev(nZ, sumY, sqDiff, pr_sd_nu, pr_sd_SS, mu);
    sd_save.row(p) = sd;
  }
  s_z = sum_ident(z, neigh, blocks);
  sum.push_back(s_z);
  ParticleComparitor pcomp(s_z, pathMx, sdMx);
  arma::vec weights = calcWeights(particles, pcomp);
  epsilon.push_back(weights.max());
  time(&timer);
  Rcpp::Rcout << "Took " << difftime(timer,timer2) << "sec to calculate S(z)=" << s_z << "\n";

  do
  {
    it++;
    particles.resize(it*nP);
    newBeta = Rcpp::runif(nP, pr_beta[0], pr_beta[1]);
    arma::vec newb(newBeta.begin(), newBeta.size(), false);
    particles.insert_rows((it-1)*nP, newb);
  time(&timer2);
  Rcpp::Rcout << "Took " << difftime(timer2,timer) << "sec to draw " << nP << " new values\n";
    pcomp.setSz(s_z);
    weights = calcWeights(particles, pcomp);
    arma::uvec wtOrder = arma::sort_index(weights);
    unsigned cutoff = wtOrder[nP-1];
    epsilon.push_back(weights[cutoff]);
  time(&timer);
  Rcpp::Rcout << "Took " << difftime(timer,timer2) << "sec to sort values according to weight\n";
    double sumBeta = 0.0;
    for (unsigned p=0; p<nP; p++)
    {
      unsigned idx = wtOrder[p];
      sumBeta += particles[idx];
      // update labels
      arma::mat alpha = dnorm(yunique, ymatch, mu, sd);
      gibbsLabels(neigh, blocks, z, alloc, particles[idx], alpha);
      updateStats(y, z, nZ, sumY, sqDiff);
      
      // update means
      mu = gibbsMeans(nZ, sumY, pr_mu, pr_mu_tau, sd);
      mu_save.row(p) = mu;
      
      // update standard deviations
      sd = gibbsStdDev(nZ, sumY, sqDiff, pr_sd_nu, pr_sd_SS, mu);
      sd_save.row(p) = sd;
    }
    s_z = sum_ident(z, neigh, blocks);
    sum.push_back(s_z);
  time(&timer2);
  Rcpp::Rcout << "Took " << difftime(timer2,timer) << "sec to calculate S(z)=" << s_z << "\n";
    double muBeta = sumBeta/((double)nP);
    double ssdBeta = 0.0;
    for (unsigned p=0; p<nP; p++)
    {
      unsigned idx = wtOrder[p];
      ssdBeta += pow(particles[idx] - muBeta, 2.0);
    }
    varBeta = ssdBeta/((double)nP-1);
    var.push_back(varBeta);
  time(&timer);
  Rcpp::Rcout << "Took " << difftime(timer,timer2) << "sec to calculate var(B)=" << varBeta << "\n";
  } while(epsilon.front() >= epsilon.back() && var.back() >= var.front()/2.0);

  return Rcpp::List::create(
    Rcpp::Named("alloc")   = alloc,     // count of allocations to each component
    Rcpp::Named("mu")      = mu_save,   // means of each mixture component
    Rcpp::Named("sigma")   = sd_save,   // standard deviations
    Rcpp::Named("beta")    = particles, // inverse temperature
    Rcpp::Named("wt")      = weights,   // importance weights
    Rcpp::Named("epsilon") = epsilon,   // ABC tolerance
    Rcpp::Named("sum")     = sum,       // sufficient statistic
    Rcpp::Named("variance")= var        // empirical variance of the particles
  );
  END_RCPP
}

SEXP testResample(SEXP vS, SEXP wS, SEXP pS)
{
BEGIN_RCPP
  time_t timer,timer2;
  time(&timer);
  Rcpp::NumericVector particles(vS), wR(wS); // creates Rcpp vector from SEXP
  Rcpp::NumericMatrix pseudo(pS);
  unsigned nP = particles.size();
  unsigned nM = pseudo.ncol();
  arma::vec logWeights(nP);
  for (unsigned i=0; i<nP; i++)
  {
    logWeights[i] = log(wR[i]);
  }

  Rcpp::IntegerVector idx = resample_resid(particles, logWeights, pseudo);
  Rcpp::NumericVector newP(nP);
  Rcpp::NumericMatrix newM(nP, nM);
  for (unsigned p=0; p<nP; p++)
  {
    int j=idx[p];
    newP[p] = particles[j];
    for (unsigned m=0; m<nM; m++)
    {
      newM(p,m) = pseudo(j,m);
    }
  }
  particles = newP;
  pseudo = newM;
  logWeights.fill(-log((double)nP));
  time(&timer2);
  Rcpp::Rcout << "Took " << difftime(timer2,timer) << "sec to resample " << nP << " particles\n";
    
  return Rcpp::List::create(
    Rcpp::Named("beta")    = particles, // inverse temperature
    Rcpp::Named("wt")      = arma::exp(logWeights), // importance weights
    Rcpp::Named("pseudo")  = pseudo,   // pseudo-data
    Rcpp::Named("idx")     = idx
  );
  END_RCPP
}

SEXP exactPotts(SEXP nS, SEXP bS, SEXP kS, SEXP betaS)
{
BEGIN_RCPP
  Rcpp::NumericVector beta(betaS);  // creates Rcpp vector from SEXP
  Rcpp::IntegerMatrix nR(nS);       // creates Rcpp matrix from SEXP
  Rcpp::List bR(bS);
  int k = Rcpp::as<int>(kS);

  // no easy conversion from IntegerMatrix to umat
  arma::umat neigh = unsignMx(nR) - 1;
  unsigned n = neigh.n_rows;

  // block index vectors are not symmetric
  std::vector<arma::uvec> blocks;
  blocks.reserve(bR.length());
  for (int b=0; b<bR.length(); b++)
  {
    Rcpp::IntegerVector block = bR[b];
    arma::uvec ublock = unsign(block - 1);
    blocks.push_back(ublock);
  }

  unsigned niter = (unsigned)pow((double)k, (double)n);
  arma::umat z = randomIndices(n, k);
  arma::vec exp_save = arma::zeros(beta.size());
  arma::vec esq_save = arma::zeros(beta.size());
  arma::vec var_save = arma::zeros(beta.size());
  arma::vec sum_save = arma::zeros(niter);
  arma::vec pl_exp = arma::zeros(beta.size());
  arma::vec pl_var = arma::zeros(beta.size());
  arma::vec pl_sum = arma::zeros(beta.size()); // pseudolikelihood is unnormalised...

  // initialization
  z.zeros();
  for (unsigned i=0; i < n; i++)
  {
    z(i,0) = 1;
  }
  sum_save(0) = sum_ident(z, neigh, blocks);
  arma::uvec e(z.n_rows-1); // allocation vector
  arma::mat ne = arma::zeros(z.n_cols, z.n_rows-1);  // counts of like neighbours
  neighbj(ne, e, z, neigh);
  for (int b=0; b < beta.size(); b++)
  {
    double log_pl = pseudolike(ne, e, beta[b], z.n_rows-1, z.n_cols);
    pl_sum(b) = exp(log_pl);
    pl_exp(b) = exp(log_pl + log(sum_save(0)));
    pl_var(b) = exp(log_pl + 2*log(sum_save(0)));
  }

  // enumerate all possible combinations of labels
  // and calculate the corresponding sufficient statistic
  for (unsigned it=1; it < niter; it++)
  {
    increment_labels(z, 0);
    sum_save(it) = sum_ident(z, neigh, blocks);
    neighbj(ne, e, z, neigh);
    for (int b=0; b < beta.size(); b++)
    {
      double log_pl = pseudolike(ne, e, beta[b], z.n_rows-1, z.n_cols);
      pl_sum(b) += exp(log_pl);
      pl_exp(b) += exp(log_pl + log(sum_save(it)));
      pl_var(b) += exp(log_pl + 2*log(sum_save(it)));
    }
    // check for interrupt every thousand iterations
    if (it % 1000 == 0) Rcpp::checkUserInterrupt();
  }

  // now calculate the expectation and the variance
  for (int b=0; b < beta.size(); b++)
  {
    arma::vec p = arma::exp(beta[b]*sum_save);
    p = p / arma::sum(p);
    // Schur (element-wise) product of 2 vectors
    exp_save[b] = arma::sum(p % sum_save);
    arma::vec psq = p % arma::pow(sum_save,2);
    esq_save[b] = arma::sum(psq);
    pl_exp(b) /= pl_sum(b);
    pl_var(b) /= pl_sum(b);
  }
  var_save = esq_save - arma::pow(exp_save, 2);
  pl_var = pl_var - arma::pow(pl_exp, 2);

  return Rcpp::List::create(
      Rcpp::Named("beta") = beta,            // inverse temperature
      Rcpp::Named("expectation") = exp_save, // expectation of the sufficient statistic
      Rcpp::Named("variance") = var_save,    // variance of the sufficient statistic
      Rcpp::Named("exp_PL") = pl_exp,        // pseudolikelihood approximation of the expectation of S(z)
      Rcpp::Named("var_PL") = pl_var         // PL approx. of the variance of the sufficient statistic
  );
END_RCPP
}
