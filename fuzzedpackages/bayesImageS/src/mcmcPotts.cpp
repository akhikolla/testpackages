// ----------------------------------------------------------------------
// This file is part of the R package bayesImageS. It contains
// implementations of Metropolis-Hastings algorithms for image
// segmentation using a hidden Potts model.
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
#include "mcmcPotts.h"
#include "PottsUtil.h"

// updates inverse temperature using approximate exchange algorithm
// Lionel Cucala, J-M Marin, C. P. Robert & D. M. Titterington (2009)
unsigned exchangeBeta(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                      const arma::umat & z, double & beta,
                      const double prior_beta[2], const unsigned aux, const bool useSW, const bool swapAux, const double bw)
{
  // random walk proposal for B' ~ N(B, 0.01^2)
  double bprime = rwmh(beta, bw, prior_beta);

  // approximate W' using Swendsen-Wang algorithm
  arma::umat alloc = arma::zeros<arma::umat>(z.n_rows-1, z.n_cols);
  arma::umat w;
  if (swapAux) w = z;
  else w = randomIndices(z.n_rows-1, z.n_cols);

  for (unsigned i=0; i<aux; i++)
  {
    if (useSW)
    {
      swLabelsNoData(neigh, blocks, bprime, w.n_cols, w, alloc);
    }
    else
    {
      gibbsLabelsNoData(neigh, blocks, w, alloc, bprime);
    }
  }

  // calculate Metropolis-Hastings ratio
  double sum_z = sum_ident(z, neigh, blocks);
  double sum_w = sum_ident(w, neigh, blocks);
  double log_ratio = (bprime-beta)*sum_z + (beta-bprime)*sum_w;

  Rcpp::Rcout << exp(log_ratio);
  // accept/reject
  if (unif_rand() < exp(log_ratio))
  {
    beta = bprime;
    Rcpp::Rcout << "\t1\t" << beta << "\n";
    return 1;
  }
  Rcpp::Rcout << "\t0\n";
  return 0;
}

// updates inverse temperature using approximate exchange algorithm
// but interpolating between precomputed values of the sufficient statistic
unsigned accelExchange(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                       const arma::mat & path, const arma::mat & sdMx,
                       const arma::umat & z, double & beta, const double prior_beta[2],
                       const unsigned accept)
{
  unsigned ix = 0;
  while(path(0,ix) <= beta) ix++;
  // initial bandwidth must cover the entire prior
  double bw = (prior_beta[1] - prior_beta[0])/3;
  if (accept > 0) {
    bw = 6 / interp(beta, ix-1, sdMx);
    Rcpp::Rcout << "(BW " << bw << ") ";
  }
  
  // random walk proposal for B' ~ N(B, 0.01^2)
  double bprime = rwmh(beta, bw, prior_beta);

  // approximate W' by interpolation
  ix = 0;
  while(path(0,ix) <= bprime) ix++;
  double sprime = interp(bprime, ix-1, path);
  double sd_prime = bw;
  if (accept > 0) {
    sd_prime = 6 / interp(bprime, ix-1, sdMx);
  }
  Rcpp::Rcout << sprime << " (" << bprime << ") ";

  // calculate Metropolis-Hastings ratio
  double sum_z = sum_ident(z, neigh, blocks);
  double log_ratio = (bprime-beta)*sum_z + (beta-bprime)*sprime;
  if (accept > 0) {
    log_ratio += ::Rf_dnorm4(beta, bprime, sd_prime, 1) - ::Rf_dnorm4(bprime, beta, bw, 1);
  }
  Rcpp::Rcout << exp(log_ratio);
  // accept/reject
  if (unif_rand() < exp(log_ratio))
  {
    beta = bprime;
    Rcpp::Rcout << "\t1\t" << beta << "\n";
    return 1;
  }
  Rcpp::Rcout << "\t0\n";
  return 0;
}

// updates inverse temperature using approximate Bayesian computation (ABC)
// A Grelaud, C P Robert, J-M Marin, F Rodolphe & J-F Taly (2009)
unsigned abcBeta(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                 const arma::umat & z, double & beta, const double prior_beta[2],
                 const unsigned aux, const bool useSW, const bool swapAux, const double bw, const double epsilon)
{
  // random walk proposal for B' ~ N(B, 0.01^2)
  double bprime = rwmh(beta, bw, prior_beta);
  
  // approximate W' using Swendsen-Wang algorithm
  arma::umat alloc = arma::zeros<arma::umat>(z.n_rows-1, z.n_cols);
  arma::umat w;
  if (swapAux) w = z;
  else w = randomIndices(z.n_rows-1, z.n_cols);

  for (unsigned i=0; i<aux; i++)
  {
    if (useSW)
    {
      swLabelsNoData(neigh, blocks, bprime, w.n_cols, w, alloc);
    }
    else
    {
      gibbsLabelsNoData(neigh, blocks, w, alloc, bprime);
    }
  }

  // accept/reject
  double sum_z = sum_ident(z, neigh, blocks);
  double sum_w = sum_ident(w, neigh, blocks);
  double delta = fabs(sum_w - sum_z);
  Rcpp::Rcout << delta;
  if (delta < epsilon)
  {
    beta = bprime;
    Rcpp::Rcout << "\t1\t" << beta << "\n";
    return 1;
  }
  Rcpp::Rcout << "\t0\n";
  return 0;
}

// updates inverse temperature using approximate Bayesian computation (ABC-MCMC)
// but interpolating between precomputed values of the sufficient statistic
unsigned accelABC_MCMC(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                  const arma::mat & path, const arma::mat & sdMx,
                  const arma::umat & z, double & beta, const double prior_beta[2],
                  double epsilon, unsigned accept)
{
  unsigned ix = 0;
  while(path(0,ix) <= beta) ix++;
  // initial bandwidth must cover the entire prior
  double bw = (prior_beta[1] - prior_beta[0])/3;
  if (accept > 10) {
    bw = 3 / interp(beta, ix-1, sdMx);
    epsilon = R::qnorm5(0.999, 0, 3/bw, 1, 0);
  }
  Rcpp::Rcout << "(BW " << bw << "; e " << epsilon << ") ";
  
  // random walk proposal for B' ~ N(B, 0.01^2)
  double bprime = rwmh(beta, bw, prior_beta);

  // approximate the sufficient statistic
  ix = 0;
  while(path(0,ix) <= bprime) ix++;
  double sprime = interp(bprime, ix-1, path);
  double sd_prime = bw;
  if (accept > 10) {
    sd_prime = 3 / interp(bprime, ix-1, sdMx);
  }
  Rcpp::Rcout << sprime << " (" << bprime << ") - ";

  // accept/reject
  double sum_z = sum_ident(z, neigh, blocks);
  double delta = fabs(sprime - sum_z);
  Rcpp::Rcout << sum_z << " (" << beta << ") = " << delta;
  double log_ratio = 0;
  if (accept > 10) {
    log_ratio = ::Rf_dnorm4(beta, bprime, sd_prime, 1) - ::Rf_dnorm4(bprime, beta, bw, 1);
  }
  if ((unif_rand() < exp(log_ratio)) && (delta < epsilon))
  {
    beta = bprime;
    Rcpp::Rcout << "\t*\n";
    return 1;
  }
  Rcpp::Rcout << "\t-\n";
  return 0;
}

// updates inverse temperature using the approximate Bayesian computation (ABC)
// rejection sampler of Pritchard et al. (1999) and Grelaud et al. (2009),
// but interpolating between precomputed values of the sufficient statistic
unsigned accelABC(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                  const arma::mat & path, const arma::mat & sdMx,
                  const arma::umat & z,  double & beta, const double prior_beta[2],
                  const double epsilon)
{
  // independent proposal for B' ~ pr(B)
  double bprime = ::Rf_runif(prior_beta[0], prior_beta[1]);
  
  // approximate E[S(z)|bprime] using linear interpolation
  unsigned ix = 0;
  while(path(0,ix) <= bprime) ix++;
  double mu_Sz = interp(bprime, ix-1, path);
  double sd_Sz = interp(bprime, ix-1, sdMx);

  // accept/reject using normalized difference
  double sum_z = sum_ident(z, neigh, blocks);
  double delta = fabs((sum_z - mu_Sz)/sd_Sz);
  Rcpp::Rcout << sum_z << " (" << bprime << ") = " << delta;
  if (delta < epsilon)
  {
    beta = bprime;
    Rcpp::Rcout << "\t*\n";
    return 1;
  }
  Rcpp::Rcout << "\t-\n";
  return 0;
}

double calcApproxVar(const double beta, const double bcrit, const double v0, const double vmax1,
                     const double vmax2, const double phi1, const double phi2)
{
  if (beta <= bcrit)
  {
    return v0 + (vmax1 - v0)*exp(-phi1*sqrt(bcrit - beta));
  }
  else
  {
    return vmax2*exp(-phi2*sqrt(beta - bcrit));
  }
}

double calcApproxExp(const double beta, const double bcrit, const double v0, const double vmax1,
                     const double vmax2, const double phi1, const double phi2, const double e0, const double ecrit)
{
  if (beta <= bcrit)
  {
    double sqrtBcritPhi = sqrt(bcrit)*phi1;
    double sqrtBdiffPhi = sqrt(bcrit - beta)*phi1;
    return e0 + beta*v0 - ((2*(vmax1-v0))/pow(phi1,2.0))*((sqrtBcritPhi + 1)/exp(sqrtBcritPhi) - (sqrtBdiffPhi + 1)/exp(sqrtBdiffPhi));
  }
  else
  {
    double sqrtBdiff = sqrt(beta - bcrit);
    return ecrit - ((2*vmax2)/phi2)*(sqrtBdiff/exp(phi2*sqrtBdiff) + (exp(-phi2*sqrtBdiff) - 1)/phi2);
  }
}

unsigned accelAuxModel(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
              const arma::umat & z, double & beta, const double prior_beta[2], const double bw,
              const double bcrit, const double ecrit, const double e0, const double v0,
              const double vmax1, const double vmax2, const double phi1, const double phi2, const double sdMult)
{
  // random walk proposal for B' ~ N(B, 0.01^2)
  double bprime = rwmh(beta, bw, prior_beta);

  // approximate E[S(z)|beta] & V[S(z)|beta] using binding functions
  double exp_b0 = calcApproxExp(beta, bcrit, v0, vmax1, vmax2, phi1, phi2, e0, ecrit);
  double exp_b1 = calcApproxExp(bprime, bcrit, v0, vmax1, vmax2, phi1, phi2, e0, ecrit);
  double sd_b0 = sdMult*sqrt(calcApproxVar(beta, bcrit, v0, vmax1, vmax2, phi1, phi2));
  double sd_b1 = sdMult*sqrt(calcApproxVar(bprime, bcrit, v0, vmax1, vmax2, phi1, phi2));

  // accept/reject using parametric auxiliary model
  double sum_z = sum_ident(z, neigh, blocks);
  double log_ratio = ::Rf_dnorm4(sum_z, exp_b1, sd_b1, 1) - ::Rf_dnorm4(sum_z, exp_b0, sd_b0, 1);
  if (log(unif_rand()) < log_ratio)
  {
    beta = bprime;
//    Rcpp::Rcout << "\t1\t" << beta << "\n";
    return 1;
  }
//  Rcpp::Rcout << "\t0\n";
  return 0;  
}

// area under a trapezium
double trapezium(double fb0, double fb1, double diff)
{
  return(diff*(fb0 + fb1)/2);
}

// numerical integration using the trapezoidal rule
double quadrature(const double beta, const double bprime, const arma::mat & path)
{
  double a = std::min(beta,bprime);
  double b = std::max(beta,bprime);
  double sum = 0;
  unsigned ix = 0;
  while(path(0,ix) <= a) ix++;
  double fb0 = interp(a, ix-1, path);

  while(path(0,ix) < b)
  {
    sum += trapezium(fb0, path(1,ix), path(0,ix) - a);
    a = path(0,ix);
    fb0 = path(1,ix);
    ix++;
  }
  sum += trapezium(fb0, interp(b, ix-1, path), b - a);
  return(beta < bprime ? sum : -sum);
}

// updates inverse temperature using path sampling (Gelman & Meng, 1998)
unsigned pathBeta(const arma::umat & neigh, const std::vector<arma::uvec> & blocks, const arma::mat & path,
                  const arma::umat & z, double & beta, const double prior_beta[2], const double bw)
{
  // random walk proposal for B' ~ N(B, 0.01^2)
  double bprime = rwmh(beta, bw, prior_beta);
  
  // approximate log(Z(bprime) / Z(beta))
  double log_ratio = quadrature(bprime, beta, path) + (bprime-beta) * sum_ident(z, neigh, blocks);

  // accept/reject
  if (unif_rand() < exp(log_ratio))
  {
    beta = bprime;
    return 1;
  }
  return 0;
}

// updates inverse temperature using point pseudo-likelihood
// Tobias Ryde'n & D. M. Titterington (1998)
unsigned pseudoBeta(const arma::umat & neigh, const std::vector<arma::uvec> & blocks,
                    const arma::umat & z, double & beta, const double prior_beta[2],
                    const double bw)
{
  // random walk proposal for B' ~ N(B, 0.01^2)
  double bprime = rwmh(beta, bw, prior_beta);
  
  arma::uvec e(z.n_rows-1); // allocation vector
  arma::mat ne = arma::zeros(z.n_cols, z.n_rows-1);  // counts of like neighbours
  neighbj(ne, e, z, neigh);

  double log_ratio = pseudolike(ne, e, bprime, z.n_rows-1, z.n_cols) - pseudolike(ne, e, beta, z.n_rows-1, z.n_cols);
  Rcpp::Rcout << " (" << bprime << ") = " << exp(log_ratio);
  // accept/reject
  if (log(unif_rand()) < log_ratio)
  {
    beta = bprime;
    Rcpp::Rcout << "\t*\n";
    return 1;
  }
  Rcpp::Rcout << "\t.\n";
  return 0;
}

// calculates the number of correctly-classified pixels
unsigned Dice(arma::umat labels, arma::umat truth)
{
  arma::umat intersection = labels % truth;
  return arma::accu(intersection);
}

SEXP mcmcPotts(SEXP yS, SEXP nS, SEXP bS, SEXP itS, SEXP biS, SEXP prS, SEXP mhS, SEXP truthS)
{
BEGIN_RCPP
  Rcpp::NumericVector yR(yS);       // creates Rcpp vector from SEXP
  Rcpp::IntegerMatrix nR(nS);       // creates Rcpp matrix from SEXP
  Rcpp::List bR(bS), prR(prS), mhR(mhS);
  unsigned niter = Rcpp::as<unsigned>(itS), nburn = Rcpp::as<unsigned>(biS);
  
  std::string alg = Rcpp::as<std::string>(mhR["algorithm"]);
  bool pseudo = (alg.compare(0, alg.length(), "pseudolikelihood", 0, alg.length()) == 0);
  bool abc = (alg.compare(0, alg.length(), "abc", 0, alg.length()) == 0);
  bool exchange = (alg.compare(0, alg.length(), "exchange", 0, alg.length()) == 0);
  bool path = (alg.compare(0, alg.length(), "path", 0, alg.length()) == 0);
  bool auxMod = (alg.compare(0, alg.length(), "auxiliary", 0, alg.length()) == 0);
  if (!pseudo && !exchange && !abc && !path && !auxMod)
  {
    throw std::invalid_argument("algorithm not supported");
  }
  bool sortMeans = false;
  if (mhR.containsElementNamed("sort"))
  {
    sortMeans = Rcpp::as<bool>(mhR["sort"]);
  }
  unsigned aux = 0;
  bool aux_sw = true, aux_swap=false;
  if (mhR.containsElementNamed("auxiliary"))
  {
    aux = Rcpp::as<unsigned>(mhR["auxiliary"]);
  }
  if (exchange || abc)
  {
    if (mhR.containsElementNamed("aux_alg"))
    {
      std::string aux_alg = Rcpp::as<std::string>(mhR["aux_alg"]);
      aux_sw = (aux_alg.compare(0, aux_alg.length(), "Swendsen-Wang", 0, aux_alg.length()) == 0);
    }
    if (mhR.containsElementNamed("aux_swap"))
    {
      aux_swap = Rcpp::as<bool>(mhR["aux_swap"]);
      Rcpp::Rcout << "Swapping auxiliary variable: " << aux_swap << "\n";
    }
  }
  double bcrit=0, ecrit=0, e0=0, v0=0, vmax1=0, vmax2=0, phi1=0, phi2=0, sdMult=1;
  if (auxMod)
  {
    bcrit = Rcpp::as<double>(mhR["bcrit"]);
    e0 = Rcpp::as<double>(mhR["E0"]);
    v0 = Rcpp::as<double>(mhR["V0"]);
    vmax1 = Rcpp::as<double>(mhR["Vmax1"]);
    vmax2 = Rcpp::as<double>(mhR["Vmax2"]);
    phi1 = Rcpp::as<double>(mhR["phi1"]);
    phi2 = Rcpp::as<double>(mhR["phi2"]);
    ecrit = Rcpp::as<double>(mhR["Ecrit"]);
    sdMult = Rcpp::as<double>(mhR["factor"]);
  }
  arma::mat pathMx, sdMx;
  if (mhR.containsElementNamed("path"))
  {
    Rcpp::NumericMatrix pR = mhR["path"];
    pathMx = arma::mat(pR.begin(), pR.nrow(), pR.ncol(), false);
    if (mhR.containsElementNamed("sd"))
    {
      Rcpp::NumericMatrix sdR = mhR["sd"];
      sdMx = arma::mat(sdR.begin(), sdR.nrow(), sdR.ncol(), false);
    }
  }
  double bw = Rcpp::as<double>(mhR["bandwidth"]);
  double epsilon = 0;
  if (mhR.containsElementNamed("epsilon"))
  {
    epsilon = Rcpp::as<double>(mhR["epsilon"]);
  }
  double init_beta = 0.0;
  if (mhR.containsElementNamed("init"))
  {
    init_beta = Rcpp::as<double>(mhR["init"]);
  }

  // whether to use the algorithm of Garthwaite, Fan & Sisson (2010)
  bool adaptGFS=true;
  double adaptTarget = 0.44; // see Roberts and Rosenthal (2001)
  if (mhR.containsElementNamed("adaptive"))
  {
    if (Rf_isNull(mhR["adaptive"]) || ISNA(mhR["adaptive"]))
    {
      adaptGFS = false;
    }
    else
    {
      std::string adapt_alg = Rcpp::as<std::string>(mhR["adaptive"]);
      adaptGFS = (adapt_alg.compare(0, adapt_alg.length(), "GFS", 0, adapt_alg.length()) == 0);
    }
    if (adaptGFS && mhR.containsElementNamed("target"))
    {
      adaptTarget = Rcpp::as<double>(mhR["target"]);
    }
  }

  if (prR.length() == 0)
  {
    throw std::invalid_argument("prior is empty");
  }
  int nvert = nR.nrow();
  if (nvert != yR.size())
  {
    throw std::invalid_argument("mismatch between observations and neighbourhood matrix");
  }
  int k = Rcpp::as<int>(prR["k"]);
  Rcpp::NumericVector prior_mu = prR["mu"];
  Rcpp::NumericVector prior_mu_sd = prR["mu.sd"];
  Rcpp::NumericVector prior_sd = prR["sigma"];
  Rcpp::NumericVector prior_sd_nu = prR["sigma.nu"];
  arma::mat lpr_xfield;
  if (prR.containsElementNamed("xfield"))
  {
    Rcpp::NumericMatrix prior_xfield = prR["xfield"];
    arma::mat pr_xfield(prior_xfield.begin(), nvert, k, false);
    lpr_xfield = arma::log(pr_xfield);
  }
  Rcpp::NumericVector prior_beta = prR["beta"];

  Rcpp::NumericVector yunique = Rcpp::unique(yR);
  Rcpp::IntegerVector ymatchR = Rcpp::match(yR, yunique);
  // no easy conversion from IntegerVector to uvec
  arma::uvec ymatch = unsign(ymatchR) - 1;
  arma::umat neigh = unsignMx(nR) - 1;
  arma::umat truth;
  if (!Rf_isNull(truthS)) {
    Rcpp::IntegerMatrix tR(truthS);
    truth = unsignMx(tR);
    truth.insert_rows(truth.n_rows,1); // append a row of zeros
  }

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
  double beta = init_beta; // use the supplied value only if it is within the prior support
  if (prior_beta[0] > init_beta) {
    beta = prior_beta[0];
  }
  
  arma::mat mu_save = arma::zeros(niter, k); // history of simulated values of mu
  arma::mat sd_save = arma::zeros(niter, k); // history of simulated values of sigma
  arma::vec beta_save = arma::zeros(niter);  // history of simulated values of beta
  arma::vec sum_save = arma::zeros(niter);   // sum of identical neighbours
  arma::vec dice_save = arma::zeros(niter);  // Dice similarity coefficient

  arma::umat alloc = arma::zeros<arma::umat>(nR.nrow(), k);
  arma::rowvec nZ(k), sumY(k), sqDiff(k);
  unsigned accept = 0, acc_old = 0;
  Rcpp::Rcout << niter << " iterations of " << alg << " discarding the first " << nburn << " as burn-in." << "\n";

  for (unsigned it=0; it<niter; it++){
    // reset labels after burn-in
    if (it == nburn) alloc.zeros();
    // check for interrupt every thousand iterations
    if (it % 1000 == 0) Rcpp::checkUserInterrupt();

    // update labels
    arma::mat alpha = dnorm(yunique, ymatch, mu, sd);
    if (!lpr_xfield.is_empty()) alpha = alpha + lpr_xfield;
    gibbsLabels(neigh, blocks, z, alloc, beta, alpha);
    updateStats(y, z, nZ, sumY, sqDiff);
    sum_save(it) = sum_ident(z, neigh, blocks);
    if (!Rf_isNull(truthS)) dice_save(it) = Dice(z, truth);

    // update means
    mu = gibbsMeans(nZ, sumY, pr_mu, pr_mu_tau, sd);
    if (sortMeans) mu = arma::sort(mu);
    mu_save.row(it) = mu;
    
    // update standard deviations
    sd = gibbsStdDev(nZ, sumY, sqDiff, pr_sd_nu, pr_sd_SS, mu);
    sd_save.row(it) = sd;
    
    // update inverse temperature
    if (path)
    {
      accept += pathBeta(neigh, blocks, pathMx, z, beta, pr_beta, bw);
    }
    else if (auxMod)
    {
      accept += accelAuxModel(neigh, blocks, z, beta, pr_beta, bw,
                              bcrit, ecrit, e0, v0, vmax1, vmax2, phi1, phi2, sdMult);
    }
    else
    {
      // slow algorithm, so print progress at each iteration
      Rcpp::Rcout << it << "of" << niter << "(bw " << bw << ")\t";
      if (pseudo)
      {
        accept += pseudoBeta(neigh, blocks, z, beta, pr_beta, bw);
      }
      else if (abc)
      {
        if (aux > 0)
        {
          accept += abcBeta(neigh, blocks, z, beta, pr_beta, aux, aux_sw, aux_swap, bw, epsilon);
        }
        else
        {
          //accept += accelABC(neigh, blocks, pathMx, sdMx, z, beta, pr_beta, epsilon);
          accept += accelABC_MCMC(neigh, blocks, pathMx, sdMx, z, beta, pr_beta, epsilon, accept);
        }
      }
      else
      {
        if (aux > 0)
        {
          accept += exchangeBeta(neigh, blocks, z, beta, pr_beta, aux, aux_sw, aux_swap, bw);
        }
        else
        {
          accept += accelExchange(neigh, blocks, pathMx, sdMx, z, beta, pr_beta, accept);
        }
      }

      // adaptive MCMC algorithm of Garthwaite, Fan & Sisson (2010)
      if (adaptGFS && accept > 10)
      {
        if (accept > acc_old)
        {
          bw = bw + bw/adaptTarget/it;
        }
        else
        {
          bw = bw - bw/(1-adaptTarget)/it;
        }
        acc_old = accept;
      }
    }
    beta_save[it] = beta;
  }
  
  Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("alloc") = alloc,    // count of allocations to each component
      Rcpp::Named("mu")    = mu_save,  // history of simulated values of mu
      Rcpp::Named("sigma") = sd_save,  // history of simulated values of sigma
      Rcpp::Named("beta")  = beta_save,// history of simulated values of beta
      Rcpp::Named("accept") = accept,  // M-H acceptance
      Rcpp::Named("sum") = sum_save    // sum of identical neighbours
  );
  if (!Rf_isNull(truthS)) {
    result["correct"] = dice_save;     // append to the Rcpp::List
  }
  if (adaptGFS) {
    result["bandwidth"] = bw;          // final M-H bandwidth produced by GFS algorithm
  }
  if (pseudo) {
    arma::uvec e(z.n_rows-1);
    arma::mat ne = arma::zeros(z.n_cols, z.n_rows-1);
    neighbj(ne, e, z, neigh);
    result["z"] = z;
    result["e"] = e;                  // allocation vector
    result["ne"] = ne;                // counts of like neighbours
  }
  return result;
END_RCPP
}

SEXP gibbsNorm(SEXP yS, SEXP itS, SEXP prS) {
BEGIN_RCPP
  unsigned niter = Rcpp::as<unsigned>(itS);
  unsigned k = 1;
  
  Rcpp::NumericVector yR(yS);       // creates Rcpp vector from SEXP
  Rcpp::NumericVector yunique = Rcpp::unique(yR);
  Rcpp::IntegerVector ymatchR = Rcpp::match(yR, yunique);
  // no easy conversion from IntegerVector to uvec
  arma::uvec ymatch = unsign(ymatchR - 1);
  arma::colvec y(yR.begin(), yR.size(), false); // reuses memory and avoids extra copy
  
  Rcpp::List prR(prS);
  double prior_mu = prR["mu"];
  double prior_mu_sd = prR["mu.sd"];
  double prior_sd = prR["sigma"];
  double prior_sd_nu = prR["sigma.nu"];
  arma::rowvec pr_mu(k);
  arma::rowvec pr_mu_tau(k);
  arma::rowvec pr_sd_nu(k);
  arma::rowvec pr_sd_SS(k);
  pr_mu[0] = prior_mu;
  pr_mu_tau[0] = pow(prior_mu_sd,-2);
  pr_sd_nu[0] = prior_sd_nu;
  pr_sd_SS[0] = prior_sd_nu * pow(prior_sd,2);
  
  Rcpp::RNGScope scope;               // initialize random number generator
  arma::rowvec mu = rnorm(pr_mu,arma::pow(pr_mu_tau,-0.5));
  arma::rowvec sd = arma::pow(rgamma(pr_sd_nu/2, pr_sd_SS/2), -0.5);
  
  arma::mat mu_save = arma::zeros(niter,k); // history of simulated values of mu
  arma::mat sd_save = arma::zeros(niter,k); // history of simulated values of sigma
  arma::rowvec nZ(k), sumY(k), sqDiff(k);
  arma::umat z = arma::ones<arma::umat>(y.n_elem,k);
  updateStats(y, z, nZ, sumY, sqDiff);

  for (unsigned it=0; it<niter; it++){
    // update means
    mu = gibbsMeans(nZ, sumY, pr_mu, pr_mu_tau, sd);
    mu_save.row(it) = mu;
    
    // update standard deviations
    sd = gibbsStdDev(nZ, sumY, sqDiff, pr_sd_nu, pr_sd_SS, mu);
    sd_save.row(it) = sd;
  }

  return Rcpp::List::create(
      Rcpp::Named("mu")    = mu_save,  // history of simulated values of mu
      Rcpp::Named("sigma") = sd_save   // history of simulated values of sigma
  );
END_RCPP
}

SEXP gibbsGMM(SEXP yS, SEXP itS, SEXP biS, SEXP prS) {
BEGIN_RCPP
  unsigned niter = Rcpp::as<unsigned>(itS), nburn = Rcpp::as<unsigned>(biS);
  
  Rcpp::NumericVector yR(yS);       // creates Rcpp vector from SEXP
  Rcpp::NumericVector yunique = Rcpp::unique(yR);
  Rcpp::IntegerVector ymatchR = Rcpp::match(yR, yunique);
  // no easy conversion from IntegerVector to uvec
  arma::uvec ymatch = unsign(ymatchR - 1);
  arma::colvec y(yR.begin(), yR.size(), false); // reuses memory and avoids extra copy
  
  Rcpp::List prR(prS);
  unsigned k = Rcpp::as<unsigned>(prR["k"]);
  Rcpp::NumericVector prior_lambda = prR["lambda"];
  Rcpp::NumericVector prior_mu = prR["mu"];
  Rcpp::NumericVector prior_mu_sd = prR["mu.sd"];
  Rcpp::NumericVector prior_sd = prR["sigma"];
  Rcpp::NumericVector prior_sd_nu = prR["sigma.nu"];
  arma::rowvec pr_lambda(prior_lambda.begin(), prior_lambda.size(), false);
  arma::rowvec pr_mu(prior_mu.begin(), prior_mu.size(), false);
  arma::rowvec pr_mu_sd(prior_mu_sd.begin(), prior_mu_sd.size(), false);
  arma::rowvec pr_mu_tau = arma::pow(pr_mu_sd, -2);
  arma::rowvec pr_sd(prior_sd.begin(), prior_sd.size(), false);
  arma::rowvec pr_sd_nu(prior_sd_nu.begin(), prior_sd_nu.size(), false);
  arma::rowvec pr_sd_SS = pr_sd_nu % arma::square(pr_sd); // Schur product
  
  Rcpp::RNGScope scope;               // initialize random number generator
  arma::rowvec mu = rnorm(pr_mu,pr_mu_sd);
  arma::rowvec sd = arma::pow(rgamma(pr_sd_nu/2, pr_sd_SS/2), -0.5);
  arma::rowvec wt(k);
  for (unsigned j=0; j<k; j++)
  {
    wt[j] = unif_rand();
  }
  wt = wt/arma::sum(wt); // normalize to sum to 1
  arma::umat z = randomIndices(y.n_elem, k);
  
  arma::mat mu_save = arma::zeros(niter,k); // history of simulated values of mu
  arma::mat sd_save = arma::zeros(niter,k); // history of simulated values of sigma
  arma::mat wt_save = arma::zeros(niter,k); // history of simulated values of lambda
  arma::umat alloc  = arma::zeros<arma::umat>(y.n_elem, k);
  arma::rowvec nZ(k), sumY(k), sqDiff(k);

  for (unsigned it=0; it<niter; it++){
    // reset labels after burn-in
    if (it == nburn) alloc.zeros();

    // update labels
    arma::mat alpha = dnorm(yunique, ymatch, mu, sd);
    classify(z, alloc, arma::log(wt), alpha);
    updateStats(y, z, nZ, sumY, sqDiff);

    // update means
    mu = gibbsMeans(nZ, sumY, pr_mu, pr_mu_tau, sd);
    mu_save.row(it) = mu;
    
    // update standard deviations
    sd = gibbsStdDev(nZ, sumY, sqDiff, pr_sd_nu, pr_sd_SS, mu);
    sd_save.row(it) = sd;
    
    // update mixture weights
    wt = gibbsDirichlet(nZ, pr_lambda);
    wt_save.row(it) = wt;
  }

  return Rcpp::List::create(
      Rcpp::Named("alloc") = alloc,    // count of allocations to each component
      Rcpp::Named("mu")    = mu_save,  // history of simulated values of mu
      Rcpp::Named("sigma") = sd_save,  // history of simulated values of sigma
      Rcpp::Named("lambda") = wt_save  // history of simulated values of lambda
  );
END_RCPP
}

SEXP mcmcPottsNoData(SEXP betaS, SEXP kS, SEXP nS, SEXP bS, SEXP itS, SEXP randS) {
BEGIN_RCPP
  Rcpp::IntegerMatrix nR(nS);       // creates Rcpp matrix from SEXP
  Rcpp::List bR(bS);
  unsigned niter = Rcpp::as<unsigned>(itS);
  int k = Rcpp::as<int>(kS);
  double beta = Rcpp::as<double>(betaS);
  bool randInit = Rcpp::as<bool>(randS);
  
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

  Rcpp::RNGScope scope;               // initialize random number generator
  arma::umat z = arma::zeros<arma::umat>(n+1,k);
  if (randInit)
  {
    z = randomIndices(n, k);
  }
  else
  {
    unsigned j = 0;
    for (unsigned i=0; i<n; i++)
    {
      z(i,j) = 1;
    }
  }
  arma::umat alloc  = arma::zeros<arma::umat>(neigh.n_rows, k);
  arma::vec sum_save = arma::zeros(niter);

  for (unsigned it=0; it<niter; it++){
    // update labels
    gibbsLabelsNoData(neigh, blocks, z, alloc, beta);
    sum_save(it) = sum_ident(z, neigh, blocks);
  }

  return Rcpp::List::create(
      Rcpp::Named("alloc") = alloc,    // count of allocations to each component
      Rcpp::Named("z")     = z,        // final sample from Gibbs distribution
      Rcpp::Named("sum") = sum_save    // sum of identical neighbours
  );
END_RCPP
}

SEXP swNoData(SEXP betaS, SEXP kS, SEXP nS, SEXP bS, SEXP itS, SEXP randS) {
BEGIN_RCPP
  Rcpp::IntegerMatrix nR(nS);       // creates Rcpp matrix from SEXP
  Rcpp::List bR(bS);
  unsigned niter = Rcpp::as<unsigned>(itS);
  int k = Rcpp::as<int>(kS);
  double beta = Rcpp::as<double>(betaS);
  bool randInit = Rcpp::as<bool>(randS);

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
  
  Rcpp::RNGScope scope;               // initialize random number generator
  arma::umat z = arma::zeros<arma::umat>(n+1,k);
  if (randInit)
  {
    z = randomIndices(n, k);
  }
  else
  {
    unsigned j = 0;
    for (unsigned i=0; i<n; i++)
    {
      z(i,j) = 1;
    }
  }
  arma::umat alloc  = arma::zeros<arma::umat>(neigh.n_rows, k);
  arma::vec sum_save = arma::zeros(niter);

  for (unsigned it=0; it<niter; it++){
    // update labels
    swLabelsNoData(neigh, blocks, beta, k, z, alloc);
    sum_save(it) = sum_ident(z, neigh, blocks);
  }

  return Rcpp::List::create(
      Rcpp::Named("alloc") = alloc,    // count of allocations to each component
      Rcpp::Named("z")     = z,        // final sample from Gibbs distribution
      Rcpp::Named("sum") = sum_save,    // sum of identical neighbours
      Rcpp::Named("neigh") = neigh
  );
END_RCPP
}

SEXP sufficientStat(SEXP zS, SEXP nS, SEXP bS, SEXP kS)
{
BEGIN_RCPP
  Rcpp::IntegerMatrix zR(zS), nR(nS);   // creates Rcpp matrix from SEXP
  Rcpp::List bR(bS);
  int k = Rcpp::as<int>(kS);

  // no easy conversion from IntegerMatrix to umat
  arma::umat labels = unsignMx(zR);
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

  unsigned niter = labels.n_cols;
  arma::umat z = randomIndices(neigh.n_rows, k);
  arma::vec sum_save = arma::zeros(niter);

  for (unsigned it=0; it < niter; it++)
  {
    z.zeros();
    for (unsigned x=0; x < neigh.n_rows; x++)
    {
      unsigned j = labels(x,it) - 1;
      z(x,j) = 1;
    }
    sum_save(it) = sum_ident(z, neigh, blocks);
  }

  return Rcpp::List::create(
      Rcpp::Named("sum") = sum_save    // sum of identical neighbours
  );
END_RCPP
}
