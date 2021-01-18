#include <boost/math/special_functions/gamma.hpp>
#include <random>
#include "RcppArmadillo.h"
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

//~ sample three integers among {0, 1, ..., n-1} --------------------------- ~//
const std::array<int, 3> choose3(const int n,
                                 std::default_random_engine& generator) {
  std::uniform_int_distribution<int> sampler1(0, n - 1);
  std::uniform_int_distribution<int> sampler2(0, n - 2);
  std::uniform_int_distribution<int> sampler3(0, n - 3);
  const int i1 = sampler1(generator);
  int i2 = sampler2(generator);
  int i3 = sampler3(generator);
  if(i3 == i2)
    i3 = n - 2;
  if(i3 == i1)
    i3 = n - 1;
  if(i2 == i1)
    i2 = n - 1;
  return {i1, i2, i3};
}

//~ sample two integers among {0, 1, ..., n-1} ----------------------------- ~//
const std::array<int, 2> choose2(const int n,
                                 std::default_random_engine& generator) {
  std::uniform_int_distribution<int> sampler1(0, n - 1);
  std::uniform_int_distribution<int> sampler2(0, n - 2);
  const int i1 = sampler1(generator);
  int i2 = sampler2(generator);
  if(i2 == i1) {
    i2 = n - 1;
  }
  return {i1, i2};
}

//~ Beta-quantiles for a vector `beta` ------------------------------------- ~//
Rcpp::NumericVector BetaQuantile(const double g,
                                 const double s,
                                 const double a,
                                 const double prob,
                                 const Rcpp::NumericVector beta) {
  Rcpp::NumericVector alpha = (1.0 - beta) / prob;
  Rcpp::NumericVector Q;
  if(g == 0.0) {
    Q = a - s * log(alpha);
  } else {
    Q = a + s / g * (pow(alpha, -g) - 1);
  }
  return Q;
}

arma::colvec BetaQuantileArma(const double g,
                              const double s,
                              const double a,
                              const double prob,
                              const arma::colvec& beta) {
  arma::colvec Q;
  if(g == 0.0) {
    Q = a - s * log((1.0 - beta) / prob);
  } else {
    Q = a + s / g * (pow((1.0 - beta) / prob, -g) - 1);
  }
  return Q;
}

//~ Calculates the Jacobian ------------------------------------------------ ~//
double Jacobian(const double g,
                const double s,
                const double a,
                const size_t Jnumb,
                const Rcpp::NumericVector X,
                std::default_random_engine& generator) {
  Rcpp::NumericMatrix Xchoose3(Jnumb, 3);
  const int n = X.size();
  for(size_t i = 0; i < Jnumb; i++) {
    const std::array<int, 3> indices = choose3(n, generator);
    Xchoose3(i, Rcpp::_) = Rcpp::NumericVector::create(
        X(indices[0]), X(indices[1]), X(indices[2]));
  }
  Rcpp::NumericMatrix Xdiff(Jnumb, 3);
  Xdiff(Rcpp::_, 0) = Xchoose3(Rcpp::_, 1) - Xchoose3(Rcpp::_, 2);
  Xdiff(Rcpp::_, 1) = Xchoose3(Rcpp::_, 2) - Xchoose3(Rcpp::_, 0);
  Xdiff(Rcpp::_, 2) = Xchoose3(Rcpp::_, 0) - Xchoose3(Rcpp::_, 1);
  double Jmean;
  if(g == 0.0) {
    const Rcpp::NumericVector Jvec =
        Xdiff(Rcpp::_, 0) * Xdiff(Rcpp::_, 1) * Xdiff(Rcpp::_, 2);
    Jmean = Rcpp::mean(Rcpp::abs(Jvec));
  } else {
    const Rcpp::NumericMatrix A = g * (Xchoose3 - a) / s;
    Rcpp::NumericVector Jmat0 = (log1p(A) * (1.0 + A) / g / g) * Xdiff;
    const Rcpp::NumericVector Jvec =
        Jmat0[Rcpp::Range(0, Jnumb - 1)] +
        Jmat0[Rcpp::Range(Jnumb, 2 * Jnumb - 1)] +
        Jmat0[Rcpp::Range(2 * Jnumb, 3 * Jnumb - 1)];
    Jmean = Rcpp::mean(Rcpp::abs(Jvec));
  }
  return Jmean;
}

double JacobianArma1(const double g,
                     const double s,
                     const arma::vec& X,
                     const int n,
                     const arma::mat& UpperTriOnes) {
  double Jmean;
  if(g == 0.0) {
    const arma::mat XiXj = (X % X) * X.t();
    Jmean =
        arma::accu(abs(XiXj - XiXj.t()) % UpperTriOnes) / (s * s * n * (n - 1));
  } else {
    const arma::vec A = g / s * X;
    const arma::mat XiXj = X * ((1 + A) % log1p(A)).t();
    Jmean = 2.0 * arma::accu(abs(XiXj - XiXj.t()) % UpperTriOnes) /
            (g * g * n * (n - 1));
  }
  return Jmean;
}

double JacobianArma2(const double g,
                     const double s,
                     const size_t Jnumb,
                     const arma::vec& X,
                     const int n,
                     std::default_random_engine& generator) {
  
  arma::mat Xchoose2(2, Jnumb);
  for(size_t i = 0; i < Jnumb; i++) {
    const std::array<int, 2> indices = choose2(n, generator);
    const arma::colvec2 col_i = {X.at(indices[0]), X.at(indices[1])};
    Xchoose2.col(i) = col_i;
  }
  arma::mat Xchoose2t = Xchoose2.t();
  
  double Jmean;
  if(g == 0.0) {
    const arma::vec Jvec = Xchoose2t.col(0) % Xchoose2t.col(1) %
                           (Xchoose2t.col(0) - Xchoose2t.col(1)) / (2.0 * s * s);
    Jmean = arma::mean(abs(Jvec));
  } else {
    const arma::mat A = g / s * Xchoose2t;
    const arma::vec Jvec =
        (Xchoose2t.col(0) % (1 + A.col(1)) % log1p(A.col(1)) -
         Xchoose2t.col(1) % (1 + A.col(0)) % log1p(A.col(0))) /
        g / g;
    Jmean = arma::mean(abs(Jvec));
  }
  return Jmean;
}

//~ log-density of fiducial distribution of (gamma,sigma,a) ---------------- ~//
const double log_gpd_dens(const double g,
                          const double s,
                          const double a,
                          Rcpp::NumericVector X,
                          const size_t Jnumb,
                          const size_t n,
                          std::default_random_engine& generator) {
  double log_density;
  X = X[X > a];
  const double Max = Rcpp::max(X - a);
  if(s > 0 && g > (-s / Max) && a > 0.0 && g > -0.5) {
    const double J = Jacobian(g, s, a, Jnumb, X, generator);
    if(g == 0.0) {
      log_density = -1 / s * Rcpp::sum(X - a) + log(J) - n * log(s + a);
    } else {
      log_density = Rcpp::sum((-1 / g - 1) * log1p(g * (X - a) / s)) + log(J) -
                    n * log(s + a);
    }
  } else {
    log_density = -INFINITY;
  }
  return log_density;
}

//~ log-density of fiducial distribution of (gamma,sigma) ------------------ ~//
const double log_gpd_densArma(const double g,
                              const double s,
                              const arma::vec& X,
                              const int n,
                              const size_t Jnumb,
                              std::default_random_engine& generator,
                              const arma::mat& UpperTriOnes) {
  double log_density;
  const double Max = arma::max(X);

  if(s > 0 && g > (-s / Max)) {
    const double J = n < 250 ? JacobianArma1(g, s, X, n, UpperTriOnes)
                             : JacobianArma2(g, s, Jnumb, X, n, generator);

    if(g == 0.0) {
      log_density = -arma::accu(X) / s + log(J) - n * log(s);
    } else {
      log_density =
          (-1 / g - 1) * arma::accu(log1p(g * X / s)) + log(J) - n * log(s);
    }
  } else {
    log_density = -INFINITY;
  }

  return log_density;
}

//~ distributions to be sampled -------------------------------------------- ~//
std::uniform_real_distribution<double> uniform(0.0, 1.0);
std::cauchy_distribution<double> cauchy(0.0, 1.0);

//~ propose a new (gamma,sigma) value or a new index for the threshold ----- ~//
std::array<double, 3> MCMCnewpoint(const double g,
                                   const double s,
                                   const double i_dbl,
                                   const double p1,
                                   const double p2,
                                   double lambda,
                                   const double sd_g,
                                   const double sd_s,
                                   const Rcpp::NumericVector X,
                                   const size_t Jnumb,
                                   const int n,
                                   std::default_random_engine& generator,
                                   std::poisson_distribution<int>& poisson1,
                                   std::poisson_distribution<int>& poisson2) {
  double MHratio, g_star, s_star;

  const int i = (int)(i_dbl + 0.5);
  double a = X(i - 1);
  int i_star;

  if(uniform(generator) > p1) {
    int plus_minus = n;
    double dens_pois_star, dens_pois;

    if(uniform(generator) < p2) {
      while(plus_minus > n - i - 10) {
        plus_minus = poisson1(generator);
      }
      i_star = i + plus_minus;
      dens_pois_star = p2 / boost::math::gamma_q(n - i - 9, lambda);
      dens_pois = (1.0 - p2) / boost::math::gamma_q(i_star, lambda);
    } else {
      lambda = lambda < i_dbl ? lambda : i_dbl;
      while(plus_minus > i - 1) {
        plus_minus = poisson2(generator);
      }
      i_star = i - plus_minus;
      dens_pois_star = (1.0 - p2) / boost::math::gamma_q(i, lambda);
      dens_pois = p2 / boost::math::gamma_q(n - i_star - 9, lambda);
    }

    const double a_star = X(i_star - 1);
    g_star = g;
    s_star = s + g * (a_star - a);
    MHratio = exp(log_gpd_dens(g_star, s_star, a_star, X, Jnumb, n, generator) -
                  log_gpd_dens(g, s, a, X, Jnumb, n, generator)) *
              dens_pois / dens_pois_star;
  } else {
    g_star = g + sd_g * cauchy(generator);
    s_star = s + sd_s * cauchy(generator);
    i_star = i;

    MHratio = exp(log_gpd_dens(g_star, s_star, a, X, Jnumb, n, generator) -
                  log_gpd_dens(g, s, a, X, Jnumb, n, generator));
  }

  std::array<double, 3> newpoint;
  if(!std::isnan(MHratio) && !std::isinf(MHratio) &&
     uniform(generator) < MHratio) {
    newpoint = {g_star, s_star, (double)i_star};
  } else {
    newpoint = {g, s, (double)i};
  }

  return newpoint;
}

arma::colvec2 MCMCnewpointArma(const double g,
                               const double s,
                               const double sd_g,
                               const double sd_s,
                               const arma::vec& X,
                               const int n,
                               const size_t Jnumb,
                               std::default_random_engine& generator,
                               const arma::mat& UpperTriOnes) {
  const double g_star = g + sd_g * cauchy(generator);
  const double s_star = s + sd_s * cauchy(generator);
  const double MHratio = exp(
      log_gpd_densArma(g_star, s_star, X, n, Jnumb, generator, UpperTriOnes) -
      log_gpd_densArma(g, s, X, n, Jnumb, generator, UpperTriOnes));
  arma::colvec2 newpoint;
  if(!std::isnan(MHratio) && !std::isinf(MHratio) &&
     uniform(generator) < MHratio) {
    newpoint = {g_star, s_star};
  } else {
    newpoint = {g, s};
  }
  return newpoint;
}

//~ helper function for MCMCchain ------------------------------------------ ~//
Rcpp::NumericVector concat(const double g,
                           const double s,
                           const double i,
                           const Rcpp::NumericVector beta,
                           const size_t lbeta) {
  Rcpp::NumericVector out(3 + lbeta);
  out(0) = g;
  out(1) = s;
  out(2) = i;
  for(size_t k = 3; k < 3 + lbeta; k++) {
    out(k) = beta(k - 3);
  }
  return out;
}

//~ function that runs the MCMC chain -------------------------------------- ~//
// [[Rcpp::export]]
Rcpp::NumericMatrix MCMCchain(Rcpp::NumericVector X,
                              const Rcpp::NumericVector beta,
                              const double g,
                              const double s,
                              const double a,
                              const int i,
                              const double p1,
                              const double p2,
                              const double lambda1,
                              const double lambda2,
                              const double sd_g,
                              const double sd_s,
                              const size_t niter,
                              const size_t nburnin,  // almost not used here
                              const size_t Jnumb,
                              const unsigned seed) {
  std::default_random_engine generator(seed);

  X = X - X(0);
  const size_t lbeta = beta.size();
  const double i_dbl = (double)i;
  const int n = X.size();

  Rcpp::NumericMatrix xt(niter, 3 + lbeta);
  xt(0, Rcpp::_) =
      concat(g, s, i_dbl, BetaQuantile(g, s, a, 1.0 - i_dbl / n, beta), lbeta);

  std::poisson_distribution<int> poisson1(lambda1);
  std::poisson_distribution<int> poisson2(lambda2);
  std::poisson_distribution<int> poisson3(i_dbl);

  double lambda;

  for(size_t j = 0; j < niter - 1; j++) {
    bool b = j % 10 == 0;
    lambda = b ? lambda2 : lambda1;
    std::array<double, 3> gsi;
    if(lambda < i_dbl) {
      if(b) {
        gsi = MCMCnewpoint(xt(j, 0), xt(j, 1), xt(j, 2), p1, p2, lambda, sd_g,
                           sd_s, X, Jnumb, n, generator, poisson2, poisson2);
      } else {
        gsi = MCMCnewpoint(xt(j, 0), xt(j, 1), xt(j, 2), p1, p2, lambda, sd_g,
                           sd_s, X, Jnumb, n, generator, poisson1, poisson1);
      }
    } else {
      if(b) {
        gsi = MCMCnewpoint(xt(j, 0), xt(j, 1), xt(j, 2), p1, p2, lambda, sd_g,
                           sd_s, X, Jnumb, n, generator, poisson2, poisson3);
      } else {
        gsi = MCMCnewpoint(xt(j, 0), xt(j, 1), xt(j, 2), p1, p2, lambda, sd_g,
                           sd_s, X, Jnumb, n, generator, poisson1, poisson3);
      }
    }
    xt(j + 1, Rcpp::_) =
        concat(gsi[0], gsi[1], gsi[2],
               BetaQuantile(gsi[0], gsi[1], X((int)(gsi[2] + 0.5) - 1),
                            1.0 - gsi[2] / n, beta),
               lbeta);
  }

  xt = xt(Rcpp::Range(nburnin, niter - 1), Rcpp::_);

  return xt;
}

// [[Rcpp::export]]
arma::mat MCMCchainArma(const arma::vec& Xfull,
                        const arma::vec& beta,
                        const double g,
                        const double s,
                        const double a,
                        const double prob,
                        const double sd_g,
                        const double sd_s,
                        const size_t niter,
                        const size_t Jnumb,
                        const unsigned seed) {
  std::default_random_engine generator(seed);
  arma::mat xt(2 + beta.size(), niter);
  arma::colvec2 gs = {g, s};
  xt.col(0) = arma::join_cols(gs, BetaQuantileArma(g, s, a, prob, beta));
  const arma::vec X = Xfull.elem(arma::find(Xfull > a)) - a;
  const int n = X.size();
  const arma::mat UpperTriOnes = arma::trimatu(arma::ones(n, n));
  for(size_t j = 0; j < niter - 1; j++) {
    arma::colvec2 newpoint =
        MCMCnewpointArma(xt.at(0, j), xt.at(1, j), sd_g, sd_s, X, n, Jnumb,
                         generator, UpperTriOnes);
    xt.col(j + 1) = arma::join_cols(
        newpoint,
        BetaQuantileArma(newpoint.at(0), newpoint.at(1), a, prob, beta));
  }
  return xt.t();
}
