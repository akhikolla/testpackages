#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// draw from multivariate normal dist using Cholesky decomp. method
arma::mat mvrnormCpp(arma::vec mu, arma::mat Sigma) {
  int nCol = Sigma.n_cols;
  arma::mat u = arma::randn(1, nCol);
  arma::mat temp = mu.t() + u*arma::chol(Sigma);
  return temp.t();
}

// compute whether to accept or reject proposal based on
// log of acceptance rate
bool logAccept(double logR) {
  if (logR > 0.0) { //accept with prob min(r, 1)
    logR = 0.0;
  }
  double logU = log(R::runif(0, 1));
  if (logU < logR) {
    return(true);
  }
  else {
    return(false);
  }
}

// template to compute mean
template <class T>
inline double GetMean(T& x ) {
  return mean(x);
}

// compute row means
arma::vec RowMean(arma::mat m) {
  int nRow = m.n_rows;
  arma::vec out(nRow);
  for (int i = 0; i < nRow; i++) {
    arma::rowvec c = m.row(i);
    out(i) = GetMean(c);
  }
  return(out);
}

// compute and sum over quadratic form a_{i}'ba_{i}
double sumQFSq(arma::mat a, arma::mat b) {
  arma::vec sum = arma::zeros(1);
  int n = a.n_cols;
  for (unsigned i = 0; i < n; i++) {
    sum = sum + a.col(i).t()*b*a.col(i);
  }
  return(sum(0));
}

// compute and sum over quadratic form a_{i}'bc_{i}
double sumQFCon(arma::mat a, arma::mat b, arma::colvec c) {
  arma::vec sum = arma::zeros(1);
  int n = a.n_cols;
  for (unsigned i = 0; i < n; i++) {
    sum = sum + a.col(i).t()*b*c;
  }
  return(sum(0));
}

// compute log of the determinant using that the product of
// the eigenvalues of a matrix equals its determinant
double logDet(arma::mat Sigma) {
  return(sum(arma::log(arma::eig_sym(Sigma))));
}



//repeate the same value
arma::colvec rep(double a, int n) {
  arma::colvec A = arma::zeros(n);
  for (int i = 0; i < n; i++) {
    A(i) = a;
  }
  return(A);
}

//' Run MCMC to Fit Contour Model
//'
//' @param n_iter   number of iterations to run the MCMC
//' @param dists symmetric matrix of the same dimension as the number of
//'                lines being used, specifying distances among starting locations
//'                or angles.
//' @param x a matrix of observed distances (y) of dimension number
//'          of vectors by number of years
//' @param xU_vecs vector giving the indices of each x value vector in each year
//'               that is not observed (vector indices and year indices are
//'               paired, so must be ordered the same as xU_years)
//' @param xU_years vector giving the indices of each year in which each x value
//'                vector is not observed (vector indices and year indices are
//'                paired, so must be ordered the same as xU_vecs)
//' @param xU_prop_sd Standard deviation for proposals for xU
//' @param xU_lb Lower bounds for xU values being sampled (order must match
//'             orded of xU_vecs and xU_years)
//' @param xU_ub Upper bounds for xU values being sampled (order must match
//'             orded of xU_vecs and xU_years)
//' @param mu vector of the same length as the number of lines which specifies
//'           the values from which each element of \code{mu} will be initialized
//'           in the MCMC.
//' @param mu0 vector of the same length as the number of lines which specifies
//'            the prior mean for \code{mu}.
//' @param lambda0 matrix of the same dimension as the number of lines which
//'                specifices the prior covariance matrix for \code{mu}.
//' @param sigma vector of the same length as the number of lines which
//'              specifies the values from which each element in \code{sigma}
//'              will be initialized from
//' @param sigma_ind_1 vector giving the first index of each section of sigma's
//'                  to be sampled together
//' @param sigma_ind_2 vector giving the last index of each section of sigma's
//'                  to be sampled together
//' @param sigma_prop_cov covariance matrix of the same length as the number of
//'                     lines that is used in sampling \code{sigma} values
//' @param rho double between 0 and 1 from which the value of \code{rho} will
//'             be initialized
//' @param rho0_lb double between 0 and 1 which gives the lower bound of the
//'               uniform prior for \code{rho}
//' @param rho0_ub double between 0 and 1 which gives the upper bound of the
//'               uniform prior for \code{rho}.
//' @param rho_prop_sd standard deviation for the normal proposal distribution used
//'                  when proposing value for \code{rho} in the sampler. Defaults
//'                  to 0.01
//' @param sigma0_lb vector of the same length as the number of lines which
//'                 specifies the lower bound of the uniform prior for each
//'                 sigma value
//' @param sigma0_ub vector of the same length as the number of lines which
//'                 specifies the upper bound of the uniform prior for each sigma
//'                 value.
//' @param w Integer specifying how many samples of the parameters will be
//'           maintained. Samples from every wth iteration is stored.
//'
//' @return List of length 7 that gives the values of the MCMC chain for
//'         \code{xU}, \code{mu}, \code{sigma} and \code{rho} along with
//'         indicators of acceptance on each iteration: \code{xURate},
//'         \code{sigmaRate}, and \code{rhoRate}.
// [[Rcpp::export]]
List RunMCMC(int n_iter,  arma::mat dists,
             arma::mat x, arma::vec xU_vecs, arma::vec xU_years,
             arma::vec xU_prop_sd, arma::vec xU_lb,  arma::vec xU_ub,
             arma::colvec mu, arma::vec mu0,  arma::mat lambda0,
             arma::vec sigma, arma::uvec sigma_ind_1,
             arma::uvec sigma_ind_2, arma::mat sigma_prop_cov,
             double rho, double rho0_lb, double rho0_ub, double rho_prop_sd,
             arma::vec sigma0_lb, arma::vec sigma0_ub, int w) {

  //constants
  int nVecs = x.n_rows;
  double nYears = x.n_cols;
  int nXU = xU_lb.size();
  int nSigmaSec = sigma_ind_1.n_elem;

  //storage vectors and matrices
  arma::mat muStore(nVecs, n_iter/w);
  arma::mat sigmaStore(nVecs, n_iter/w);
  arma::vec rhoStore(n_iter/w);
  arma::mat xUStore(nXU, n_iter/w);

  //acceptance rate storage vectors and matrices
  arma::vec xURate = arma::zeros(nXU);
  arma::vec sigmaRate = arma::zeros(nSigmaSec);
  arma::vec rhoRate = arma::zeros(1);

  //derived initial values
  arma::mat Sigma(nVecs, nVecs);
  for (unsigned l = 0; l < nVecs; l++) {
    for (unsigned k = 0; k <= l; k++) {
      Sigma(l, k) = Sigma(k, l) = sigma(l)*sigma(k)*pow(rho, dists(l, k));
    }
  }
  arma::mat SigmaInv = inv(Sigma);
  arma::mat lambda0Inv = inv(lambda0);

  // // /////////main loop/////////////
  for (unsigned j = 0; j < n_iter; j++) {
    ///// Metropolis for xU's//////
    if (!xU_vecs.has_nan()) {
      for (unsigned i = 0; i < nXU; i++) {
        double xUProp = R::rnorm(x(xU_vecs(i), xU_years(i)), xU_prop_sd(i));
        bool acceptXU = false;
        //acceptance possible only if within prior bounds
        if ((xUProp > xU_lb(i)) && (xUProp < xU_ub(i))) {
          arma::vec xI = x.col(xU_years(i));
          arma::vec xPropI = xI; xPropI(xU_vecs(i)) = xUProp;
          arma::vec logRXU = -.5*xPropI.t()*SigmaInv*xPropI
                             + xPropI.t()*SigmaInv*mu
                             +.5*xI.t()*SigmaInv*xI
                             - xI.t()*SigmaInv*mu;
          acceptXU = logAccept(logRXU(0));
        }
        if (acceptXU) {
          xURate(i) = xURate(i) + 1;
          x(xU_vecs(i), xU_years(i)) = xUProp;
          if (j%w == 0) {
            xUStore(i, j/w) = xUProp;
          }
        } else {
          if (j%w == 0) {
            xUStore(i, j/w) =  x(xU_vecs(i), xU_years(i));
          }
        }
      }
    }

    //////////Gibbs for mu/////////////
    arma::vec xMean = RowMean(x);
    arma::mat lambdaN = inv(lambda0Inv + nYears*SigmaInv);
    arma::vec muN = (inv(lambda0Inv + nYears*SigmaInv)*
                    (lambda0Inv*mu0 + nYears*SigmaInv*xMean));
    arma::mat mu = mvrnormCpp(muN, lambdaN);
    if (j%w == 0) {
      muStore.col(j/w) = mu;
    }

    ///////Metropolis for sigma///////
    for (int g = 0; g < nSigmaSec; g++) {
      arma::uvec gInd = arma::linspace<arma::uvec>(sigma_ind_1(g), sigma_ind_2(g),
                                                   sigma_ind_2(g) - sigma_ind_1(g)
                                                   + 1);
      arma::mat gProp = mvrnormCpp(sigma(gInd), sigma_prop_cov(gInd, gInd));
      arma::vec sigmaGProp = sigma; sigmaGProp(gInd) = gProp;
      arma::mat SigmaGProp(nVecs, nVecs); arma::mat SigmaGPropInv(nVecs, nVecs);
      bool acceptSigmaG = false;
      //acceptance possible only if within prior bounds
      if (all(gProp > sigma0_lb(gInd)) && all(gProp < sigma0_ub(gInd))) {
        for (unsigned l = 0; l < nVecs; l++) {
          for (unsigned k = 0; k <= l; k++) {
            SigmaGProp(l, k) = SigmaGProp(k, l) = sigmaGProp(l)*sigmaGProp(k)*
                                                  pow(rho, dists(l, k));
          }
        }
        SigmaGPropInv = inv(SigmaGProp);
        double logRSigmaG = - (nYears/2)*logDet(SigmaGProp)
                            -.5*sumQFSq(x, SigmaGPropInv)
                            + sumQFCon(x, SigmaGPropInv, mu)
                            - (nYears/2)*sumQFSq(mu, SigmaGPropInv)
                            + (nYears/2)*logDet(Sigma)
                            + .5*sumQFSq(x, SigmaInv)
                            - sumQFCon(x, SigmaInv, mu)
                            + (nYears/2)*sumQFSq(mu, SigmaInv);
        acceptSigmaG = logAccept(logRSigmaG);
      }
      if (acceptSigmaG) {
        sigmaRate(g) = sigmaRate(g) + 1;
        sigma = sigmaGProp;
        Sigma = SigmaGProp;
        SigmaInv = SigmaGPropInv;
        if (j%w == 0) {
          sigmaStore(arma::span(sigma_ind_1(g), sigma_ind_2(g)), j/w) = sigmaGProp(gInd);
        }

      } else {
        if (j%w == 0) {
          sigmaStore(arma::span(sigma_ind_1(g), sigma_ind_2(g)), j/w) =  sigma(gInd);
        }
      }
    }

    /////////Metropolis for rho///////
    bool acceptRho = false;
    double rhoProp = R::rnorm(rho, rho_prop_sd);
    arma::mat SigmaRhoProp(nVecs, nVecs);
    arma::mat SigmaRhoPropInv(nVecs, nVecs);
    //acceptance possible only if within prior bounds
    if ((rhoProp > rho0_lb) && (rhoProp < rho0_ub)) {
      for (unsigned l = 0; l < nVecs; l++) {
        for (unsigned k = 0; k <= l; k++) {
          SigmaRhoProp(l, k) = SigmaRhoProp(k, l) = sigma(l)*sigma(k)*
                                                    pow(rhoProp, dists(l, k));
        }
      }
      SigmaRhoPropInv = inv(SigmaRhoProp);
      double logRRho = - (nYears/2)*logDet(SigmaRhoProp)
                       -.5*sumQFSq(x, SigmaRhoPropInv)
                       + sumQFCon(x, SigmaRhoPropInv, mu)
                       - (nYears/2)*sumQFSq(mu, SigmaRhoPropInv)
                       + (nYears/2)*logDet(Sigma)
                       + .5*sumQFSq(x, SigmaInv)
                       - sumQFCon(x, SigmaInv, mu)
                       + (nYears/2)*sumQFSq(mu, SigmaInv);
      acceptRho = logAccept(logRRho);
    }
    if (acceptRho) {
      rhoRate(0) = rhoRate(0) + 1;
      rho = rhoProp;
      Sigma = SigmaRhoProp;
      SigmaInv = SigmaRhoPropInv;
      if (j%w == 0) {
        rhoStore(j/w) = rhoProp;
      }
    } else {
      if (j%w == 0) {
        rhoStore(j/w) =  rho;
      }
    }
  }

  //return values
  List res;
  res["xU"] = xUStore; res["xU_rate"] = xURate/n_iter;
  res["xU_lb"] = xU_lb; res["xU_ub"] = xU_ub;
  res["mu"] = muStore;
  res["sigma"] = sigmaStore; res["sigma_rate"] = sigmaRate/n_iter;
  res["sigma_ind_1"] = sigma_ind_1; res["sigma_ind_2"] = sigma_ind_2;
  res["sigma0_lb"] = sigma0_lb; res["sigma0_ub"] = sigma0_ub;
  res["rho"] = rhoStore; res["rho_rate"] = rhoRate/n_iter;
  res["dists"] = dists;
  res["w"] = w;
  return(res);
}

