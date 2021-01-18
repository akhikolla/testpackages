#ifndef _SpatialTools_UTILS_H
#define _SpatialTools_UTILS_H

#include <RcppArmadillo.h>

arma::mat rcondsim(int nsim, arma::vec y, arma::mat w, arma::mat Vediag, arma::mat dV, int method, double m);
arma::mat decomp_V(const arma::mat &V, int method);
arma::mat rmvnorm(int nsim, const arma::mat &mu, const arma::mat &V, int method);
arma::mat rcondnorm(int nsim, const arma::mat &y,
                    const arma::mat &mu, const arma::mat &mup,
                    const arma::mat &V, const arma::mat &Vp, const arma::mat &Vop, int method);
arma::mat dist1(const arma::mat &coords);
arma::mat dist2(const arma::mat &coords, const arma::mat &pcoords);
arma::mat cov_spBayes(const arma::mat &D, int sp_type, double sigmasq,
                      double phi, double nu, double ev, double fv);

#endif
