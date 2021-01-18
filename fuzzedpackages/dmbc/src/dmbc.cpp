// dmbc.cpp

#include "dmbc.h"

// Note: RcppExport is an alias for extern "C"

//' Internal functions for MCMC simulation.
//'
//' For internal use only.
//'
//' @param raiD internal SEXP data structure
//' @param raix internal SEXP data structure
//' @param raing internal SEXP data structure
//' @param radalpha internal SEXP data structure
//' @param rn internal SEXP data structure
//' @param rp internal SEXP data structure
//' @param rG internal SEXP data structure
//' @param rS internal SEXP data structure
//' @param rtotiter internal SEXP data structure
//' @param radZ internal SEXP data structure
//' @param rgamma_z internal SEXP data structure
//' @param reta internal SEXP data structure
//' @param rgamma_alpha internal SEXP data structure
//' @param rsigma2 internal SEXP data structure
//' @param rlambda internal SEXP data structure
//' @param rhyper_eta_a internal SEXP data structure
//' @param rhyper_eta_b internal SEXP data structure
//' @param rhyper_sigma2_a internal SEXP data structure
//' @param rhyper_sigma2_b internal SEXP data structure
//' @param rhyper_lambda internal SEXP data structure
//' @param rfamily internal SEXP data structure
//' @param rverbose internal SEXP data structure
//'
//' @aliases dmbc-internal
//' @aliases dmbc_internal
//'
// [[Rcpp::export]]
RcppExport SEXP dmbc_mcmc(
  SEXP raiD,
  SEXP raix,
  SEXP raing,
  SEXP radalpha,
  SEXP rn,
  SEXP rp,
  SEXP rG,
  SEXP rS,
  SEXP rtotiter,
  SEXP radZ,
  SEXP rgamma_z,
  SEXP reta,
  SEXP rgamma_alpha,
  SEXP rsigma2,
  SEXP rlambda,
  SEXP rhyper_eta_a,
  SEXP rhyper_eta_b,
  SEXP rhyper_sigma2_a,
  SEXP rhyper_sigma2_b,
  SEXP rhyper_lambda,
  SEXP rfamily,
  SEXP rverbose){
  SEXP rAns = NULL;
  int rAnsItems = 12;

  int n = Rf_asInteger(rn);
  int p = Rf_asInteger(rp);
  int G = Rf_asInteger(rG);
  int S = Rf_asInteger(rS);
  int totiter = Rf_asInteger(rtotiter);
  double gamma_z = Rf_asReal(rgamma_z);
  double gamma_alpha = Rf_asReal(rgamma_alpha);
  double hyper_sigma2_a = Rf_asReal(rhyper_sigma2_a);
  double hyper_sigma2_b = Rf_asReal(rhyper_sigma2_b);
  int verbose = INTEGER(rverbose)[0];
  // const char* family = CHAR(STRING_ELT(rfamily, 0)); // [for future developments]

  SEXP rz_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*n*p*G)));
  SEXP ralpha_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  SEXP reta_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  SEXP rsigma2_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  SEXP rlambda_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  SEXP rprob_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*S*G)));
  SEXP rx_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*S)));
  SEXP rx_ind_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*S*G)));
  SEXP raccept = PROTECT(Rf_allocVector(REALSXP, (2*G)));
  SEXP rloglik = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogprior = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogpost = PROTECT(Rf_allocVector(REALSXP, totiter));

  dmbc_mcmc_binom(REAL(rz_chain), REAL(ralpha_chain), REAL(reta_chain), REAL(rsigma2_chain), REAL(rlambda_chain),
    REAL(rprob_chain), REAL(rx_chain), REAL(rx_ind_chain), REAL(raccept), REAL(rloglik), REAL(rlogprior),
    REAL(rlogpost), INTEGER(raiD), REAL(radZ), INTEGER(raix), INTEGER(raing), REAL(radalpha), REAL(reta),
    REAL(rsigma2), REAL(rlambda), REAL(rhyper_eta_a), REAL(rhyper_eta_b), REAL(rhyper_lambda), gamma_z,
    gamma_alpha, hyper_sigma2_a, hyper_sigma2_b, totiter, n, p, S, G, verbose);

  // packing results
  PROTECT(rAns = Rf_allocVector(VECSXP, rAnsItems));

  SET_VECTOR_ELT(rAns, 0, rz_chain);
  SET_VECTOR_ELT(rAns, 1, ralpha_chain);
  SET_VECTOR_ELT(rAns, 2, reta_chain);
  SET_VECTOR_ELT(rAns, 3, rsigma2_chain);
  SET_VECTOR_ELT(rAns, 4, rlambda_chain);
  SET_VECTOR_ELT(rAns, 5, rprob_chain);
  SET_VECTOR_ELT(rAns, 6, rx_chain);
  SET_VECTOR_ELT(rAns, 7, rx_ind_chain);
  SET_VECTOR_ELT(rAns, 8, raccept);
  SET_VECTOR_ELT(rAns, 9, rloglik);
  SET_VECTOR_ELT(rAns, 10, rlogprior);
  SET_VECTOR_ELT(rAns, 11, rlogpost);

  // cleanup and return
  UNPROTECT(13);  // rAns, rz_chain, ralpha_chain, reta_chain, rsigma2_chain, rlambda_chain, rprob_chain,
                  // rx_ind_chain, rx_chain, raccept, rloglik, rlogprior, rlogpost

  return rAns;
}

//' Function for relabeling the parameter chain
//'
//' @param radtheta internal SEXP data structure
//' @param radz internal SEXP data structure
//' @param radeta internal SEXP data structure
//' @param radsigma2 internal SEXP data structure
//' @param radlambda internal SEXP data structure
//' @param radprob internal SEXP data structure
//' @param raix_ind internal SEXP data structure
//' @param rinit internal SEXP data structure
//' @param rM internal SEXP data structure
//' @param rR internal SEXP data structure
//'
//' @aliases dmbc-internal
//' @aliases dmbc_internal
//'
//' @rdname dmbc_mcmc
//'
// [[Rcpp::export]]
RcppExport SEXP dmbc_relabel(
  SEXP radtheta,
  SEXP radz,
  SEXP radalpha,
  SEXP radeta,
  SEXP radsigma2,
  SEXP radlambda,
  SEXP radprob,
  SEXP raix_ind,
  SEXP rinit,
  SEXP rn,
  SEXP rp,
  SEXP rS,
  SEXP rM,
  SEXP rR,
  SEXP rG,
  SEXP rverbose){
  SEXP rAns = NULL;
  const int rAnsItems = 8;

  int init = Rf_asInteger(rinit);
  int n = Rf_asInteger(rn);
  int p = Rf_asInteger(rp);
  int S = Rf_asInteger(rS);
  int M = Rf_asInteger(rM);
  int R = Rf_asInteger(rR);
  int G = Rf_asInteger(rG);
  int verbose = INTEGER(rverbose)[0];

  relabel_celeux(REAL(radtheta), REAL(radz), REAL(radalpha), REAL(radeta), REAL(radsigma2), REAL(radlambda),
      REAL(radprob), INTEGER(raix_ind), init, n, p, S, M, R, G, verbose);

  // packing results
  PROTECT(rAns = Rf_allocVector(VECSXP, rAnsItems));

  SET_VECTOR_ELT(rAns, 0, radtheta);
  SET_VECTOR_ELT(rAns, 1, radz);
  SET_VECTOR_ELT(rAns, 2, radalpha);
  SET_VECTOR_ELT(rAns, 3, radeta);
  SET_VECTOR_ELT(rAns, 4, radsigma2);
  SET_VECTOR_ELT(rAns, 5, radlambda);
  SET_VECTOR_ELT(rAns, 6, radprob);
  SET_VECTOR_ELT(rAns, 7, raix_ind);

  // cleanup and return
  UNPROTECT(1);  // rAns

  return rAns;
}

//' Packing of the parameter chain to be run before relabeling
//'
//' @rdname dmbc_mcmc
//'
// [[Rcpp::export]]
RcppExport SEXP dmbc_pack_par(
  SEXP radz,
  SEXP radalpha,
  SEXP radlambda,
  SEXP rn,
  SEXP rp,
  SEXP rM,
  SEXP rG){
  int n = Rf_asInteger(rn);
  int p = Rf_asInteger(rp);
  int M = Rf_asInteger(rM);
  int G = Rf_asInteger(rG);
  int r = n*(n - 1)/2;

  SEXP rAns = PROTECT(Rf_allocVector(REALSXP, M*(r + 1)*G));

  pack_par(REAL(rAns), REAL(radz), REAL(radalpha), REAL(radlambda), n, p, M, G);

  // cleanup and return
  UNPROTECT(1);  // rAns

  return rAns;
}
