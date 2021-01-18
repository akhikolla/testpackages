.onLoad <- function(libname, pkgname) { Rcpp::loadModule("stan_fit4cdf_reg_dtox_mod", TRUE)
										Rcpp::loadModule("stan_fit4cdf_reg_pktox_mod", TRUE)
										Rcpp::loadModule("stan_fit4logit_reg_pkcov_mod", TRUE)
										Rcpp::loadModule("stan_fit4logit_reg_pklogit_mod", TRUE)
										Rcpp::loadModule("stan_fit4logit_reg_pkpop_mod", TRUE)
										Rcpp::loadModule("stan_fit4reg_auc_mod", TRUE) 
}
