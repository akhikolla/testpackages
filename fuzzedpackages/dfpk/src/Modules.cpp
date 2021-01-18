#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4cdf_reg_dtox_mod) {


    class_<rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> >("model_cdf_reg_dtox")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_cdf_reg_dtox_namespace::model_cdf_reg_dtox, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4cdf_reg_pktox_mod) {


    class_<rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> >("model_cdf_reg_pktox")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_cdf_reg_pktox_namespace::model_cdf_reg_pktox, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4logit_reg_pkcov_mod) {


    class_<rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> >("model_logit_reg_pkcov")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_logit_reg_pkcov_namespace::model_logit_reg_pkcov, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4logit_reg_pklogit_mod) {


    class_<rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> >("model_logit_reg_pklogit")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_logit_reg_pklogit_namespace::model_logit_reg_pklogit, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4logit_reg_pkpop_mod) {


    class_<rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> >("model_logit_reg_pkpop")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_logit_reg_pkpop_namespace::model_logit_reg_pkpop, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4reg_auc_mod) {


    class_<rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> >("model_reg_auc")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_reg_auc_namespace::model_reg_auc, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
