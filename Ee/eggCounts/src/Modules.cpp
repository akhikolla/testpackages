#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4indefficacy_mod) {


    class_<rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> >("model_indefficacy")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_indefficacy_namespace::model_indefficacy, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4nb_mod) {


    class_<rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> >("model_nb")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_nb_namespace::model_nb, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4paired_mod) {


    class_<rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> >("model_paired")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_paired_namespace::model_paired, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4simple_mod) {


    class_<rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> >("model_simple")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_simple_namespace::model_simple, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4unpaired_mod) {


    class_<rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> >("model_unpaired")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_unpaired_namespace::model_unpaired, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4zinb_mod) {


    class_<rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> >("model_zinb")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_zinb_namespace::model_zinb, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4zipaired_mod) {


    class_<rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> >("model_zipaired")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_zipaired_namespace::model_zipaired, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4ziunpaired_mod) {


    class_<rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> >("model_ziunpaired")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_ziunpaired_namespace::model_ziunpaired, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
