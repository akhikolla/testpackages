
#include "conf.h"

#include <RcppArmadillo.h>
#include <string>

using namespace arma;

#include "misc_functions.h"
#include "memberships/sbm.h"
#include "memberships/lbm.h"
#include "generic_functions.h"
#include "model_base.h"
#include "general_em.h"
#include "models/naive_bernoulli.h"
#include "models/bernoulli.h"
#include "models/poisson.h"
#include "models/poisson_covariates.h"
#include "models/gaussian.h"
#include "models/gaussian_multivariate.h"
#include "models/gaussian_multivariate_independent.h"
#include "models/gaussian_multivariate_independent_homoscedastic.h"
#include "models/bernoulli_covariates.h"
#include "models/bernoulli_covariates_fast.h"
#include "models/bernoulli_multiplex.h"
#include "models/gaussian_covariates.h"


template<class membership_type, class model_type, class network_type, bool real_EM>
Rcpp::List estim(membership_type & membership_init,
           Rcpp::List & network_from_R)
{
    network_type net(network_from_R);

    result<membership_type, model_type> res =
        em<membership_type, model_type, network_type, real_EM>( membership_init, net );

    return res.export_to_R();
}

template<class membership_type, bool real_EM>
Rcpp::List dispatcher_model(membership_type & membership_init,
                      std::string & model_name,
                      Rcpp::List & network_from_R)
{

    /* Each model must have a entry here */

    if(model_name == "naive_bernoulli")
        return estim<membership_type, naive_bernoulli, naive_bernoulli::network, real_EM>(
                membership_init, network_from_R);

    if(model_name == "bernoulli")
        return estim<membership_type, bernoulli, bernoulli::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "bernoulli_covariates")
        return estim<membership_type, bernoulli_covariates, bernoulli_covariates::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "bernoulli_covariates_fast")
        return estim<membership_type, bernoulli_covariates_fast, bernoulli_covariates_fast::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "bernoulli_multiplex")
        return estim<membership_type, bernoulli_multiplex, bernoulli_multiplex::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "poisson")
        return estim<membership_type, poisson, poisson::network, real_EM>(
                membership_init, network_from_R);

    if(model_name == "poisson_covariates")
        return estim<membership_type, poisson_covariates, poisson_covariates::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "gaussian")
        return estim<membership_type, gaussian, gaussian::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "gaussian_multivariate")
        return estim<membership_type, gaussian_multivariate, gaussian_multivariate::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "gaussian_multivariate_independent")
        return estim<membership_type, gaussian_multivariate_independent, gaussian_multivariate_independent::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "gaussian_multivariate_independent_homoscedastic")
        return estim<membership_type, gaussian_multivariate_independent_homoscedastic, gaussian_multivariate_independent_homoscedastic::network, real_EM>(
                membership_init, network_from_R);
    
    if(model_name == "gaussian_covariates")
        return estim<membership_type, gaussian_covariates, gaussian_covariates::network, real_EM>(
                membership_init, network_from_R);


    // What we are doing here ? Nothing.
    return Rcpp::List();
}

template<class membership_type, bool real_EM>
Rcpp::List init_membership_and_dispatcher_model(Rcpp::List & membership_from_R,
                                          std::string & model_name,
                                          Rcpp::List & network_from_R)
{
    membership_type membership_init(membership_from_R);
    return dispatcher_model<membership_type,real_EM>(
            membership_init, model_name, network_from_R);
}

template<bool real_EM>
Rcpp::List distpatcher_membership_model(std::string & membership_name,
                                  Rcpp::List & membership_init_from_R,
                                  std::string & model_name,
                                  Rcpp::List & network_from_R)
{
    /* Each membership must have a entry here */

    if(membership_name == "SBM")
        return init_membership_and_dispatcher_model<SBM,real_EM>(membership_init_from_R,
                                                                 model_name,
                                                                 network_from_R);
    
    if(membership_name == "SBM_sym")
        return init_membership_and_dispatcher_model<SBM_sym,real_EM>(membership_init_from_R,
                                                                 model_name,
                                                                 network_from_R);
    
    if(membership_name == "LBM")
        return init_membership_and_dispatcher_model<LBM,real_EM>(membership_init_from_R,
                                                                 model_name,
                                                                 network_from_R);

    // What we are doing here ? Nothing.
    return Rcpp::List();

}
                                          

RcppExport
SEXP dispatcher(SEXP membership_name,
                SEXP membership_init_from_R,
                SEXP model_name,
                SEXP network_from_R,
                SEXP real_EM);

