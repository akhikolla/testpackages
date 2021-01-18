#include "dispatcher.h"

SEXP dispatcher(SEXP sexp_membership_name,
                SEXP sexp_membership_init_from_R,
                SEXP sexp_model_name,
                SEXP sexp_network_from_R,
                SEXP sexp_real_EM)
{
    std::string membership_name = Rcpp::as<std::string>(sexp_membership_name);
    Rcpp::List membership_init_from_R(sexp_membership_init_from_R);
    std::string model_name = Rcpp::as<std::string>(sexp_model_name);
    Rcpp::List network_from_R(sexp_network_from_R);
    bool real_EM = Rcpp::as<bool>(sexp_real_EM);

    if(real_EM)
        return distpatcher_membership_model<true>(membership_name,
                                                  membership_init_from_R,
                                                  model_name,
                                                  network_from_R);
    else
        return distpatcher_membership_model<false>(membership_name,
                                                  membership_init_from_R,
                                                  model_name,
                                                  network_from_R);
}
    

