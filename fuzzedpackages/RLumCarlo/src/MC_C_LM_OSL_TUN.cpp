// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Title:   MC_C_CW_OSL_DELOC.cpp
// Author:  Johannes Friedrich, University of Bayreuth (Germany), Sebastian Kreutzer, Geography & Earth Sciences, Aberystwyth University (United Kingdom), based on
// equations provided by Vasilis Pagonis
// Contact: sebastian.kreutzer@aber.ac.uk
// Date:    Sun Feb 24 14:59:39 2019
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// [[Rcpp::depends(RcppArmadillo)]]
#include "RLumCarlo.h"
using namespace Rcpp;

// [[Rcpp::export("MC_C_LM_OSL_TUN")]]
List MC_C_LM_OSL_TUN(arma::vec times, int N_e, arma::vec r, double rho, double A) {

  //determine delta_t which allows to have delta t != 1
  double delta_t = calc_deltat(times);

  // set output matrices
  NumericMatrix signal (times.size(), r.size());
  NumericMatrix remaining_e (times.size(), r.size());
  NumericVector r_num;

    for(std::size_t k = 0; k < r.size(); ++k){

      std::size_t n_filled = N_e;

      for(std::size_t t = 0; t < times.size(); ++t){

        double P = A * (times[t]/max(times)) * exp(-(pow(rho,-1.0/3.0)) * r[k]);

        for(std::size_t j = 0; j < n_filled; ++j){

          //draw random number
          r_num = runif(1);

          if (r_num[0] < (P * delta_t))
            n_filled = n_filled - 1;

          if (n_filled == 0)
            break;

        } // end n_filled
        signal(t,k) = n_filled * P * 3 * pow((double)r[k],2.0) * exp(-(pow(r[k],3.0)));
        remaining_e(t,k) = n_filled;

        if (n_filled == 0)
          break;

      } // end t-loop

    } // end r-loop

    return(Rcpp::List::create(Named("signal") = signal,
                              Named("remaining_e") = remaining_e));
}

