/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2020
 *     Darjus Hosszejni Copyright (C) 2019-2020
 *  
 *  This file is part of the R package stochvol: Efficient Bayesian
 *  Inference for Stochastic Volatility Models.
 *  
 *  The R package stochvol is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *  
 *  The R package stochvol is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with the R package stochvol. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */

/*
 * utils_latent_states.cc
 * 
 * Definitions of the functions declared in utils_latent_states.h.
 * Documentation can also be found in utils_latent_states.h.
 */

#include <RcppArmadillo.h>
#include "utils.h"
#include "utils_latent_states.h"
#include "densities.h"
#include <cmath>

namespace stochvol {

namespace fast_sv {

CholeskyTridiagonal cholesky_tridiagonal(
    const arma::vec& omega_diag,
    const double omega_offdiag) {
  const int T = omega_diag.n_elem - 1;
  arma::vec chol_diag(T+1);
  arma::vec chol_offdiag(T+1);
  chol_diag[0] = std::sqrt(omega_diag[0]);
  for (int j = 1; j < T+1; j++) {
    chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
    chol_diag[j] = std::sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
  }
  return {std::move(chol_diag), std::move(chol_offdiag)};
}

arma::vec forward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector) {
  const int T = chol_diag.n_elem - 1;
  arma::vec htmp(T+1);
  htmp[0] = covector[0]/chol_diag[0];
  for (int j = 1; j < T+1; j++) {
    htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
  }
  return htmp;
}

arma::vec backward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp) {
  const int T = chol_diag.size() - 1;
  arma::vec h(T+1);
  h[T] = htmp[T] / chol_diag[T];
  for (int j = T-1; j >= 0; j--) {
    h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
  }
  return h;
}

arma::uvec inverse_transform_sampling(
    const arma::vec& mixprob,
    const int T) {
  arma::uvec r(T);
  for (int j = 0; j < T; j++) {
    int index = (10-1)/2;  // start searching in the middle
    const double unnorm_cdf_value = R::unif_rand()*mixprob[9 + 10*j];  // current (non-normalized) value
    bool larger = false;  // indicates that we already went up
    bool smaller = false; // indicates that we already went down
    while(true) {
      if (unnorm_cdf_value > mixprob[index +  10*j]) {
        index++;
        if (smaller) {
          break;
        } else {
          larger = true;
        }
      } else if (larger || index == 0) {
        break;
      } else {
        index--;
        smaller = true;
      }
    }
    r[j] = index;
  }
  return r;
}

arma::vec find_mixture_indicator_cdf(
    const arma::vec& datanorm)  {
  const int T = datanorm.n_elem;
  arma::vec mixprob(10 * T);
  for (int j = 0; j < T; j++) {  // TODO slow (10*T calls to exp)!
    const int first_index = 10*j;
    mixprob[first_index] = std::exp(mix_pre[0]-(datanorm[j]-mix_mean[0])*(datanorm[j]-mix_mean[0])*mix_2varinv[0]);
    for (int r = 1; r < 10; r++) {
      mixprob[first_index+r] = mixprob[first_index+r-1] + std::exp(mix_pre[r]-(datanorm[j]-mix_mean[r])*(datanorm[j]-mix_mean[r])*mix_2varinv[r]);
    }
  }
  return mixprob;
}

}  // END namespace fast_sv

namespace general_sv {

double h_log_posterior(
    const arma::vec& h,  // centered
    const arma::vec& y,
    const double phi,
    const double rho,
    const double sigma,
    const double mu,
    const double h0) {
  const double rho_const = std::sqrt(1 - rho * rho);
  const int n = y.size();
  const arma::vec exp_h_half = arma::exp(0.5 * h);  // TODO cached?
  double result = logdnorm2(h[0], mu + phi * (h0 - mu), sigma);  // log p(h_1 | theta)
  for (int t = 0; t < n - 1; t++) {  // TODO parallel?
    result += logdnorm2(h[t + 1], mu + phi * (h[t] - mu), sigma);
    result += logdnorm2(y[t], exp_h_half[t] * rho * (h[t + 1] - mu - phi * (h[t] - mu)) / sigma, exp_h_half[t] * rho_const, .5 * h[t]);
  }
  result += logdnorm2(y[n - 1], 0, exp_h_half[n - 1], .5 * h[n - 1]);
  return result;
}

double h_aux_log_posterior(
    const arma::vec& h,  // centered
    const arma::vec& y_star,
    const arma::ivec& d,
    const double phi,
    const double rho,
    const double sigma,
    const double mu,
    const double h0) {
  const int n = y_star.size();
  static const int mix_count = mix_a.n_elem;
  const double sigma2 = std::pow(sigma, 2);
  static const arma::vec::fixed<10> exp_m_half = arma::exp(mix_mean * .5);
  const arma::vec C = rho * sigma * exp_m_half;  // re-used constant

  double result = logdnorm2(h[0], mu + phi * (h0 - mu), sigma);  // log p(h_1 | theta)
  for (int t = 0; t < n; t++) {  // TODO parallel?
    double subresult = 0;
    if (t < n - 1) {
      for (int j = 0; j < mix_count; j++) {
        const double h_mean = mu + phi * (h[t] - mu) + d[t] * C[j] * mix_a[j];
        const double h_var = mix_var[j] * std::pow(C[j] * mix_b[j], 2) + sigma2 * (1 - rho * rho);
        const double yh_cov = d[t] * C[j] * mix_b[j] * mix_var[j];
        const double y_mean = h[t] + mix_mean[j] + yh_cov / h_var * (h[t + 1] - h_mean);
        const double y_var = (1 - std::pow(yh_cov, 2) / (mix_var[j] * h_var)) * mix_var[j];
        subresult += std::exp(-.5 * (std::pow(y_star[t] - y_mean, 2) / y_var + std::pow(h[t + 1] - h_mean, 2) / h_var)) / std::sqrt(y_var * h_var) * mix_prob[j];
      }
    } else {
      for (int j = 0; j < mix_count; j++) {
        subresult += std::exp(-.5 * std::pow(y_star[t] - (h[t] + mix_mean[j]), 2) / mix_var[j]) / std::sqrt(mix_var[j]) * mix_prob[j];
      }
    }
    result += std::log(subresult);
  }

  return result;
}

}  // END namespace general_sv

}

