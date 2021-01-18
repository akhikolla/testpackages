///////////////////////////////////////////////////////////////////////////
// Copyright (C) 2011 Whit Armstrong                                     //
//                                                                       //
// This program is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published by  //
// the Free Software Foundation, either version 3 of the License, or     //
// (at your option) any later version.                                   //
//                                                                       //
// This program is distributed in the hope that it will be useful,       //
// but WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
// GNU General Public License for more details.                          //
//                                                                       //
// You should have received a copy of the GNU General Public License     //
// along with this program.  If not, see <http://www.gnu.org/licenses/>. //
///////////////////////////////////////////////////////////////////////////

#ifndef MCMC_MATH_HPP
#define MCMC_MATH_HPP

#include <stdexcept>
#include <cmath>
#include <armadillo>

namespace cppbugs {

  inline float log_approx(float val) {
    union { float f; int i; } valu;
    valu.f = val;
    float exp = valu.i >> 23;

    float addcst = val > 0 ? -89.76031805111f : -INFINITY;

    valu.i = (valu.i & 0x7FFFFF) | 0x3F800000;

    float x = valu.f;

    return x*(2.79224171f+x*(-1.44246914f+x*(0.43585738f+x*-0.05486225f)))
      + (addcst + 0.69314718055995f*exp);
  }

  inline double log_approx(const double x) {
    return cppbugs::log_approx((float)x);
  }

  inline float log_approx(const int x) {
    return cppbugs::log_approx((float)x);
  }

  inline arma::vec log_approx(const arma::vec& x) {
    arma::vec res(size(x));
    for(int i = 0; i < x.size(); i++) {
      res[i] = cppbugs::log_approx(x[i]);
    }
    return res;
  }

  // We do not inline these constants, because that makes GCC
  // not able to recognize max/min pattern, and then the code is
  // not vectorized.
  float exp_cst1 = 2139095040.f;
  float exp_cst2 = 0.f;

  inline float exp_approx(float val) {
    union { int i; float f; } xu;
    float val2 = 12102203.1615614f*val+1065353216.f;
    float val3 = val2 < exp_cst1 ? val2 : exp_cst1;
    float val4 = val3 > exp_cst2 ? val3 : exp_cst2;
    int val4i = (int) val4;
    xu.i = val4i & 0x7F800000;
    union { int i; float f; } xu2;
    xu2.i = (val4i & 0x7FFFFF) | 0x3F800000;
    float b = xu2.f;
    float ret = xu.f * (0.51079604f+b*(0.30980503f+b*(0.16876894f+b*(-0.00303925f+b*0.01367652f))));
    return ret;
  }

  inline arma::vec exp_approx(const arma::vec& x) {
    arma::vec res(size(x));
    for(int i = 0; i < x.size(); i++) {
      res[i] = cppbugs::exp_approx(x[i]);
    }
    return res;
  }

  // factln
  inline double factln(const int i) {
    static std::vector<double> factln_table;

    if(i < 0) {
      return -std::numeric_limits<double>::infinity();
    }

    if(i > 100) {
      return ::lgamma(double(i) + 1);
    }

    if(factln_table.size() < static_cast<size_t>(i+1)) {
      for(int j = factln_table.size(); j < (i+1); j++) {
        factln_table.push_back(::lgamma(double(j) + 1));
      }
    }
    return factln_table[i];
  }

  template<typename T>
  inline bool arma_all(const T& x) { return arma::all(x); }
  inline bool arma_all(const bool x) { return x; }

  // Stochastic/Math related functions
  template<typename T, typename U>
  inline auto schur_product(T&& x, U&& y) ->
    decltype (arma::operator%(std::forward<T>(x), std::forward<U>(y))){
    return arma::operator%(std::forward<T>(x), std::forward<U>(y));
  }

  template<typename T>
  inline auto schur_product(const typename T::elem_type& x, T&& y) ->
    decltype (arma::operator*(x, std::forward<T>(y))){
    return arma::operator*(x, std::forward<T>(y));
  }

  template<typename T>
  inline auto schur_product(T&& x, const typename T::elem_type& y) ->
    decltype (arma::operator*(std::forward<T>(x), y)){
    return arma::operator*(std::forward<T>(x), y);
  }

  double schur_product(const int x, const double y) { return x * y; }
  double schur_product(const double x, const int y) { return x * y; }
  double schur_product(const double x, const double y) { return x * y; }
  double schur_product(const float x, const double y) { return x * y; }
  double schur_product(const double x, const float y) { return x * y; }
  float schur_product(const int x, const float y) { return x * y; }
  float schur_product(const float x, const int y) { return x * y; }
  float schur_product(const float x, const float y) { return x * y; }
  int schur_product(const int x, const int y) { return x * y; }

  double dim_size(const double) {
    return 1;
  }

  double dim_size(const float) {
    return 1;
  }

  double dim_size(const int) {
    return 1;
  }

  double dim_size(const bool) {
    return 1;
  }

  template<typename T>
  double dim_size(const T& x) {
    return x.n_elem;
  }

  static inline double square(double x) {
    return x*x;
  }

  static inline float square(float x) {
    return x*x;
  }

  static inline int square(int x) {
    return x*x;
  }

  double mahalanobis(const arma::vec& x, const arma::vec& mu, const arma::mat& sigma) {
    const arma::vec err = x - mu;
    return arma::as_scalar(err.t() * sigma.i() * err);
  }

  double mahalanobis(const arma::rowvec& x, const arma::rowvec& mu, const arma::mat& sigma) {
    const arma::rowvec err = x - mu;
    return arma::as_scalar(err * sigma.i() * err.t());
  }

  template<typename T, typename U, typename V>
  double normal_logp(const T& x, const U& mu, const V& tau) {
    return arma::accu(0.5f*cppbugs::log_approx(0.5f*tau/M_PI)
                      - 0.5f * schur_product(tau, square(x - mu)));
  }

  template<typename T, typename U, typename V>
  double uniform_logp(const T& x, const U& lower, const V& upper) {
    if(!arma_all(x > lower) || !arma_all(x < upper))
      return -std::numeric_limits<double>::infinity();
    return -arma::accu(cppbugs::log_approx(upper - lower));
  }

  template<typename T, typename U, typename V>
  double gamma_logp(const T& x, const U& alpha, const V& beta) {
    if(!arma_all(x > 0))
      return -std::numeric_limits<double>::infinity();
    return
      arma::accu(schur_product((alpha - 1.0f),cppbugs::log_approx(x))
                 - schur_product(beta,x) - lgamma(alpha)
                 + schur_product(alpha,cppbugs::log_approx(beta)));
  }

  template<typename T, typename U, typename V>
  double beta_logp(const T& x, const U& alpha, const V& beta) {
    if(!arma_all(x > 0) || !arma_all(x < 1) ||
       !arma_all(alpha > 0) || !arma_all(beta > 0))
      return -std::numeric_limits<double>::infinity();
    return arma::accu(lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta)
                      + schur_product(alpha - 1.0f, cppbugs::log_approx(x))
                      + schur_product(beta - 1.0f, cppbugs::log_approx(1.0f - x)));
  }

  template<typename T, typename U, typename V>
  double binom_logp(const T& x, const U& n, const V& p) {
    if(!arma_all(x >= 0) || !arma_all(x <= n))
      return -std::numeric_limits<double>::infinity();
    return arma::accu(schur_product(x,cppbugs::log_approx(p))
                      + schur_product((n-x),cppbugs::log_approx(1-p)) +
                                            factln(n) - factln(x) - factln(n-x));
  }

  template<typename T, typename U>
  double bernoulli_logp(const T& x, const U& p) {
    if(!arma_all(x >= 0) || !arma_all(x <= 1))
      return -std::numeric_limits<double>::infinity();
    return arma::accu(schur_product(x, cppbugs::log_approx(p))
                      + schur_product((1-x), cppbugs::log_approx(1-p)));
  }

  // sigma denotes cov matrix rather than precision matrix
  double multivariate_normal_sigma_logp(const arma::rowvec& x, const arma::rowvec& mu, const arma::mat& sigma) {
    const double log_2pi = log(2 * M_PI);
    arma::mat R(arma::zeros<arma::mat>(sigma.n_cols,sigma.n_cols));

    // non-positive definite test via chol
    if(chol(R,sigma) == false) { return -std::numeric_limits<double>::infinity(); }

    // otherwise calc logp
    return -(x.n_elem * log_2pi + cppbugs::log_approx(arma::det(sigma)) +
             mahalanobis(x,mu,sigma))/2;
  }

  // sigma denotes cov matrix rather than precision matrix
  double multivariate_normal_sigma_logp(const arma::vec& x, const arma::vec& mu, const arma::mat& sigma) {
    const double log_2pi = log(2 * M_PI);
    arma::mat R(arma::zeros<arma::mat>(sigma.n_cols,sigma.n_cols));

    // non-positive definite test via chol
    if(chol(R,sigma) == false) { return -std::numeric_limits<double>::infinity(); }

    // otherwise calc logp
    return -(x.n_elem * log_2pi + cppbugs::log_approx(arma::det(sigma)) +
             mahalanobis(x,mu,sigma))/2;
  }

  template<typename T, typename U, typename V>
  void dimension_check(const T& x, const U& hyper1, const V& hyper2) {
    if(dim_size(hyper1) > dim_size(x) || dim_size(hyper2) > dim_size(x)) {
      throw std::logic_error("ERROR: dimensions of hyperparmeters are larger than the stochastic variable itself (is this really what you wanted to do?)");
    }
  }

  template<typename T, typename U>
  void dimension_check(const T& x, const U& hyper1) {
    if(dim_size(hyper1) > dim_size(x)) {
      throw std::logic_error("ERROR: dimensions of hyperparmeters are larger than the stochastic variable itself (is this really what you wanted to do?)");
    }
  }

} // namespace cppbugs
#endif // MCMC_STOCHASTIC_HPP
