///////////////////////////////////////////////////////////////////////////
// Copyright (C) 2011 Jacques-Henri Jourdan                              //
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

#ifndef MCMC_EXPONENTIAL_CENSORED_HPP
#define MCMC_EXPONENTIAL_CENSORED_HPP


#include <cmath>
#include <armadillo>
#include <cppbugs/mcmc.stochastic.hpp>
#include <cstdio>
namespace cppbugs {

  template <typename T,typename U, typename V>
  class ExponentialCensoredLikelihiood : public Likelihiood {
    const T& x_;
    const U lambda_;
    const V delta_;
  public:
    ExponentialCensoredLikelihiood(const T& x, const U& lambda, const V& delta): x_(x), lambda_(lambda), delta_(delta) { dimension_check(x_, lambda_, delta_); }
    inline double calc() const {
      if(!arma_all(x_ > 0))
        return -std::numeric_limits<double>::infinity();
      return arma::accu(schur_product(delta_, cppbugs::log_approx(lambda_)) - schur_product(lambda_, x_));
    }
  };

  template<typename T>
  class ObservedExponentialCensored : public Observed<T> {
  public:
    ObservedExponentialCensored(const T& value): Observed<T>(value) {}

    template<typename U>
    ObservedExponentialCensored<T>& dexpcens(/*const*/ U&& lambda) {
      Stochastic::likelihood_functor =
        new ExponentialCensoredLikelihiood<typename std::remove_reference<T>::type::first_type, U, typename std::remove_reference<T>::type::second_type>
           (Observed<T>::value.first,
            lambda,
            Observed<T>::value.second);
      return *this;
    }
  };

} // namespace cppbugs
#endif // MCMC_EXPONENTIAL_CENSORED_HPP
