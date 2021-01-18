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

#ifndef MCMC_DISCRETE_HPP
#define MCMC_DISCRETE_HPP


#include <cmath>
#include <armadillo>
#include <cppbugs/mcmc.stochastic.hpp>

namespace cppbugs {

  template<typename T, typename U, typename Enable = void>
  class DiscreteLikelihiood;

  template<typename T, typename U>
  class DiscreteLikelihiood<T, U, typename std::enable_if<std::is_integral<typename std::remove_reference<T>::type>::value>::type> : public Likelihiood {
    const T& x_;
    const U p_;
  public:
    DiscreteLikelihiood(const T& x, const U& p): x_(x), p_(p) { }
    inline double calc() const {
      if(x_ < 0 || x_ >= (int)p_.n_elem)
        return -std::numeric_limits<double>::infinity();
      return cppbugs::log_approx(p_[x_]) - cppbugs::log_approx(arma::accu(p_));
    }
  };

  template<typename T, typename U>
  class DiscreteLikelihiood<T, U, typename std::enable_if<std::is_integral<typename std::remove_reference<T>::type::elem_type>::value>::type> : public Likelihiood {
    const T& x_;
    const U p_;
  public:
    DiscreteLikelihiood(const T& x, const U& p): x_(x), p_(p) { }
    inline double calc() const {
      if(!arma_all(x_ >= 0) || !arma_all(x_ < (int)p_.n_elem))
        return -std::numeric_limits<double>::infinity();
      double sum = 0;
      for(unsigned i = 0; i < x_.n_elem; i++)
        sum += cppbugs::log_approx(p_[x_[i]]);
      return sum - x_.n_elem * cppbugs::log_approx(arma::accu(p_));
    }
  };

  template<typename T>
  class Discrete : public DynamicStochastic<T> {
  public:
    Discrete(T value): DynamicStochastic<T>(value) {}

    template<typename U>
    Discrete<T>& ddiscr(/*const*/ U&& distr) {
      Stochastic::likelihood_functor = new DiscreteLikelihiood<T, U>(DynamicStochastic<T>::value,distr);
      return *this;
    }
  };

  template<typename T>
  class ObservedDiscrete : public Observed<T> {
  public:
    ObservedDiscrete(const T& value): Observed<T>(value) {}

    template<typename U>
    ObservedDiscrete<int>& ddiscr(/*const*/ U&& distr) {
      Stochastic::likelihood_functor = new DiscreteLikelihiood<T, U>(Observed<T>::value,distr);
      return *this;
    }
  };
} // namespace cppbugs
#endif // MCMC_BERNOULLI_HPP
