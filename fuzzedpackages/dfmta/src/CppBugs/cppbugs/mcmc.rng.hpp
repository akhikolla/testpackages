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

#ifndef MCMC_RNG_HPP
#define MCMC_RNG_HPP

#include <random>
#include <cppbugs/mcmc.rng.base.hpp>

namespace cppbugs {

  template<typename T>
  class SpecializedRng : public RngBase {
  private:
    T generator_;
    std::uniform_real_distribution<double> uniform_rng_;
    double next_norm_;
  public:
    SpecializedRng(long seed): RngBase(),
                               generator_(seed),
                               uniform_rng_(0, 1) {
      next_norm_ = NAN;
    }

    double normal() {
      if(next_norm_ != next_norm_) {
        double x, y, s;
        do {
          x = uniform_rng_(generator_)-0.5;
          y = uniform_rng_(generator_)-0.5;
          s = x*x+y*y;
        } while(s > 0.25 || s == 0.);
        double coef = sqrt(-2*cppbugs::log_approx(s*4)/s);
        next_norm_ = coef*y;
        return coef*x;
      } else {
        double r = next_norm_;
        next_norm_ = NAN;
        return r;
      }
    }

    double uniform() { return uniform_rng_(generator_); }
  };

} // namespace cppbugs
#endif // MCMC_RNG_HPP
