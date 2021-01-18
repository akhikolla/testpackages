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

#ifndef MCMC_DYNAMIC_HPP
#define MCMC_DYNAMIC_HPP

#include <vector>
#include <cppbugs/mcmc.math.hpp>
#include <cppbugs/mcmc.object.hpp>

namespace cppbugs {

  template<typename T>
  class Dynamic;

  template<typename T>
  class Dynamic<T&> : public MCMCObject {
    bool save_history_;

  public:
    std::vector<T> history;
    T& value;
    T old_value;

    Dynamic(T& shape): MCMCObject(), save_history_(true), value(shape), old_value(shape) {}

    static void fill(arma::ivec& x) { x.fill(0); }
    static void fill(arma::mat& x) { x.fill(0); }
    static void fill(double& x) { x = 0; }
    static void fill(int& x) { x = 0; }

    T mean() const {
      if(history.size() == 0) {
        return T();
      }

      T ans(*history.begin());
      fill(ans);
      for(T v : history)
        ans += v;
      ans /= static_cast<double>(history.size());
      return ans;
    }

    void setSaveHistory(const bool save_history) {
      save_history_ = save_history;
    }

    void preserve() { old_value = value; }
    void revert() { value = old_value; }
    void tally() { if(save_history_) { history.push_back(value); } }
    double size() const { return dim_size(value); }
  };

} // namespace cppbugs
#endif //MCMC_DYNAMIC_HPP
