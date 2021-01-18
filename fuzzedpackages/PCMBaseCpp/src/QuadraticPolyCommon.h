/*
 *  QuadraticPolyCommon.h
 *  PCMBaseCpp
 *
 * Copyright 2017,2018 Venelin Mitov
 *
 * This file is part of PCMBaseCpp: A C++ backend for calculating the likelihood of phylogenetic comparative models.
 *
 * PCMBaseCpp is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU General Public License as
 * published by the Free Software Foundation.
 *
 * PCMBaseCpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with PCMBaseCpp.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */
#ifndef PCMBase_QuadraticPolyCommon_H_
#define PCMBase_QuadraticPolyCommon_H_

#include "SPLITT.h"
#include <armadillo>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <exception>


namespace PCMBaseCpp {

typedef unsigned int uint;

struct LengthAndRegime {
  double length_;
  arma::uword regime_;

  LengthAndRegime() {}

  LengthAndRegime(double length, arma::uword regime):
    length_(length), regime_(regime) {}
};

struct LengthRegimeAndJump {
  double length_;
  arma::uword regime_;
  arma::u8 jump_;
  
  LengthRegimeAndJump() {}
  
  LengthRegimeAndJump(double length, arma::uword regime, arma::u8 jump):
    length_(length), regime_(regime), jump_(jump) {}
};

template<class RegimeType>
std::vector<arma::uword> mapRegimesToIndices(
    std::vector<RegimeType> const& regimes,
    std::vector<RegimeType> const& regimes_unique) {
  
  if(regimes_unique.size() == 0) {
    throw std::logic_error("QuadraticPolyCommon.h:mapRegimesToIndices:: regimes_unique has 0 length but should have at least one regime.");
  }
  std::unordered_map<RegimeType, arma::uword> map_regimes;
  arma::uword next_regime = 0;
  for(auto r: regimes_unique) {
    auto it = map_regimes.insert(std::pair<RegimeType, arma::uword>(r, next_regime));
    if(!it.second) {
      std::ostringstream os;
      os<<"QuadraticPolyCommon.h:mapRegimesToIndices:: The regime named '"<<r<<"' is dupliclated. Remove duplicates from regimes_unique.";
      throw std::logic_error(os.str());
    } else {
      ++next_regime;
    }
  }
  std::vector<arma::uword> regimeIndices;
  for(auto r: regimes) {
    auto it = map_regimes.find(r);
    if(it == map_regimes.end()) {
      std::ostringstream os;
      os<<"QuadraticPolyCommon.h:mapRegimesToIndices:: The regime named '"<<r<<"' was not found in regimes_unique.";
      throw std::logic_error(os.str());
    } else {
      regimeIndices.push_back(it->second);
    }
  }
  return regimeIndices;
}

// Wraps the call TraverseTree with an exception handler. 
// Needed for https://github.com/venelin/PCMBaseCpp/issues/1#issue-507813590.
template<class TraversalSpec> 
class TraversalTaskWrapper {
public:
  typedef typename SPLITT::TraversalTask<TraversalSpec> TraversalTaskType;
  
  typedef typename TraversalTaskType::StateType StateType;
  typedef typename TraversalTaskType::ParameterType ParameterType;
  typedef typename TraversalTaskType::NodeType NodeType;
  typedef typename TraversalTaskType::LengthType LengthType;
  typedef typename TraversalTaskType::DataType DataType;
  
  typedef typename TraversalTaskType::TreeType TreeType;
  typedef typename TraversalTaskType::AlgorithmType AlgorithmType;
  
  TraversalTaskType taskObject_;

  TraversalTaskWrapper(
    std::vector<NodeType> const& branch_start_nodes,
    std::vector<NodeType> const& branch_end_nodes,
    std::vector<LengthType> const& branch_lengths,
    DataType const& data): 
  taskObject_(branch_start_nodes, branch_end_nodes, branch_lengths, data) {}
  
  std::string TraverseTree(ParameterType const& par, uint mode) {
    try {
      taskObject_.spec().ResetError();
      taskObject_.TraverseTree(par, mode);
    } catch(std::logic_error& e) {
      return std::string("logic_error: ") + e.what();
    } catch(std::exception& e) {
      return std::string("exception: ") + e.what();
    } catch(...) {
      return std::string("unknown error.");
    }
    return taskObject_.spec().GetError();
  }
  
  StateType StateAtNode(uint i) {
    return taskObject_.StateAtNode(i);
  }
  
  TreeType & tree() {
    return taskObject_.tree();
  }
  AlgorithmType & algorithm() {
    return taskObject_.algorithm();
  }
};
}

#endif // PCMBase_QuadraticPolyCommon_H_
