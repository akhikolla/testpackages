/*
 *  QuadraticPolyMixedGaussian.h
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
#ifndef QuadraticPoly_MixedGaussian_H_
#define QuadraticPoly_MixedGaussian_H_

#include "QuadraticPoly.h"
#include "QuadraticPolyBM.h"
#include "QuadraticPolyJOU.h"
#include "QuadraticPolyOU.h"
#include "QuadraticPolyDOU.h"
#include <armadillo>
#include <sstream>

namespace PCMBaseCpp {


typedef SPLITT::OrderedTree<SPLITT::uint, LengthRegimeAndJump> MixedGaussianTreeType;

class MixedGaussian: public QuadraticPoly<MixedGaussianTreeType> {
public:
  typedef MixedGaussianTreeType TreeType;
  typedef QuadraticPoly<TreeType> BaseType;
  typedef MixedGaussian MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;
  
  MixedGaussian(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data) {
    
    for(uint i = 0; i < input_data.regime_models_.size(); i++) {
      std::string modelName = input_data.regime_models_[i];
      if(modelName == "BM") {
        BaseType::ptr_cond_dist_.push_back(new CondGaussianBM<TreeType, DataType>(tree, input_data, 1));  
      } else if(modelName == "OU") {
        BaseType::ptr_cond_dist_.push_back(new CondGaussianOU<TreeType, DataType>(tree, input_data, 1));  
      } else if(modelName == "JOU") {
        BaseType::ptr_cond_dist_.push_back(new CondGaussianJOU<TreeType, DataType>(tree, input_data, 1));  
      } else if(modelName == "DOU") {
        BaseType::ptr_cond_dist_.push_back(new CondGaussianDOU<TreeType, DataType>(tree, input_data, 1));  
      } 
    }
    
  }
  
  void SetParameter(ParameterType const& par) {
    uint offset = 0;
    for(uint i = 0; i < ptr_cond_dist_.size(); i++) {
      offset += ptr_cond_dist_[i]->SetParameter(par, offset);
    }
  }
  
  ~MixedGaussian() {
    for(uint i = 0; i < BaseType::ptr_cond_dist_.size(); i++) {
      delete BaseType::ptr_cond_dist_[i];
    }
  }
};


typedef TraversalTaskWrapper<MixedGaussian> QuadraticPolyMixedGaussian;
}

#endif // QuadraticPoly_MixedGaussian_H_
