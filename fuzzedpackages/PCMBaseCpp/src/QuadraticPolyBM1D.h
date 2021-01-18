/*
 *  QuadraticPolyBM1D.h
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
#ifndef QuadraticPoly_BM1D_H_
#define QuadraticPoly_BM1D_H_

#include "QuadraticPoly1D.h"
#include <armadillo>
#include <sstream>
#include <iostream>

namespace PCMBaseCpp {

typedef SPLITT::OrderedTree<SPLITT::uint, LengthAndRegime> BM1DTreeType;

template<class TreeType, class DataType>
struct CondGaussianBM1D: public CondGaussianOmegaPhiV1D {
  
  TreeType const& ref_tree_;
  
  //
  // model parameters
  //
  
  // number of regimes;
  uint R_; 
  
  // Each element correponds to one regime
  arma::vec X0;
  
  // Each element of the following vectors correponds to one regime
  arma::vec Sigma;
  arma::vec Sigmae;
  
  CondGaussianBM1D(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->R_ = R;
  }
  
  CondGaussianBM1D(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->R_ = ref_data.R_;
  }
  
  arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) {
    using namespace arma;
    
    uint npar = R_ * 3;
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os<<"QuadraticPolyBM1D.h:CondBM1D.SetParameter:: The length of the parameter vector minus offset ("<<par.size() - offset<<
        ") should be at least of R*3, where R="<<R_<<" is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    X0 = vec(&par[offset], R_);
    Sigma = vec(&par[offset + R_], R_);
    Sigmae = vec(&par[offset + 2*R_], R_);
    
    // By convention the parameters Sigma and Sigmae are square roots.
    for(uword r = 0; r < R_; r++) {
      Sigma(r) = Sigma(r) * Sigma(r);
      Sigmae(r) = Sigmae(r) * Sigmae(r);  
    }
    
    return npar;
  }
  
  void CalculateOmegaPhiV(uint i, arma::uword ri, arma::vec& omega, arma::vec& Phi, arma::vec& V) {
    using namespace arma;
    
    double ti = this->ref_tree_.LengthOfBranch(i).length_;
    
    omega(i) = 0.0;
    Phi(i) = 1.0;
    
    V(i) = ti * Sigma(ri); 
    
    if(i < this->ref_tree_.num_tips()) {
      V(i) += Sigmae(ri);
    }
  }
};


class BM1D: public QuadraticPoly1D<BM1DTreeType> {
public:
  typedef BM1DTreeType TreeType;
  typedef QuadraticPoly1D<TreeType> BaseType;
  typedef BM1D MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData1D<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;
  
  CondGaussianBM1D<TreeType, DataType> cond_dist_;
  
  BM1D(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }
  
  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};


typedef TraversalTaskWrapper<BM1D> QuadraticPolyBM1D;
}

#endif // QuadraticPoly_BM1D_H_
