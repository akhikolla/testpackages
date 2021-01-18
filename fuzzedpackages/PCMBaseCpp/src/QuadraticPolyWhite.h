/*
 *  QuadraticPolyWhite.h
 *  PCMBaseCpp
 *
 * Copyright 2017,2018 Venelin Mitov
 *
 * This file is part of PCMBaseCpp: A C++ backend for calculating the likelihood
 *  of phylogenetic comparative models.
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
#ifndef QuadraticPoly_WHITE_H_
#define QuadraticPoly_WHITE_H_

#include "QuadraticPoly.h"
#include <armadillo>
#include <sstream>

namespace PCMBaseCpp {

typedef SPLITT::OrderedTree<SPLITT::uint, LengthAndRegime> WhiteTreeType;

template<class TreeType, class DataType>
struct CondGaussianWhite: public CondGaussianOmegaPhiV {
  
  TreeType const& ref_tree_;
  
  //
  // model parameters
  //
  
  // number of traits
  uint k_;
  
  // number of regimes;
  uint R_; 
  
  // Each slice or column of the following cubes or matrices correponds to one regime
  arma::mat X0;
  
  arma::cube Sigmae;
  
  // matrices of sums of pairs of eigenvalues lambda_i+lambda_j for each regime
  
  arma::mat I;
  
  CondGaussianWhite(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = R;
    this->I = arma::eye(k_, k_);
  }
  
  CondGaussianWhite(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = ref_data.R_;
    this->I = arma::eye(k_, k_);
  }
  
  arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) {
    using namespace arma;
    
    uint npar = R_*(k_*k_ + k_);
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os<<"QuadraticPolyWhite.h:CondGaussianWhite.SetParameter:: The length of the parameter vector minus offset ("<<par.size() - offset<<
        ") should be at least of R*(k^2+k), where k="<<k_<<" is the number of traits and "<<
          " R="<<R_<<" is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    X0 = mat(&par[offset], k_, R_);
    Sigmae = cube(&par[offset + k_*R_], k_, k_, R_);
    
    for(uword r = 0; r < R_; r++) {
      Sigmae.slice(r) = Sigmae.slice(r) * Sigmae.slice(r).t();  
    }
    
    return npar;
  }
  
  void CalculateOmegaPhiV(uint i, arma::uword ri, arma::mat& omega, arma::cube& Phi, arma::cube& V) {
    using namespace arma;
    
    omega.col(i).fill(0);
    Phi.slice(i) = I;
    
    V.slice(i).fill(0);
    
    if(i < this->ref_tree_.num_tips()) {
      V.slice(i) += Sigmae.slice(ri);
    }
  }
};


class White: public QuadraticPoly<WhiteTreeType> {
public:
  typedef WhiteTreeType TreeType;
  typedef QuadraticPoly<TreeType> BaseType;
  typedef White MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;

  CondGaussianWhite<TreeType, DataType> cond_dist_;
  
  White(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }
  
  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};


typedef TraversalTaskWrapper<White> QuadraticPolyWhite;
}

#endif // QuadraticPoly_WHITE_H_
