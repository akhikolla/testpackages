/*
 *  QuadraticPolyOU1D.h
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
#ifndef QuadraticPoly_OU1D_H_
#define QuadraticPoly_OU1D_H_

#include "QuadraticPoly1D.h"
#include <armadillo>
#include <sstream>

namespace PCMBaseCpp {

typedef SPLITT::OrderedTree<SPLITT::uint, LengthAndRegime> OU1DTreeType;

template<class TreeType, class DataType>
struct CondGaussianOU1D: public CondGaussianOmegaPhiV1D {
  
  TreeType const& ref_tree_;
  
  // a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
  //  eigenvalues of the parameter matrix H. This threshold-values is used as a condition to
  // take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
  //   `(Lambda_i+Lambda_j) --> 0`.
  double threshold_Lambda_ij_ = 1e-8;
  double threshold_SV_ = 1e-6;
  
  //
  // model parameters
  //
  
  // number of regimes;
  uint R_; 
  
  // Each element of the following vectors correponds to one regime
  arma::vec X0;
  
  arma::vec H;
  arma::vec Theta;
  arma::vec Sigma;
  arma::vec Sigmae;
  
  CondGaussianOU1D(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->R_ = R;
    this->threshold_Lambda_ij_ = ref_data.threshold_Lambda_ij_;
    this->threshold_SV_ = ref_data.threshold_SV_;
    InitInternal();
  }
  
  CondGaussianOU1D(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->R_ = ref_data.R_;
    this->threshold_Lambda_ij_ = ref_data.threshold_Lambda_ij_;
    this->threshold_SV_ = ref_data.threshold_SV_;
    InitInternal();
  }
  
  void InitInternal() {
  }
  
  arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) {
    using namespace arma;

    uint npar = R_ * 5;
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os<<"QuadraticPolyOU1D.h:CondOU1D.SetParameter:: The length of the parameter vector minus offset ("<<par.size() - offset<<
        ") should be at least of R*5, where R="<<R_<<" is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    X0 = vec(&par[offset], R_);
    H = vec(&par[offset + R_], R_);
    Theta = vec(&par[offset + 2*R_], R_);
    Sigma = vec(&par[offset + 3*R_], R_);
    Sigmae = vec(&par[offset + 4*R_], R_);
    
    for(uword r = 0; r < R_; r++) {
      Sigma(r) = Sigma(r) * Sigma(r);
      Sigmae(r) = Sigmae(r) * Sigmae(r);  
    }
    
    InitInternal();
  
    return npar;
  }
 
 
  void CalculateOmegaPhiV(uint i, arma::uword ri, arma::vec& omega, arma::vec& Phi, arma::vec& V) {
    using namespace arma;
      
      double ti = this->ref_tree_.LengthOfBranch(i).length_;
      
      Phi(i) = exp(-ti * H(ri));
      omega(i) = (1.0 - Phi(i)) * Theta(ri);
      
      double fLambda_ij;
      CDFExpDivLambda1D(fLambda_ij, 2 * H(ri), ti, threshold_Lambda_ij_);
      V(i) = fLambda_ij * Sigma(ri);
      
      if(i < this->ref_tree_.num_tips()) {
        V(i) += Sigmae(ri);
      }
  }
};

class OU1D: public QuadraticPoly1D<OU1DTreeType> {
public:
  typedef OU1DTreeType TreeType;
  typedef QuadraticPoly1D<TreeType> BaseType;
  typedef OU1D MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData1D<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;

  CondGaussianOU1D<TreeType, DataType> cond_dist_;
  
  OU1D(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }

  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};

typedef TraversalTaskWrapper<OU1D> QuadraticPolyOU1D;
}

#endif // QuadraticPoly_OU1D_H_
