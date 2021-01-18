/*
 *  QuadraticPolyOU.h
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
#ifndef QuadraticPoly_OU_H_
#define QuadraticPoly_OU_H_

#include "QuadraticPoly.h"
#include <armadillo>
#include <sstream>
#include <iostream>

namespace PCMBaseCpp {

typedef SPLITT::OrderedTree<SPLITT::uint, LengthAndRegime> OUTreeType;

template<class TreeType, class DataType>
struct CondGaussianOU: public CondGaussianOmegaPhiV {
  
  TreeType const& ref_tree_;
  
  // a 0-threshold for abs(Lambda_i + Lambda_j), where Lambda_i and Lambda_j are
  //  eigenvalues of the parameter matrix H. This threshold-values is used as a condition to
  // take the limit time of the expression `(1-exp(-Lambda_ij*time))/Lambda_ij` as
  //   `(Lambda_i+Lambda_j) --> 0`.
  double threshold_Lambda_ij_ = 1e-8;
  double threshold_SV_ = 1e-6;
  bool transpose_Sigma_x = false;
  
  //
  // model parameters
  //
  
  // number of traits
  uint k_;
  
  // number of regimes;
  uint R_; 
  
  // Each slice or column of the following cubes or matrices correponds to one regime
  arma::mat X0;
  
  arma::cube H;
  arma::mat Theta;
  arma::cube Sigma;
  arma::cube Sigmae;
  
  arma::cx_cube P;
  arma::cx_cube P_1;
  arma::cx_cube P_1SigmaP_1_t;
  
  // k-vectors of eigenvalues for each regime
  arma::cx_mat lambda;
  
  // matrices of sums of pairs of eigenvalues lambda_i+lambda_j for each regime
  arma::cx_cube Lambda_ij;
  
  // identity matrix
  arma::mat I;
  
  CondGaussianOU(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = R;
    this->threshold_Lambda_ij_ = ref_data.threshold_Lambda_ij_;
    this->threshold_SV_ = ref_data.threshold_SV_;
    this->transpose_Sigma_x = ref_data.transpose_Sigma_x;
    InitInternal();
  }
  
  CondGaussianOU(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = ref_data.R_;
    this->threshold_Lambda_ij_ = ref_data.threshold_Lambda_ij_;
    this->threshold_SV_ = ref_data.threshold_SV_;
    this->transpose_Sigma_x = ref_data.transpose_Sigma_x;
    InitInternal();
  }
  
  void InitInternal() {
    using namespace arma;
    this->I = eye(k_, k_);
    this->P = cx_cube(k_, k_, R_);
    this->P_1 = cx_cube(k_, k_, R_);
    this->P_1SigmaP_1_t = cx_cube(k_, k_, R_);
    this->lambda = cx_mat(k_, R_);
    this->Lambda_ij = cx_cube(k_, k_, R_);
  }
  
  arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) {
    using namespace arma;

    // for (const auto& i: par)
    //   std::cout << i << ' ';
    // std::cout<<"; "<<offset<<";"<<R_<<"; "<<k_<<"\n";
    
    uint npar = R_*(3*k_*k_ + 2*k_);
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os<<"QuadraticPolyOU.h:CondOU.SetParameter:: The length of the parameter vector minus offset ("<<par.size() - offset<<
        ") should be at least of R*(3k^2+2k), where k="<<k_<<" is the number of traits and "<<
          " R="<<R_<<" is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    X0 = mat(&par[offset], k_, R_);
    H = cube(&par[offset + k_*R_], k_, k_, R_);
    Theta = mat(&par[offset + (k_ + k_*k_)*R_], k_, R_);
    Sigma = cube(&par[offset + (k_ + k_*k_ + k_)*R_], k_, k_, R_);
    Sigmae = cube(&par[offset + (k_ + k_*k_ + k_ + k_*k_)*R_], k_, k_, R_);
    
    // cout<<"X0:\n"<<X0<<"\n";
    // cout<<"H:\n"<<H<<"\n";
    // cout<<"Theta:\n"<<Theta<<"\n";
    // cout<<"Sigma:\n"<<Sigma<<"\n";
    // cout<<"Sigmae\n"<<Sigmae<<"\n";
    
    if(transpose_Sigma_x) {
      for(uword r = 0; r < R_; r++) {
        Sigma.slice(r) = Sigma.slice(r).t() * Sigma.slice(r);
        Sigmae.slice(r) = Sigmae.slice(r).t() * Sigmae.slice(r);  
      }
    } else {
      for(uword r = 0; r < R_; r++) {
        Sigma.slice(r) = Sigma.slice(r) * Sigma.slice(r).t();
        Sigmae.slice(r) = Sigmae.slice(r) * Sigmae.slice(r).t();  
      }
    }
    InitInternal();
      
    for(uword r = 0; r < R_; ++r) {
      using namespace std;
      
      DecomposeH(lambda, P, P_1, H, r, threshold_SV_);
      
      // notice that we use st() instead of t() for P.slice(ri) to avoid conjugate (Hermitian) transpose.
      P_1SigmaP_1_t.slice(r) = P_1.slice(r) * Sigma.slice(r) * P_1.slice(r).st();
      
      PairSums(Lambda_ij.slice(r), lambda.col(r));
    }
    
    return npar;
  }
 
 
  void CalculateOmegaPhiV(uint i, arma::uword ri, arma::mat& omega, arma::cube& Phi, arma::cube& V) {
    using namespace arma;
      
      double ti = this->ref_tree_.LengthOfBranch(i).length_;
      
      Phi.slice(i) = real(P.slice(ri) * diagmat(exp(-ti * lambda.col(ri))) * P_1.slice(ri));
      omega.col(i) = (I - Phi.slice(i)) * Theta.col(ri);
      
      cx_mat fLambda_ij(k_, k_);
      CDFExpDivLambda(fLambda_ij, Lambda_ij.slice(ri), ti, threshold_Lambda_ij_);
      // notice that we use st() instead of t() for P.slice(ri) to avoid conjugate (Hermitian) transpose.
      V.slice(i) = real(P.slice(ri) * (fLambda_ij % P_1SigmaP_1_t.slice(ri)) * P.slice(ri).st());
      
      if(i < this->ref_tree_.num_tips()) {
        V.slice(i) += Sigmae.slice(ri);
      }
  }
};


class OU: public QuadraticPoly<OUTreeType> {
public:
  typedef OUTreeType TreeType;
  typedef QuadraticPoly<TreeType> BaseType;
  typedef OU MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;

  CondGaussianOU<TreeType, DataType> cond_dist_;
  
  OU(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }

  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};

typedef TraversalTaskWrapper<OU> QuadraticPolyOU;
}

#endif // QuadraticPoly_OU_H_
