/*
 *  QuadraticPoly1D.h
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
#ifndef PCMBase_QuadraticPoly1D_H_
#define PCMBase_QuadraticPoly1D_H_

#include "SPLITT.h"
#include "QuadraticPolyCommon.h"
#include <armadillo>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>
#include <iostream>


namespace PCMBaseCpp {

template<class NameType>
struct NumericTraitData1D {
  // use const references to avoid copying of big data objects
  
  // tip-names correponding to the columns in Pc_
  std::vector<NameType> const& names_;
  
  arma::mat const& X_;
  arma::cube const& VE_;
  
  uint R_;
  
  std::vector<std::string> regime_models_;
  double threshold_SV_;
  double threshold_EV_;
  double threshold_skip_singular_;
  double threshold_Lambda_ij_;
  bool skip_singular_;
  bool transpose_Sigma_x;
  double NA_double_;
  
  NumericTraitData1D(
    std::vector<NameType> const& names,
    arma::mat const& X,
    arma::cube const& VE,
    uint R,
    std::vector<std::string> regime_models,
    double threshold_SV,
    double threshold_EV,
    double threshold_skip_singular,
    bool skip_singular,
    bool transpose_Sigma_x,
    double threshold_Lambda_ij, 
    double NA_double_): names_(names), X_(X), VE_(VE), 
    R_(R), regime_models_(regime_models),
    threshold_SV_(threshold_SV), 
    threshold_EV_(threshold_EV), 
    threshold_skip_singular_(threshold_skip_singular),
    threshold_Lambda_ij_(threshold_Lambda_ij),
    skip_singular_(skip_singular),
    transpose_Sigma_x(transpose_Sigma_x),
    NA_double_(NA_double_) {}
};

template<class MatType>
inline void CDFExpDivLambda1D(MatType& fLambda_ij, MatType const& Lambda_ij, double time, double threshold_Lambda_ij) {
  using namespace arma;
  if(fabs(Lambda_ij) < threshold_Lambda_ij) {
    fLambda_ij = time;
  } else {
    fLambda_ij = (1.0 - exp(-Lambda_ij * time)) / Lambda_ij;
  }
}

// Conditional Gaussian distribution of trait vector at a daughter 
// node, Xi, given trait vector at its parent, Xj, assuming that the conditional mean
// depends linearly in Xj, i.e. Mean(Xi) = omega + Phi %*% Xj, and the variance does
// not depend on Xj. This is an abstract class implemented by specific models
class CondGaussianOmegaPhiV1D {
public:
  virtual arma::uword SetParameter(
      std::vector<double> const& par, 
      arma::uword offset) = 0;
  
  virtual void CalculateOmegaPhiV(
      uint i, 
      arma::uword ri, 
      arma::vec& omega, 
      arma::vec& Phi, 
      arma::vec& V) = 0;
  
  virtual ~CondGaussianOmegaPhiV1D() {}
};

template<class Tree>
class QuadraticPoly1D: public SPLITT::TraversalSpecification<Tree> {
public:
  typedef SPLITT::TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef std::vector<double> StateType;
  
  // singular value threshold for the determinant of V_i
  double threshold_SV_ = 1e-6;
  // positive eigenvalue threshold for V_i. 
  double threshold_EV_ = 1e-4;
  
  // threshold specifying the maximum allowed branch length for skipping the branch
  // if the corresponding matrix V is singular. This option matters only if
  // skip_singular_ is true.
  double threshold_skip_singular_ = 1e-4;
  
  // denotes branches for which V is singular
  // the treatment of these branches depends on the option PCMBase.Singular.Skip
  // If this option is set to TRUE (default), then these branches are treated as
  // 0-length and the L,m,r values accumulated from their children are added 
  // up to their parent branches without modification.
  // ATTENTION: using std::vector<bool> instead of std::vector<int> is causing a 
  // bug in parallel mode (omp for simd). This was found in the unit test for 
  // the White model, test-White.R in PCMBase, using the intel compiler icpc 
  // (icpc version 17.0.5 (gcc version 4.9.0 compatibility)), command line:
  // icpc -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/usr/local/include -I/usr/local/include/freetype2 -I/opt/X11/include -I"/Users/vmitov/Library/R/3.3/library/Rcpp/include" -I"/Users/vmitov/Library/R/3.3/library/RcppArmadillo/include"   -fPIC  -Wall -mtune=core2 -g -O2  -std=c++11 -fopenmp -Wall -O2 -march=native -c Rcpp.cpp -o Rcpp.o
  // the problem appears to be solved when using std::vector<int>.
  std::vector<int> singular_branch_;
  
  bool skip_singular_;
  
  
  //
  // Input data consists of multiple trait values for each tip. Each
  // column corresponds to a tip.
  //
  arma::vec X;
  
  // 
  // Standard error variance for each entry in the input data;
  //
  arma::vec VE;
  
  //
  // Coefficients used to calculate L, m, r, which are calculated in the
  // InitNode method of derived classes.
  //
  arma::vec A;
  arma::vec b;
  arma::vec C;
  arma::vec d;
  arma::vec E;
  arma::vec f;
  
  //
  // linear intercept and coefficient for the conditional mean: 
  // E[X_i|X_j] = omega + Phi X_j
  arma::vec omega;
  arma::vec Phi;
  
  // 
  // Variance - covariance matrices
  //
  arma::vec V;
  arma::vec V_1;
  
  
  //
  // Coefficients of the quadratic polynomial for each node
  // X.col(i).t() * L.slice(i) * X.col(i) + m.col(i) * X.col(i) + r(i)
  //
  arma::vec L;
  arma::vec m;
  arma::vec r;
  
  std::vector<CondGaussianOmegaPhiV1D*> ptr_cond_dist_;
  
  QuadraticPoly1D(
    TreeType const& tree,
    NumericTraitData1D<typename TreeType::NodeType> const& input_data):
    
    BaseType(tree),
    
    threshold_SV_(input_data.threshold_SV_),
    threshold_EV_(input_data.threshold_EV_),
    threshold_skip_singular_(input_data.threshold_skip_singular_),
    singular_branch_(tree.num_nodes(), false),
    skip_singular_(input_data.skip_singular_),
    
    X(input_data.X_.row(0).t()),
    VE(input_data.VE_.tube(0, 0)),
    
    // all these fields have to be initialized with 0 during SetParameter.
    A(tree.num_nodes()),
    b(tree.num_nodes()),
    C(tree.num_nodes()),
    d(tree.num_nodes()),
    E(tree.num_nodes()),
    f(tree.num_nodes()),
    
    omega(tree.num_nodes()),
    Phi(tree.num_nodes()),
    
    V(tree.num_nodes()),
    V_1(tree.num_nodes()),
    
    L(tree.num_nodes()),
    m(tree.num_nodes()),
    r(tree.num_nodes()) {
    
    A.fill(input_data.NA_double_);
    b.fill(input_data.NA_double_);
    C.fill(input_data.NA_double_);
    d.fill(input_data.NA_double_);
    E.fill(input_data.NA_double_);
    f.fill(input_data.NA_double_);
    omega.fill(input_data.NA_double_);
    Phi.fill(input_data.NA_double_);
    V.fill(input_data.NA_double_);
    V_1.fill(input_data.NA_double_);
    L.fill(input_data.NA_double_);
    m.fill(input_data.NA_double_);
    r.fill(input_data.NA_double_);
    
    arma::uvec ordTips(
        this->ref_tree_.OrderNodesPosType(
            input_data.names_, static_cast<arma::uword>(SPLITT::G_NA_UINT)));
    
    this->X(arma::span(0, this->ref_tree_.num_tips() - 1)) = X(ordTips);
    this->VE(arma::span(0, this->ref_tree_.num_tips() - 1)) = VE(ordTips);
    
    SPLITT::uvec node_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), tree.num_nodes());
    arma::uvec ordNodes(
        this->ref_tree_.OrderNodesPosType(
            node_names, static_cast<arma::uword>(SPLITT::G_NA_UINT)));
  }
  
  StateType StateAtNode(arma::uword i) const {
    using namespace std;
    
    StateType res(13);
    
    arma::uword j = 0;
    res[j++] = L(i); 
    res[j++] = m(i);
    res[j++] = r(i);
    res[j++] = A(i);
    res[j++] = b(i);
    res[j++] = C(i);
    res[j++] = d(i);
    res[j++] = E(i);
    res[j++] = f(i);
    res[j++] = omega(i);
    res[j++] = Phi(i);
    res[j++] = V(i);
    res[j++] = V_1(i);
  
    return res;
  }
  
  StateType StateAtRoot() const {
    return StateAtNode(this->ref_tree_.num_nodes() - 1);
  }
  
  // this function is to be called by daughter classes only after they have
  // initialized omega, Phi and V
  inline void CalculateAbCdEf(uint i) {
    using namespace arma;
    
    A(i) = -0.5 * V_1(i);
    E(i) = Phi(i) * V_1(i);
    b(i) = V_1(i) * omega(i);
    C(i) = -0.5 * E(i) * Phi(i);
    d(i) = -E(i) * omega(i);
    f(i) = -0.5*(M_LN_2PI + log(det(V(i))) + omega(i) * V_1(i) * omega(i) );
  }
  
  inline void InitLmr(uint i) {
    L(i) = 0.0;
    m(i) = 0.0;
    r(i) = 0.0;
    singular_branch_[i] = 0;
  }
  
  inline void InitNode(uint i) {
    using namespace arma;
    using namespace std;
    InitLmr(i);
    
    if(i < this->ref_tree_.num_nodes() - 1) {
      
      auto ri = this->ref_tree_.LengthOfBranch(i).regime_;
      auto ti = this->ref_tree_.LengthOfBranch(i).length_;
      
      if(ptr_cond_dist_.size() == 1) {
        ptr_cond_dist_[0]->CalculateOmegaPhiV(i, ri, omega, Phi, V);
      } else {
        ptr_cond_dist_[ri]->CalculateOmegaPhiV(i, 0, omega, Phi, V);
      }
      
      // handle measurement error
      if(i < this->ref_tree_.num_tips()) {
        // tip node
        V(i) += VE(i);
      }
      
      if( V(i) < threshold_SV_ ) {
        singular_branch_[i] = 1;
        if(!skip_singular_ || ti > threshold_skip_singular_) {
          ostringstream oss;
          oss<<"QuadraticPoly1D.h:InitNode:: V for node "<<
            this->ref_tree_.FindNodeWithId(i)<<" is smaller than threshold_SV_:"<<
              V(i)<<"<"<<threshold_SV_<<
              ". Check the model parameters, the length of the branch leading"<<
                "to the node, and the PCMBase.Threshold.SV option"<<
                  " For details on this error, read the User Guide.";
          this->SetError(oss.str()); // throw logic_error(oss.str());  
        } 
      } 
      
      if(!singular_branch_[i]) {
        // Check V is positive definite: all eigen-values must be strictly positive
        if( V(i) < threshold_EV_ ) {
          ostringstream oss;
          oss<<"QuadraticPoly1D.h:InitNode:: V for node "<<
            this->ref_tree_.FindNodeWithId(i)<<
              " is nearly 0 or negative: "<<V(i)<<"<"<<threshold_EV_<<
                ". Check the model parameters and the PCMBase.Threshold.EV option.";
          this->SetError(oss.str()); //  throw logic_error(oss.str());
        }
      
        V_1(i) = 1.0 / V(i);
        CalculateAbCdEf(i);  
      }
    }  
    
  }
  
  inline void VisitNode(uint i) {
    using namespace arma;
    using namespace std;
    
    if(!singular_branch_[i]) {
      
      if(i < this->ref_tree_.num_tips()) {
        // ensure symmetry of L.slice(i)
        L(i) = C(i);
        r(i) = X(i) * A(i) * X(i) + X(i) * b(i) + f(i);
        m(i) = d(i) + E(i) * X(i);
      } else {
        double AplusL = A(i) + L(i);
        
        double AplusL_1 = 1.0 / AplusL;
        double EAplusL_1 = E(i) * AplusL_1;
        double logDetVNode = log(-2.0*AplusL);
        
        r(i) = f(i) + r(i) + 0.5 * M_LN_2PI - 0.5 * logDetVNode -
          0.25 * (b(i) + m(i)) * AplusL_1 * (b(i) + m(i));
        m(i) = d(i) - 0.5*EAplusL_1 * (b(i) + m(i));
        L(i) = C(i) - 0.25 * EAplusL_1 * E(i);
      }
    }
  }
  
  inline void PruneNode(uint i, uint i_parent) {
    L(i_parent) += L(i);
    m(i_parent) += m(i);
    r(i_parent) += r(i);
  }
};
}

#endif // PCMBase_QuadraticPoly1D_H_
