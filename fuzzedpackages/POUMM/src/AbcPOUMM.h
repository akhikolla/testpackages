/**
 *  AbcPOUMM.h
 *  POUMM
 *
 * Copyright 2015-2019 Venelin Mitov
 *
 * This file is part of the R-package POUMM: The Phylogenetic Ornstein-Uhlenbeck Mixed Model
 *
 * POUMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * POUMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with POUMM.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */
#ifndef ABC_POUMM_H_
#define ABC_POUMM_H_

#include "./SPLITT.h"
#include <iostream>

namespace SPLITT {
template<class NameType>
struct NumericTraitData {
  // use const references to avoid copying of long vectors
  std::vector<NameType> const& names_;
  vec const& z_;
  vec const& se_;
  NumericTraitData(
    std::vector<NameType> const& names,
    vec const& z, vec const& se): names_(names), z_(z), se_(se) {}
};

template<class Tree>
class AbcPOUMM: public TraversalSpecification<Tree> {

public:
  typedef AbcPOUMM<Tree> MyType;
  typedef TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  double alpha, theta, sigmae2, sigma2;
  vec z, se;
  vec a, b, c;
  uvec visitDone;

  AbcPOUMM(TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.z_.size() != this->ref_tree_.num_tips() ||
       input_data.se_.size() != this->ref_tree_.num_tips()) {
      std::ostringstream oss;
      oss<<"The vectors z and se must be the same length as the number of tips ("<<
        this->ref_tree_.num_tips()<<"), but were"<<input_data.z_.size()<<
          " and "<<input_data.se_.size()<<" respectively.";
      throw std::invalid_argument(oss.str());
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->z = At(input_data.z_, ordNodes);
      this->se = At(input_data.se_, ordNodes);
      this->a = vec(this->ref_tree_.num_nodes());
      this->b = vec(this->ref_tree_.num_nodes());
      this->c = vec(this->ref_tree_.num_nodes());
      this->visitDone = uvec(this->ref_tree_.num_nodes());
    }
  };

  StateType StateAtRoot() const {
    vec res(3);
    res[0] = a[this->ref_tree_.num_nodes() - 1];
    res[1] = b[this->ref_tree_.num_nodes() - 1];
    res[2] = c[this->ref_tree_.num_nodes() - 1];
    return res;
  };

  void SetParameter(ParameterType const& par) {
    if(par.size() != 4) {
      throw std::invalid_argument(
      "The par vector should be of length 4 with \
      elements corresponding to alpha, theta, sigma and sigmae.");
    }
    if(par[0] < 0 || par[2] < 0 || par[3] < 0) {
      throw std::logic_error("The parameters alpha, sigma and sigmae should be non-negative.");
    }
    this->alpha = par[0];
    this->theta = par[1];
    this->sigma2 = par[2]*par[2];
    this->sigmae2 = par[3]*par[3];
  }

  inline void InitNode(uint i) {
    if(i < this->ref_tree_.num_tips()) {
      double sum_se2_sigmae2 = sigmae2 + se[i]*se[i];
      double z1 = z[i] - theta;
      
      if(sum_se2_sigmae2 != 0) {
        a[i] = -0.5 / sum_se2_sigmae2;
        b[i] = z1 / sum_se2_sigmae2;
        c[i] = -0.5 * (M_LN_2PI  + z1 * b[i] + log(sum_se2_sigmae2));
        
        visitDone[i] = 0;
        
          // gutalphasigma2[i] = e2talpha[i] + ((-0.5 / sum_se2_sigmae2) * sigma2) / fe2talpha;
          // 
          // // integration over g1 including e1 = z1 - g1
          // c[i] = -0.5 * log(gutalphasigma2[i]) -
          //   0.25 * sigma2 * z1z1[i - eFirst] / (sum_se2_sigmae2[i]*sum_se2_sigmae2[i]) /
          //     (fe2talpha[i] - alpha + (-0.5 / sum_se2_sigmae2[i]) * sigma2) +
          //       talpha[i] + (-0.5 * (M_LN_2PI  + z1z1[i-eFirst] / sum_se2_sigmae2[i]) - log_se_total[i]);
          // b[i] = (etalpha[i] * (z1[i] / sum_se2_sigmae2[i])) / gutalphasigma2[i];
          // a[i] = (-0.5 / sum_se2_sigmae2[i]) / gutalphasigma2[i];  
        
      } else {
        double z1z1 = z1 * z1;
        double t = this->ref_tree_.LengthOfBranch(i);
        double talpha = t * alpha;
        double etalpha = exp(talpha);
        double e2talpha = etalpha * etalpha;
        double fe2talpha;
        if(alpha != 0) {
          fe2talpha = alpha / (1 - e2talpha);
        } else {
          fe2talpha = -0.5 / t;
        }
        
        // integration over g1 including e1 = 0
        a[i] = fe2talpha / sigma2;  
        b[i] = -2 * etalpha * z1 * a[i];
        c[i] = talpha + 0.5 * log(-fe2talpha) -
          M_LN_SQRT_PI - log(sqrt(sigma2)) + e2talpha * z1z1 * a[i];
        visitDone[i] = 1;
      }
    } else {
      a[i] = b[i] = c[i] = 0;
      visitDone[i] = 0;
    }
  }

  // x / (e2talpha + (x * sigma2) / fe2talpha) = (fe2talpha / sigma2)
  //   x = e2talpha/(fe2talpha / sigma2) + (x * sigma2)/fe2talpha * sigma2/fe2talpha
  //   x = e2talpha*sigma2/fe2talpha +
  //     
  inline void VisitNode(uint i) {
    if(!visitDone[i]) {
      double t = this->ref_tree_.LengthOfBranch(i);
      double talpha = t * alpha;
      double etalpha = exp(talpha);
      double e2talpha = etalpha * etalpha;
      double fe2talpha;
      if(alpha != 0) {
        fe2talpha = alpha / (1 - e2talpha);
      } else {
        fe2talpha = -0.5 / t;
      }
      double gutalphasigma2 = e2talpha + (a[i] * sigma2) / fe2talpha;
      
      c[i] = -0.5 * log(gutalphasigma2) - 0.25 * sigma2 * b[i] * b[i] /
        (fe2talpha - alpha + a[i] * sigma2) + talpha + c[i];
      b[i] = (etalpha * b[i]) / gutalphasigma2;
      a[i] /= gutalphasigma2;
    }
    
  }

  inline void PruneNode(uint i, uint i_parent) {
    a[i_parent] += a[i];
    b[i_parent] += b[i];
    c[i_parent] += c[i];
  }

};

typedef TraversalTask<
  AbcPOUMM<OrderedTree<uint, double>> > ParallelPruningAbcPOUMM;
}
#endif //ABC_POUMM_H_
