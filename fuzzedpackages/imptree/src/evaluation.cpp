// evaluation.cpp: Imprecise Classification Trees
//
// Copyright (C) 2018  Paul Fink
//
// This file is part of imptree.
//
// imptree is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// imptree is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with imptree.  If not, see <https://www.gnu.org/licenses/>.

#include "evaluation.h"

Evaluation::Evaluation(const double utility, const Dominance dominance,
                       const std::vector<ProbInterval> &probInts, Data observations) : 
  utility_(utility), dominance_(dominance), probInts_(probInts), observations_(observations)
{
  boolPrediction_ = Rcpp::LogicalMatrix(probInts_.at(0).freq.size(), probInts_.size());
  evaluate();
}

void Evaluation::evaluate() {
  
  int nObs = probInts_.size();
  for (int idx = 0; idx < nObs; ++idx) {
    updateCredalStatistics(idx);
  }
  finalizeCredalStatistics();
}

Rcpp::List Evaluation::summary() const {
  return Rcpp::List::create(
    Rcpp::Named("nObs") = probInts_.size(),
    Rcpp::Named("deter") = obs_det_,
    Rcpp::Named("nObsIndet") = obs_indet_,
    Rcpp::Named("indetSize") = size_indet_,
    Rcpp::Named("acc_single") = acc_single_,
    Rcpp::Named("acc_set") = acc_set_,
    Rcpp::Named("acc_disc") = acc_disc_,
    Rcpp::Named("acc_util") = acc_util_
  );
}

Rcpp::List Evaluation::probIntervalList() const {
  
  std::vector<Rcpp::NumericMatrix> res;
  for(ProbInterval probint : probInts_) {
    res.push_back(probint.toMatrix());
  }
  return Rcpp::wrap(res);
}

/**
 * Quadratic utility u(acc_disc)=(1+4*util*(1-acc_disc))*acc_disc
 *
 * @param acc_disc raw discounted accuracy
 * @return u(x0.5)
 */
double Evaluation::quadratic_utility(double acc_disc) {
  return ((1.0 + 4.0 * utility_ * (1.0 - acc_disc)) * acc_disc);
}

// true means it is not dominated
std::vector<bool> Evaluation::computeNonDominatedSet(const ProbInterval &probint) {
  
  int n = probint.freq.size();
  std::vector<bool> set(n, false);
  if (dominance_ == Dominance::maximality) {
    
    int idx = std::distance(probint.upper.begin(), std::max_element(probint.upper.begin(), probint.upper.end()));
    set[idx] = true;
    
  } else if (dominance_ == Dominance::interval) {
    
    for (int k = 0; k < n; ++k) {
      for (int i = 0; i < n; ++i) {
        if(i == k || set[i]) {
          continue;
        }
        if(probint.upper[k] < probint.lower[i]) {
          set[k] = true;
        }
      }
    }
    std::transform (set.begin(), set.end(), set.begin(), std::logical_not<bool>());
  }
  
  return set;
}

void Evaluation::updateCredalStatistics(int obsIdx) {
  
  Rcpp::LogicalVector ndset = Rcpp::wrap(computeNonDominatedSet(probInts_.at(obsIdx)));
  bool inset = ndset[observations_.data(obsIdx, observations_.classidx)];
  
  int ndset_size = sum(ndset);
  // determinant prediction
  if (ndset_size == 1) {
    obs_det_++;
    if (inset) {
      acc_single_++;
      acc_disc_++;
      acc_util_++;
    }
  } else {
    obs_indet_++;
    size_indet_ += ndset_size;
    if (inset) {
      acc_set_++;
      acc_disc_ += 1.0 / ndset_size;
      acc_util_ += quadratic_utility(1.0 / ndset_size);
    }
  }
  boolPrediction_(Rcpp::_,obsIdx) = ndset;
}

void Evaluation::finalizeCredalStatistics() {
  
  int nObs = observations_.data.nrow();
  
  if(obs_det_ > 0) {
    acc_single_ /= obs_det_;
  } else {
    acc_single_ = NA_REAL;
  }
  if (obs_indet_ > 0) {
    acc_set_ /= obs_indet_;
    size_indet_ /= obs_indet_;
  } else {
    acc_set_= NA_REAL;
    size_indet_ = NA_REAL;
  }
  acc_disc_ /= nObs;
  acc_util_ /= nObs;
  obs_det_ /=  nObs;
  
}
