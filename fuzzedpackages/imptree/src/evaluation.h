// evaluation.h: Imprecise Classification Trees
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

#ifndef RCPP_IMPTREE_EVALUATION_H
#define RCPP_IMPTREE_EVALUATION_H

#include <Rcpp.h>
#include "enums.h"
#include "structs.h"

class Evaluation {
  
  double utility_;
  Dominance dominance_;
  std::vector<ProbInterval> probInts_;
  Rcpp::LogicalMatrix boolPrediction_;
  Data observations_;
  
  double obs_det_ = 0.0;
  int obs_indet_ = 0;
  
  double size_indet_ = 0.0;
  
  double acc_disc_ = 0.0;
  double acc_util_ = 0.0;
  
  double acc_single_ = 0.0;
  double acc_set_ = 0.0;
  
  double quadratic_utility(double acc_disc);
  std::vector<bool> computeNonDominatedSet(const ProbInterval &probint);
  void updateCredalStatistics(int obsIdx);
  void finalizeCredalStatistics();
  void evaluate();
  
public:
  Evaluation(const double utility, const Dominance dominance,
             const std::vector<ProbInterval> &probInts, Data observations);
  
  Rcpp::List summary() const;
  Rcpp::List probIntervalList() const;
  Rcpp::LogicalMatrix predictions() const {return boolPrediction_;} 
  
};

#endif /*RCPP_IMPTREE_EVALUATION_H*/
