// node.h: Imprecise Classification Trees
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

// [[Rcpp::plugins(cpp11)]]
#ifndef RCPP_IMPTREE_NODE_H
#define RCPP_IMPTREE_NODE_H


#include <Rcpp.h>
#include <memory>
#include "translation.h"
#include "structs.h"
#include "enums.h"
#include "evaluation.h"

class Node {

  Node *parent_ = nullptr;
  int depth_;
  std::vector<Node*> children_;
  ProbInterval probInt_;
  
  std::vector<int> obsidxs_;
  
  int splitvaridx_ = -1;
  std::vector<int> splitset_;
  
  Rcpp::IntegerVector getNodeObservations(const int variableIndex);
  void calculateProbinterval();
  int calcSplitVariable();
  double calcT(const double maxE, const double minE,
                      const double maxEbase, const double minEbase, 
                      const double maxEposs, const double gamma) const;
  
protected:
  const std::shared_ptr<Data> datap_;
  const std::shared_ptr<Config> configp_;
  double entropy(std::vector<double> x) const;
  
public:
  Node(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp, 
       int depth = 0, Node* parent = nullptr);
  virtual ~Node();
  
  inline bool hasParent() const {return nullptr != this->parent_;}
  
  const std::shared_ptr<Data> getData() const;

  inline Node* getChild(const size_t i) {return children_.at(i);}
  inline Node* getChild(const size_t i) const {return children_.at(i);}
  inline int size() const {return children_.size();}
  
  virtual ProbInterval probabilityInterval(const std::vector<int>& classtable) = 0;
  virtual std::vector<double> maxEntropyDist(const ProbInterval& probint, const bool exact) = 0;
  virtual std::vector<double> minEntropyDist(const ProbInterval& probint) = 0;
  virtual double correctionEntropy(const std::vector<double>& probs, const int n) = 0;
  
  void setSplitVariable(const int idx);
  void setSplitSet(std::vector<int> splitset);
  inline void addSplitObs(const int obsidx) {obsidxs_.push_back(obsidx);}
  
  Evaluation evaluate(const Rcpp::IntegerMatrix & newdata, const Rcpp::List & evalconfig);
  ProbInterval classify(Rcpp::IntegerVector observation);
  
  // Factory Method
  static Node* createNode(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp,
                          int depth, Node* parent = 0);
  
  void makeChildren();
  
  // Summary methods
  int numLeaves() const;
  int numNodes() const;
  void addDepth(std::vector<int> * depths) const;
  
  void printNode(const int parentIdx, const int nsmall, const std::string &sep) const;
  Rcpp::List getNodeByIndex(std::vector<int>& idxs) const;
};

class IDMNode : public Node {
  
  std::vector<double> minVals(const std::vector<double>& array);
  
public:
  IDMNode(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp,
          int depth, Node* parent = nullptr);
  
  ProbInterval probabilityInterval(const std::vector<int>& classtable);
  std::vector<double> maxEntropyDist(const ProbInterval& probint, const bool exact = true);
  std::vector<double> minEntropyDist(const ProbInterval& probint);
  double correctionEntropy(const std::vector<double>& probs, const int n);
};

class NPINode : public Node {
  
  std::vector<double> maxEntropyDistApprox(const ProbInterval& probint);
  std::vector<double> maxEntropyDistExact(const ProbInterval& probint);
  int maxIndexInSet(std::vector<int> array, std::vector<bool> set);
  
public:
  NPINode(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp,
          int depth, Node* parent = nullptr);
  
  ProbInterval probabilityInterval(const std::vector<int>& classtable);
  std::vector<double> maxEntropyDist(const ProbInterval& probint, const bool exact = true);
  std::vector<double> minEntropyDist(const ProbInterval& probint);
  double correctionEntropy(const std::vector<double>& probs, const int n);
};


#endif /* RCPP_IMPTREE_NODE_H */
