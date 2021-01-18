// Rinterface.cpp: Imprecise Classification Trees
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

#include <Rcpp.h>
#include "translation.h"
#include "node.h"
#include "structs.h"
#include "evaluation.h"

// [[Rcpp::export]]
Rcpp::XPtr<Node> treebuilder_cpp(const Rcpp::IntegerMatrix & data, const Rcpp::List & config) {
  
  std::shared_ptr<Data> datap = std::make_shared<Data>(data);
  std::shared_ptr<Config> configp = std::make_shared<Config>();
  configp->s = Rcpp::as<double>(config["s"]);
  configp->gamma = Rcpp::as<double>(config["gamma"]);
  configp->tbase = Rcpp::as<double>(config["tbase"]);
  configp->minbucket = Rcpp::as<int>(config["minbucket"]);
  configp->maxdepth = Rcpp::as<int>(config["depth"]);
  configp->ec = static_cast<EntropyCorrection>(Rcpp::as<int>(config["correction"]));
  configp->sm = static_cast<SplitMetric>(Rcpp::as<int>(config["splitmetric"]));
  configp->ip = static_cast<IpType>(Rcpp::as<int>(config["iptype"]));
  
  Node* root = Node::createNode(datap, configp, 0, 0);
  std::vector<int> splitset;
  for(int i = 0; i < data.ncol(); ++i) {
    if(i != datap->classidx) {
      splitset.push_back(i);
    }
  }
  root->setSplitSet(splitset);
  for(int i = 0; i < data.nrow(); ++i) {
    root->addSplitObs(i);
  }
  root->makeChildren();
  
  Rcpp::XPtr<Node> proot( root, true );
  return proot;
}

// [[Rcpp::export]]
bool hasRoot_cpp(Rcpp::XPtr<Node> prootnode) {
  return !prootnode->hasParent();
}

// [[Rcpp::export]]
Rcpp::List predict_cpp(Rcpp::XPtr<Node> prootnode, 
                   const Rcpp::IntegerMatrix & newdata,
                   const Rcpp::List & evalconfig) {
  
  Evaluation evalresult = prootnode->evaluate(newdata, evalconfig);
  
  return Rcpp::List::create(Rcpp::Named("probintlist") = evalresult.probIntervalList(),
                            Rcpp::Named("classes") = evalresult.predictions(),
                            Rcpp::Named("evaluation") = evalresult.summary()
  );
}

// [[Rcpp::export]]
Rcpp::IntegerVector treeInformation_cpp(Rcpp::XPtr<Node> prootnode) {
  
  
  std::vector<int> * depths = new std::vector<int>();
  prootnode->addDepth(depths);
  std::vector<int>::iterator mresult = std::max_element(depths->begin(), depths->end());
  Rcpp::IntegerVector result = Rcpp::IntegerVector::create(
    Rcpp::Named("depth") = *mresult,
    Rcpp::Named("nleaves") = prootnode->numLeaves(),
    Rcpp::Named("nnodes") = prootnode->numNodes()
  );
  delete depths;
  return result;
}

// [[Rcpp::export]]
void treePrint_cpp(Rcpp::XPtr<Node> prootnode, const int nsmall, const std::string &sep) {
   prootnode->printNode(-1, nsmall, sep);
}

// [[Rcpp::export]]
Rcpp::List getNode_cpp(Rcpp::XPtr<Node> prootnode, Rcpp::IntegerVector idxs) {
  
  std::vector<int> stdidxs = Rcpp::as< std::vector<int> >(idxs);
  std::reverse(stdidxs.begin(), stdidxs.end());
  stdidxs.pop_back();
  return prootnode->getNodeByIndex(stdidxs);
}

// [[Rcpp::export]]
Rcpp::List createProbIntInformation_cpp(const Rcpp::IntegerVector & vec, 
                                        const Rcpp::List & config, 
                                        const bool minentropy, 
                                        const bool maxentropy) {
  std::shared_ptr<Config> configp = std::make_shared<Config>();
  configp->s = Rcpp::as<double>(config["s"]);
  configp->ec = static_cast<EntropyCorrection>(Rcpp::as<int>(config["correction"]));
  configp->ip = static_cast<IpType>(Rcpp::as<int>(config["iptype"]));
  Node* node = Node::createNode(nullptr, configp, 0, nullptr);
 
  Rcpp::List result;
  
  ProbInterval probInt = node->probabilityInterval(Rcpp::as< std::vector<int> >(vec));
  result["probint"] = probInt.toMatrix();
  
  if(maxentropy) {
    std::vector<double> maxEntDist = node->maxEntropyDist(probInt, configp->ip != IpType::npiapprox);
    result["maxEntDist"] = Rcpp::wrap(maxEntDist);
    result["maxEntCorr"] = node->correctionEntropy(maxEntDist, probInt.obs);
  }
  if(minentropy) {
    std::vector<double> minEntDist = node->minEntropyDist(probInt);
    result["minEntDist"] = Rcpp::wrap(minEntDist);
    result["minEntCorr"] = node->correctionEntropy(minEntDist, probInt.obs);
  }
  delete node;
  
  return result;
}
