// node.cpp: Imprecise Classification Trees
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

#include "node.h"

Node::Node(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp, int depth, Node* parent) 
  : parent_(parent), depth_(depth), datap_(datap), configp_(configp)
{
}

Node::~Node() {
  for(int i = 0; i < size(); ++i) {
    delete (Node *) children_[i];
  }
}

Node* Node::createNode(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp, int depth, Node *parent) {
  IpType ipt = configp->ip;
  switch(ipt) {
  case IpType::idm:
    return new IDMNode(datap, configp, depth, parent);
  case IpType::npi:
  case IpType::npiapprox:
    return new NPINode(datap, configp, depth, parent);
  default:
    Rcpp::warning(_("other IPType is not implemented"));
  }
  return nullptr;
}

const std::shared_ptr<Data> Node::getData() const {
  return datap_;
}

void Node::setSplitSet(std::vector<int> splitset) {
  this->splitset_ = splitset;
}

void Node::setSplitVariable(const int idx) {
  splitvaridx_ = idx;
  std::vector<int>::iterator position = std::find(splitset_.begin(), splitset_.end(), idx);
  if (position != splitset_.end()) {
    splitset_.erase(position);
  }
}

void Node::makeChildren() {

  //#Rcpp::Rcout << "current depth reached: " << depth_ << std::endl;
  // Calculate the probability interval of the node
  calculateProbinterval();
  
  // Abandon if there are already children present or max depth was reached
  if(!children_.empty() || depth_ >= configp_->maxdepth) {
    //#Rcpp::Rcout << "max depth reached" << std::endl;
    return; 
  }
  
  // Calculate the index of the splitVariable
  int splitIdx = calcSplitVariable();
  //#Rcpp::Rcout << "splitvaribale index: " <<  splitIdx << std::endl;
  
  // If there is none found, abandon
  if(splitIdx < 0) {
    return;
  }
  // set the splitvarindex and remove from the list of possible split vars
  setSplitVariable(splitIdx);
  
  // Get all possible split points
  int childs = datap_->nlevels[splitvaridx_];
  // Create all childs
  for(int i = 0; i < childs; ++i) {
    Node* child = Node::createNode(datap_ , configp_, depth_+1, this);
    child->setSplitSet(splitset_);
    children_.push_back(child);
  }
  
  // Assign the observations to the childs
  for(unsigned int i = 0; i < obsidxs_.size(); ++i) {
    int obsidx = obsidxs_[i];
    int splitval = datap_->data(obsidx, splitvaridx_);
    getChild(splitval)->addSplitObs(obsidx);

  }
  // Go into recursion
  for(int i = 0; i < size(); ++i) {
    getChild(i)->makeChildren();
  }
}

Rcpp::IntegerVector Node::getNodeObservations(const int variableIndex) {
  Rcpp::IntegerMatrix allObsAllVar = datap_->data;
  Rcpp::IntegerVector allObsByVar = datap_->data(Rcpp::_,variableIndex);
  Rcpp::IntegerVector obsidxRcpp = Rcpp::wrap(obsidxs_);
  Rcpp::IntegerVector nodeObsByVar = allObsByVar[obsidxRcpp];
  return nodeObsByVar;
}

void Node::calculateProbinterval() {

  int classnlevels = datap_->nlevels[datap_->classidx];
  
  Rcpp::IntegerVector classvals = getNodeObservations(datap_->classidx);
  std::vector<int> nodevals(classnlevels, 0);
  for(int i = 0; i < classvals.size(); ++i) {
    nodevals[classvals[i]]++;
  }
  probInt_ = probabilityInterval(nodevals);
  
}

double Node::calcT(const double maxE, const double minE,
                   const double maxEbase, const double minEbase, 
                   const double maxEposs, const double gamma) const {
  
  double Tvalue = 2.0;
  
  if(maxE < maxEbase) {
    
    double pessimism = (maxEposs - maxEbase) / (maxEposs - maxE);
    double optimism = (minEbase - minE) / (maxE + fabs(minE - minEbase));
    // Minus in Hurvitz-like value due to partial negative entropy values in fraction
    Tvalue =  gamma * pessimism - (1.0 - gamma) * optimism;
    
    if(maxE < minEbase) {
      // T may get at minimum -1 so with -3 we are guaranteed to have it as splitting candidate
      Tvalue -= 3.0;
    }
  }
  return Tvalue;
}

double Node::entropy(std::vector<double> x) const {
  
  std::transform(x.begin(), x.end(), x.begin(),
                 [](double d) -> double {
                   return ((d > 0.0) ? (d * std::log2(d)) : 0.0);
                 });
  return -std::accumulate(x.begin(), x.end(), 0.0);
}

int Node::calcSplitVariable() {
  
  //#Rcpp:: Rcout << "Begin: calcSplitVariable" << std::endl;
  // access all config variables
  double gamma = configp_->gamma;
  int minbucket = configp_->minbucket;
  bool exact = (configp_->ip != IpType::npiapprox);
  
  //#Rcpp:: Rcout << "Begin: entropy in base" << std::endl << probInt_.to_string(6) << std::endl;
  // calculate the entropies in the root node
  double maxEntBase = correctionEntropy(maxEntropyDist(probInt_, exact), probInt_.obs);
  double minEntBase = correctionEntropy(minEntropyDist(probInt_), probInt_.obs);
  std::vector<double> unifDist(probInt_.lower.size(), 1.0 / probInt_.lower.size());
  double maxEntPoss = correctionEntropy(unifDist, probInt_.obs);
  
  // Number of levels
  //#Rcpp::Rcout << "splits to consider: " << splitset_.size() << std::endl;
  
  // Initialise vectors
  std::vector<double> tvalues;
  std::vector<int> toRemoveFromVarList;
  
  int classnlevels = datap_->nlevels[datap_->classidx];
  Rcpp::IntegerVector classvals = getNodeObservations(datap_->classidx);
  
  for(const int splitVarIndex : splitset_) {
    
    //#Rcpp::Rcout << "Consider split: " << splitVarIndex << std::endl;
    
    // Observations for that variable in the node
    Rcpp::IntegerVector col = getNodeObservations(splitVarIndex);
    
    // Unique values
    Rcpp::IntegerVector uc = unique(col);
    // number of levels of variable
    int splitVarNlevels = datap_->nlevels[splitVarIndex];
    
    //#Rcpp::Rcout << "subsplits to consider: " << splitVarNlevels << std::endl;
    
    // if not all levels are present, then variable is
    // no longer a split candidate and added to 
    // list to remove from further candidates
    
    //#Rcpp::Rcout << "uc: " << uc << std::endl << "size uc: " << uc.size() << "; splitvarnlevels: " << splitVarNlevels << std::endl;
    
    if(uc.size() < splitVarNlevels) {
      toRemoveFromVarList.push_back(splitVarIndex);
      //#Rcpp::Rcout << "not all present -> next var" << std::endl;
      continue;
    }
    
    // Initialise vectors for entropies in subnodes
    std::vector<double> maxEnt;
    std::vector<double> minEnt;
    maxEnt.reserve(splitVarNlevels);
    minEnt.reserve(splitVarNlevels);
    std::vector<int> splitVarLevels(splitVarNlevels);
    std::iota(splitVarLevels.begin(), splitVarLevels.end(), 0);
    
    // Bool storing if sufficient obs are in subnodes
    bool sufficientObs = true;
    
    // Iterate over all subnodes
    for(int j = 0; j < splitVarNlevels; ++j) {
      
      //#Rcpp::Rcout << "subnode with index: " << j << std::endl;
      // create the vector of class values in subnode
      std::vector<int> vals(classnlevels, 0);
      for(int k = 0; k < col.size(); ++k) {
        if(col[k] == splitVarLevels[j]) {
          ++vals[classvals[k]];
        }
      }
      // Disqualify variable as candidate if subnode has
      // too few observations: Added it to list to remove 
      // and signal to outer loop via 'sufficientObs'
      int nodeObs = std::accumulate(vals.begin(), vals.end(), 0);
      if(nodeObs < minbucket) {
        //#Rcpp::Rcout << nodeObs << " vs " << minbucket << ": remove index: " << splitVarIndex << std::endl;
        toRemoveFromVarList.push_back(splitVarIndex);
        sufficientObs = false;
        break;
      }
      
      //#Rcpp::Rcout << "subnode has passed minbucket test" << std::endl;
      //#Rcpp::Rcout << "vals: " << vals << std::endl;
      // Perform entropy calculation in the subnode
      ProbInterval lprobi = probabilityInterval(vals);
      
      //#Rcpp::Rcout << "constructed probint" << std::endl << lprobi.to_string() << std::endl;
      
      // weight the entropy values with the relative frequencies of falling into node
      maxEnt.push_back(correctionEntropy(maxEntropyDist(lprobi, exact), lprobi.obs) * lprobi.obs / probInt_.obs);
      //#Rcpp::Rcout << "constructed maxEnt" << std::endl;
      minEnt.push_back(correctionEntropy(minEntropyDist(lprobi), lprobi.obs) * lprobi.obs / probInt_.obs);
      //#Rcpp::Rcout << "constructed minEnt" << std::endl;
    }
    
    // Proceed if still enough observations in subnode
    if(sufficientObs) {
      
      //#Rcpp::Rcout << "enough observations -> calculate min max" << std::endl;
      double lowerEnt = std::accumulate(minEnt.begin(), minEnt.end(), 0.0);
      double upperEnt = std::accumulate(maxEnt.begin(), maxEnt.end(), 0.0);
      tvalues.push_back(calcT(upperEnt, lowerEnt, maxEntBase, minEntBase, maxEntPoss, gamma));
    }
  }
  
  //#Rcpp::IntegerVector vec1 = Rcpp::wrap(splitset_);
  //#Rcpp::Rcout << vec1 << std::endl;
  
  // Remove already disqualified variables
  if(!toRemoveFromVarList.empty()) {
    //#Rcpp::IntegerVector vec2 = Rcpp::wrap(toRemoveFromVarList);
    //#Rcpp::Rcout << vec2 << std::endl;
    std::vector<int> splitset_new(splitset_.size());
    std::vector<int>::iterator it = std::set_difference (splitset_.begin(), splitset_.end(),
                                                         toRemoveFromVarList.begin(), toRemoveFromVarList.end(),
                                                         splitset_new.begin());
    splitset_new.resize(it-splitset_new.begin());
    splitset_ = splitset_new;
    tvalues.shrink_to_fit();
    
    //#Rcpp::IntegerVector vec3 = Rcpp::wrap(splitset_);
    //#Rcpp::Rcout << vec3 << " ; " << splitset_.size()<< std::endl;
    
    //#Rcpp::NumericVector vec4 = Rcpp::wrap(tvalues);
    //#Rcpp::Rcout << vec4 << " ; " << tvalues.size()<< std::endl;
  }
  
  int splitIdx = -1;
  
  // If no variable qualifies for splitting return directly
  if(splitset_.size() == 1 && tvalues[0] < configp_->tbase) {
    splitIdx = splitset_[0];
  } else if(splitset_.size() > 1) {
    // get the minimal t-value 
    std::vector<double>::iterator minTvalueIt = std::min_element(tvalues.begin(), tvalues.end());
    // do nothing if minimal t-value is in range of nogo-split
    if (*minTvalueIt < configp_->tbase) {
      if(std::count(tvalues.begin(), tvalues.end(), *minTvalueIt) > 1) {
        // get the index with lowest t-value; draw one randomly if multiple
        Rcpp::NumericVector tvaluesr = Rcpp::wrap(tvalues);
        //#Rcpp::Rcout << "tvalues: " << tvaluesr << std::endl;
        Rcpp::IntegerVector ids = Rcpp::seq_along(tvaluesr) - 1;
        Rcpp::IntegerVector minids = ids[tvaluesr == *minTvalueIt];
        splitIdx = splitset_[Rcpp::sample(minids, 1)[0]];
        //#Rcpp::Rcout<< "multiple minimal tvalues" << std::endl;
      } else {
        //#Rcpp::Rcout<< "one minimal tvalue" << std::endl;
        splitIdx = splitset_[minTvalueIt-tvalues.begin()];
      }
    }
  }
  //#Rcpp::Rcout << "splitidx: " << splitIdx << std::endl;
  return splitIdx;
}

Evaluation Node::evaluate(const Rcpp::IntegerMatrix & newdata, const Rcpp::List & evalconfig) {
  
  if(hasParent()) {
    Rcpp::stop(_("evaluation only allowed from root node"));
  }
  int nObs = newdata.rows();
  std::vector<ProbInterval> probInts;
  probInts.reserve(nObs);
  
  for(int i = 0; i < nObs; ++i) {
    Rcpp::IntegerMatrix::ConstRow row = newdata.row(i);
    probInts.push_back(this->classify(row));
  }
  double utility = Rcpp::as<double>(evalconfig["utility"]);
  Dominance dom = static_cast<Dominance>(Rcpp::as<int>(evalconfig["dominance"]));
  Data observation = Data(newdata);
  
  Evaluation eval = Evaluation(utility, dom, probInts, newdata);
  return eval;
}

ProbInterval Node::classify(Rcpp::IntegerVector observation) {
  if(splitvaridx_ < 0) {
    return probInt_;
  }
  return children_.at(observation[splitvaridx_])->classify(observation);
}


/**
 * Returns number of leaves in tree structure below node
 *
 * @return number of leaves as positive int
 */
int Node::numLeaves() const {
  
  int num = 0;
  if (splitvaridx_ < 0) {
    return 1;
  } else {
    for (Node * child : children_) {
      num += child->numLeaves();
    }
  }
  return num;
}

/**
 * Returns number of nodes in tree structure below node
 *
 * @return number of nodes as positive int
 */
int Node::numNodes() const {
  
  int no = 1;
  if (splitvaridx_ > -1) {
    for (Node * child : children_) {
      no += child->numNodes();
    }
  }
  return no;
}


void Node::addDepth(std::vector<int> * depths) const {
  
  if(splitvaridx_ > -1) {
    for (Node * child : children_) {
      child->addDepth(depths);
    }
  } else {
    depths->push_back(depth_);
  }
}


void Node::printNode(int parentIdx, const int nsmall, const std::string & sep) const {
  
  //Print the depth
  Rcpp::Rcout << "(" << depth_ << ") ";
  
  int nodedepth = depth_;  
  while(nodedepth > 0) {
    Rcpp::Rcout << "  ";
    --nodedepth;
  }
  
  // Print type of node
  if(nullptr == parent_) {
    Rcpp::Rcout << _("root: ");
  } else {
    int parent_splitIdx = parent_->splitvaridx_;
    Rcpp::CharacterVector parent_labels = Rcpp::as<Rcpp::CharacterVector>(datap_->labels.at(parent_splitIdx));
    Rcpp::Rcout << datap_->varnames[parent_splitIdx] << "=" << parent_labels[parentIdx] << ": ";
  }
  
  //Print size of observations and probabilities
  Rcpp::Rcout << "n=" << obsidxs_.size() << " (" << probInt_.to_string(nsmall, sep) << ")";
  
  // Print end of line (stared for leave) and recurse
  if(splitvaridx_ > -1) {
    Rcpp::Rcout << std::endl;
    for (int i = 0; i < size(); ++i) {
      children_[i]->printNode(i, nsmall, sep);
    }
  } else {
    Rcpp::Rcout << " *" << std::endl;
  }
}

Rcpp::List Node::getNodeByIndex(std::vector<int>& idxs) const {
  
  if(this->size() > 0 && idxs.size() > 0) {
    int idxc = idxs.back();
    idxs.pop_back();
    if(idxc >= this->size()) {
      Rcpp::stop(Rcpp::sprintf<75>(_("Queried index (%d) > child size (%d)"), idxc + 1, this->size()));
    }
    return this->getChild(idxc)->getNodeByIndex(idxs);
  } else if (idxs.size() > 0 && this->size() == 0) {
    Rcpp::stop(_("Too deep recursion: No nodes available further down!"));
  } else {
    Rcpp::NumericMatrix probint = this->probInt_.toMatrix();
    Rcpp::colnames(probint) = Rcpp::wrap<Rcpp::CharacterVector>(this->getData()->labels[this->getData()->classidx]);
    Rcpp::List result;
    result["probint"] = probint;
    result["depth"] = this->depth_;
    if(this->splitvaridx_ > -1) {
      result["splitter"] = static_cast<const char*>(this->getData()->varnames.at(this->splitvaridx_));
    } else {
      result["splitter"] = Rcpp::CharacterVector::create(NA_STRING);
    }
    result["children"] = this->size();
    result["traindataIdx"] = Rcpp::wrap(this->obsidxs_);
    IpType ipt = this->configp_->ip;
    Rcpp::List iplist;
    iplist["iptype"] = IpTypeLookup::toString(ipt);
    if(ipt == IpType::idm) {
      iplist["s"] = this->configp_->s;
    }
    result["ipmodel"] = iplist;
    return result;
  }
}
