// nodeidm.cpp: Imprecise Classification Trees
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

IDMNode::IDMNode(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp, int depth, Node* parent)
  : Node(datap, configp, depth, parent)
{
}

std::vector<double> IDMNode::minVals(const std::vector<double>& array) {
	
/* Initianlizing the minimal values with a not reachable one in the array */
  double nmin = 0.0;
  double min1 = 2.0; 
  double min2 = 2.0;
	bool samemin = true;
	for (const double val : array) {
/* In case a new minimum is identified, the old minimum in 'min1'
** is passed to 'min2' and 'min1' gets the new one */
		if (val < min1) {
			min2 = min1;
			min1 = val;
			nmin = 1;
		} else if(val == min1) {
		  ++nmin;
		} else if(val < min2) {
/* In case a the value is greater than the minimum, 
** but smaller than the second minimal value, then 
** it is assigned to 'min2' */
			min2 = val;
		  samemin = false;
		}
	}
	if(samemin) {
	  min2 = min1;
	}
/* return the minimum and second minimal */
  std::vector<double> res{min1, min2, nmin};
	return res;
}


/* Functions called by R, see description above */

std::vector<double> IDMNode::maxEntropyDist(const ProbInterval &probint, const bool /*exact*/) {
  
  std::vector<double> lower = probint.lower;
  int lsize = lower.size();
  
  double nmin, minval, sminval;
	
  // Due to the nature of IDM the intial free mass is s/(N+s)
  double assignMass, freeMass = configp_->s / (static_cast<double>(probint.obs) + configp_->s);
	
  bool hasFree = true;
	
  // Keep iterating till all assigned
	while(hasFree) {
	  // Vector of minimal values
	  std::vector<double> minvals = minVals(lower);
	  minval = minvals[0];
	  sminval = minvals[1];
	  nmin = minvals[2];

	  if(minval == sminval || !((sminval - minval) < (freeMass / nmin))) {
      // All values in lower have the same value
      //   OR
      // not enough free mass to lift the minimum value(s) to the second minimal value
	    assignMass = (freeMass / nmin);
	    // We have used up all free mass, so signal end
	    hasFree = false;
	  } else {
	    // assign as much mass to lift the minimum value(s) to the second minimal value
	    assignMass = (sminval - minval);
	    // substract the lifts from the free mass
	    freeMass -= (assignMass * nmin);
	  }
	  // Update lower
	  for(int i = 0; i < lsize; ++i) {
	    if(lower[i] == minval) {
	      lower[i] += assignMass;
	    }
	  }
	}
	return lower;
}


std::vector<double> IDMNode::minEntropyDist(const ProbInterval &probint) {

  std::vector<double> lower = probint.lower;
  // get the index with the (first) maximum
  size_t index = std::distance(lower.begin(), std::max_element(lower.begin(), lower.end()));

  // if valid index then set the value of the lower boundary to the 
  // value of the upper one
	if (index < lower.size()) {
		lower[index] = probint.upper[index];
	}
	return lower;
}

double IDMNode::correctionEntropy(const std::vector<double>& probs, const int n) {
  double s = configp_->s;
  if(s > 0 && n > 0) {
    double ent = entropy(probs);
    EntropyCorrection ec = configp_->ec;
    switch(ec) {
    case EntropyCorrection::abellan:
      ent += (s * log2(probs.size())) / (n + s);
      break;
    case EntropyCorrection::strobl:
      ent += ((probs.size() + 1.0) / (2.0 * n + s));
      break;
    default:;
    }
    return ent;
  }
  return -1;
}

ProbInterval IDMNode::probabilityInterval(const std::vector<int>& classtable) {
  
  double s = configp_->s;
  ProbInterval prob;
  prob.obs = std::accumulate(classtable.begin(), classtable.end(), 0);
  prob.freq.clear();
  prob.upper.clear();
  prob.lower.clear();
  double dobs = static_cast<double>(prob.obs);
  for(int classObs : classtable) {
    prob.freq.push_back(classObs);
    prob.upper.push_back((classObs + s) / (dobs + s));
    prob.lower.push_back(static_cast<double>(classObs) / (dobs + s));
  }
  return prob;
}
