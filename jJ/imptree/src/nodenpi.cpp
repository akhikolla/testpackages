// nodenpi.cpp: Imprecise Classification Trees
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

NPINode::NPINode(const std::shared_ptr<Data> datap, const std::shared_ptr<Config> configp, int depth, Node* parent)
  : Node(datap, configp, depth, parent)
{
}

/*
 ** get the index with the (first) maximum value
 ** in an array based on a selection defined in set
 */
int NPINode::maxIndexInSet(std::vector<int> array, std::vector<bool> set) {
  
  int index = -1;
  int max = -1;
  int asize = array.size();
  
  for(int i = 0; i < asize; i++) {
    if(set[i] && array[i] > max) {
      max = array[i];
      index = i;
    }
  }
  return index;
}

std::vector<double> NPINode::maxEntropyDist(const ProbInterval &probint, const bool exact) {
  if(exact) return maxEntropyDistExact(probint);
  return maxEntropyDistApprox(probint);
}


std::vector<double> NPINode::maxEntropyDistApprox(const ProbInterval &probint) {
  
  int ksize = probint.freq.size();
  int maxobs = (*std::max_element(probint.freq.begin(), probint.freq.end())) + 1;
  
  std::vector<int> ks(maxobs, 0);
  std::vector<double> prob(ksize, 0.0);
  for(int i = 0; i < ksize; ++i) {
    ++ks[probint.freq[i]];
  }
  int k0 = ks[0];
  int k1 = ks[1];
  int krem = ksize - k0 - k1;
  int mass = 0;

	double dnobs = static_cast<double>(probint.obs);

	int ni;
	if( krem < k0) {
		for(int i = 0; i < ksize; ++i) {
			ni = probint.freq[i];
			if(ni == 0 || ni == 1) {
				prob[i] = ((krem + (double)k1) / (dnobs * (k0 + k1)));
			} else {
				prob[i] = ((ni - 1.0) / dnobs);
			}
		}
	} else {
		mass = krem - k0;
	  for(int i = 0; i < ksize; ++i) {
	    prob[i] = std::max(probint.freq[i] - 1.0, 1.0) / dnobs;
	  }
	  
		int j = 1, kj, kj1;
		while (mass > 0) {
			
			kj = ks[j];
			kj1 = ks[j + 1];
			if((kj + kj1) < mass) {
				for(int i = 0; i < ksize; ++i){
					ni = probint.freq[i];
					if(ni == j || ni == j + 1) {
						prob[i] += 1.0/dnobs;
						--mass;
					}
				}
			} else {
				for(int i = 0; i < ksize; ++i) {
					ni = probint.freq[i];
					if(ni == j || ni == j + 1) {
						prob[i] += (mass/(dnobs * (kj + kj1)));
					}
				}
				mass = 0;
			  break;
			}
			++j;
			if(j == maxobs) {
			  Rcpp::stop(_("MaxAppox: After all iterations (%i) not all mass has been assigned!\n Remaining mass is: %f\n"), 
                ksize, mass / dnobs);	
			  break;
			}
		}
	}
	return prob;
}

std::vector<double> NPINode::maxEntropyDistExact(const ProbInterval &probint) {

  int ksize = probint.freq.size();
  int maxobs = (*std::max_element(probint.freq.begin(), probint.freq.end())) + 1;
  
  std::vector<int> ks(maxobs, 0);
  std::vector<double> prob(ksize, 0.0);
  
  for(int i = 0; i < ksize; ++i) {
    ++ks[probint.freq[i]];
  }
  int k0 = ks[0];
  int k1 = ks[1];
  int k01 = k0 + k1;
  int krem = ksize - k01;
  
	double dnobs = static_cast<double>(probint.obs);
	

	if(krem == 0) {
	  std::fill(prob.begin(), prob.end(), 1.0 / ksize);
	} else if (k0 > krem) {
    prob = probint.lower;
		int beta = (k01) / (krem + k1);
		int h = (k01) % (krem + k1);
		int j = 0;
		
		if( h < (k1 + 1)) {
			for(int i = 0; i < ksize; ++i) {
			  
			  if(probint.freq[i] <= 1) {
				  if(j < (beta * (krem - 1))) {
					  prob[i] = 1.0 / (dnobs * beta);
				    ++j;
				  } else if(j < k01) {
					  prob[i] = (k1 + 1.0) / (dnobs * (beta * (k1 + 1.0) + h));
				    ++j;
				  } else {
				    Rcpp::stop(_("Something is wrong in calculation"));
				  }
				}
			}	
		} else {
			for(int i = 0; i < ksize; ++i) {
			  if(probint.freq[i] <= 1) {
				  if(j < (h * (beta + 1))) {
					  prob[i] = 1.0 / (dnobs * (beta + 1.0));
				    ++j;
				  } else if(j < k01) {
					  prob[i] = 1.0 / (dnobs * beta);
				    ++j;
				  } else {
				    Rcpp::stop(_("Something is wrong in calculation"));
				  }
				}
			}
		}
	} else {
		int mass = krem - k0;
	  for(int i = 0; i < ksize; ++i) {
	    prob[i] = std::max(probint.freq[i] - 1.0, 1.0) / dnobs;
	  }
		
		int j = 1, Acc, W, ni;
		while(mass > 0) {
			if(ks[j] + ks[j + 1] < mass) {
				for(int i = 0; i < ksize; ++i) {
					ni = probint.freq[i];
					if((ni == j) || (ni == (j + 1))) {
						prob[i] += 1.0 / dnobs;
					}
				}
				mass -= (ks[j] + ks[j + 1]);
			} else {
				W = std::min((mass + 1 + ks[j]), (ks[j] + ks[j + 1]));
				Acc = W;
				for(int i = 0; i < ksize && Acc > 0; ++i){
					ni = probint.freq[i];
					if((ni == j) || (ni == (j + 1))) {
						prob[i] += (mass / (dnobs * W));
						--Acc;
					}
				}
				mass = 0;
				break;
			}
			++j;
			if(j == maxobs) {
			  Rcpp::stop(_("Max: After all iterations (%i) not all mass has been assigned!\n Remaining mass is: %f\n"),
                ksize, static_cast<double>(mass)/dnobs);	
			  break;
			}
		}
		return prob;
	}
	return prob;
}

std::vector<double> NPINode::minEntropyDist(const ProbInterval &probint) {
	
	int ksize = probint.freq.size();
	int nobs = probint.obs;
	
	std::vector<int> lower;
	lower.reserve(ksize);
	std::vector<int> upper;
	upper.reserve(ksize);
	std::vector<double> dlower;
	dlower.reserve(ksize);
	for(int const &freqi : probint.freq) {
	  lower.push_back(std::max(freqi - 1, 0));
	  upper.push_back(std::min(freqi + 1, nobs));
	}
	
	std::vector<bool> set = std::vector<bool>(ksize, true);
	
	int diff, idx, j = 0, mass = nobs - std::accumulate(lower.rbegin(), lower.rend(), 0);
	
	while(mass > 0) {
	  
		idx = maxIndexInSet(lower, set);
		diff = upper[idx] - lower[idx];
		if(diff < mass) {
			lower[idx] = upper[idx];
			set[idx] = false;
			mass -= diff;
		} else {
			lower[idx] = lower[idx] + mass;
			mass = 0;
			break;
		}
		++j;
		if(j == ksize) {
		  Rcpp::stop(_("Min: After all iterations (%i) not all mass has been assigned!\n Remaining mass is: %f\n"), 
               ksize, static_cast<double>(mass) / static_cast<double>(nobs));	
		  break;
		}
	}
	
	for(int i = 0; i < ksize; ++i) {
	  dlower.push_back(lower[i] * 1.0 / static_cast<double>(nobs));
	}
	
  return dlower;
}

double NPINode::correctionEntropy(const std::vector<double>& probs, const int n) {
  if(n > 0) {
    double ent = entropy(probs);
    EntropyCorrection ec = configp_->ec;
    switch(ec) {
    case EntropyCorrection::strobl:
      ent += ((probs.size() - 1.0) / (2.0 * n));
      break;
    case EntropyCorrection::abellan:
      Rcpp::stop(_("Entropy correction 'abellan' not permitted for NPI\n"));
      break;
    default:;
    }
    return ent;
  }
  return -1;
}

ProbInterval NPINode::probabilityInterval(const std::vector<int>& classtable) {
  
  ProbInterval prob;
  prob.obs = std::accumulate(classtable.begin(), classtable.end(), 0);
  prob.freq.clear();
  prob.upper.clear();
  prob.lower.clear();
  double dobs = static_cast<double>(prob.obs);
  for(int classObs : classtable) {
    prob.freq.push_back(classObs);
    prob.upper.push_back((std::min(classObs + 1.0, dobs) / dobs));
    prob.lower.push_back((std::max(classObs - 1.0, 0.0) / dobs));
  }
  //Debug
  //#Rcpp::Rcout << prob.to_string();
  return prob;
}
