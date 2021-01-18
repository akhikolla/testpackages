#include <Rcpp.h>

#include <vector>
#include <algorithm>

using namespace Rcpp;

// node class for the values of the statistics and their initial indices
class NodeStat {
  private:
    double value_;
	  unsigned int index_;

  public:	
	  NodeStat(void) {}
	  NodeStat(const double &value, const unsigned int &index) : value_(value), index_(index) {}
    
    double value() const {
      return value_;
    }
    
    unsigned int index() const {
      return index_;
    }
    
    // definition of operators for sort()
    bool operator > (NodeStat const& rhs) const {
      return value_ > rhs.value();
    }
    
    bool operator >= (NodeStat const& rhs) const {
      return value_ >= rhs.value();
    }
    
    bool operator == (NodeStat const& rhs) const {
      return value_ == rhs.value();
    }
    
    bool operator != (NodeStat const& rhs) const {
      return value_ != rhs.value();
    }
    
    bool operator <= (NodeStat const& rhs) const {
      return value_ <= rhs.value();
    }
    
    bool operator < (NodeStat const& rhs) const {
      return value_ < rhs.value();
    }
};


// [[Rcpp::export(name = ".criticalValuesWeights")]]
NumericVector criticalValuesWeights(const NumericVector &stat, const NumericVector &beta,
                                    const double &alpha) {
	unsigned int d = beta.size();
	unsigned int M = stat.size() / d;

	// creation of a d times M matrix of NodeStats and sorting rowwise
  std::vector< std::vector<NodeStat> > S;
  S.reserve(d);
  for (unsigned int j = 0u; j < d; ++j) {
    std::vector<NodeStat> current;
    current.reserve(M);
    for (unsigned int i = 0u; i < M; ++i) {
      current.push_back(NodeStat(stat[i * d + j], i));
		}
    std::sort(current.begin(), current.end());
    S.push_back(current);
  }
  
	// determination of the starting values
	std::vector<unsigned int> index(d);     // index of the column of the current quantile S[j][index[j]].value
  std::vector<bool> rejection(M, false);  // true if in the column at least one test rejects
  unsigned int level = 0u;                // number of columns in which at least one test rejects

	for (unsigned int j = 0u; j < d; ++j) {
		index[j] = M - static_cast<unsigned int>(alpha * beta[j] * M) - 1u;	// -1 for C-style
		for (unsigned int i = M - 1u; i > index[j]; --i) {
			if (rejection[S[j][i].index()] == false) {
				rejection[S[j][i].index()] = true;
				++level;
			}
		}
	}

	// iteration step
	unsigned int indexMin;  // index with the lowest ratio indexMin:= argmin_{j} (M - index[j]) / beta[j]
	                        // index[indexMin] will be reduced
	do {
    checkUserInterrupt();
		double minimum = (M - index[0u]) / beta[0u]; 
		indexMin = 0u;
		for (unsigned int j = 1u; j < d; ++j) {
			if ((M - index[j]) / beta[j] < minimum) {
				minimum = (M - index[j]) / beta[j];
			  indexMin = j;
			}
		}

		if (rejection[S[indexMin][index[indexMin]].index()] == false) {
			rejection[S[indexMin][index[indexMin]].index()] = true;
			++level;
		}
		--index[indexMin];
	}
	while (level <= static_cast<unsigned int>(M * alpha + 1e-6));
  // revert the last step, since level / M has to be smaller or equal than alpha
	++index[indexMin];

	NumericVector ret = NumericVector(d);
	for (unsigned int j = 0u; j < d; ++j) {
    ret[j] = S[j][index[j]].value();
	}

	return ret;
}
