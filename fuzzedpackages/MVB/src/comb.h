#ifndef COMB_H_
#define COMB_H_

#include <vector>

namespace lps {
  // class to generate nchoosek combinations
  class comb {
  private:
    // vector to get combination
    std::vector<int> vec;
    // position of pointer
    std::vector<std::vector<int> >::iterator iter;
    const int n, k; // parameters
    // result vector of vectors
    std::vector<std::vector<int> > result;
    // function for recursive calls
    void runComb(int n, int k, int start, std::vector<int> current) {
      using namespace std;
      if (k == 0) {
	result.push_back(current);
	return;
      }
      if (n - start == k) {
	vector<int> tmp(current);
	for (int i = start; i < n; i++)
	  tmp.push_back(vec[i]);
	result.push_back(tmp);
	return;
      }
      current.push_back(vec[start]);
      runComb(n, k - 1, start + 1, current);
      current.pop_back();
      runComb(n, k, start + 1, current);
    }
  public:

  comb(int inN, int inK) : n(inN), k(inK) {
      vec.resize(n);
      for (int i = 0; i < n; i++)
	vec[i] = i;
      std::vector<int> tmp;
      runComb(n, k, 0, tmp);
      iter = result.begin();
    }
  comb(int inN, int inK, const std::vector<int>& input) : 
    n(inN), k(inK) {
      vec = input;
      std::vector<int> tmp;
      runComb(n, k, 0, tmp);
      iter = result.begin();
    }
    bool empty() {
      return iter == result.end();
    }
    std::vector<int>& getNext() {
      return *iter++;
    }
  };
}

#endif // COMB_H_
