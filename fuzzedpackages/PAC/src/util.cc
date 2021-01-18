#include "util.h"

int randsample(const int l, const int u, const vector<double> & p) {
  GetRNGstate();
  int r = unif_rand()*(u-l);
  PutRNGstate();
  return r;
}

uint whichMax(const vector<double>& vec){
  return distance(vec.begin(), max_element(vec.begin(), vec.end()));
}

double vecPartialSum(const vector<double> &v, const uint end){
  double sum = 0.0;
  for(uint i = 0; i < end; i++){
    sum += v[i];
  }
  return sum;
}
double reclnArea(const matrix &r){
  double lnArea = 0;
  for(uint i = 0; i < r.size(); i++){
    lnArea += log(r[i][1] - r[i][0]);
  }
  return lnArea;
}

int nwithin(const matrix &range, const matrix *data){
  int count = 0;
  int n = (*data).size();
  for(int i = 0; i < n; i++){
    if (within(range, (*data)[i])) {
        count++;
    }
  }
  return count;
}

bool within(const matrix &range, const vector<double> & data){
  int dim = data.size();
  for (int d = 0; d < dim; d++) {
      if (data[d] < range[d][0] || data[d] > range[d][1]) {
        return false;
      }
  }
  
  return true;
}

bool int_cmp (int i, int j) {
  return (i==j);
}
uint num_unique(vector<int> vec){
  vector<int>::iterator it;
  std::sort(vec.begin(), vec.end());
  it = std::unique(vec.begin(), vec.end(), int_cmp);
  return distance(vec.begin(), it);
}

