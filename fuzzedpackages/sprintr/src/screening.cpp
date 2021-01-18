#include <cmath>
#include <queue>
#include <vector>
// we only include RcppArmadillo.h which automatically pulls Rcpp.h in for us
#include "RcppArmadillo.h"

class Triplet
{
private:
  // j, k are the indices of the interaction
  int j, k;
  // the absolute correlation of the j,k interaction with some target
  double abs_cor;
public:
  Triplet(int _j, int _k, double _abs_cor)
  {
    j = _j;
    k = _k;
    abs_cor = _abs_cor;
  }
  int get_j() const{
    return j;
  }
  int get_k() const{
    return k;
  }
  double get_cor() const{
    return abs_cor;
  }
};

// user specified comparator
class Triplet_Comparator{
public:
  bool operator() (const Triplet & p1, const Triplet & p2) const {
    // note that here we use > instead of <
    // so that the resulting heap is a min-heap
    // (by default the C++ priority_queue is a max-heap)
    return p1.get_cor() > p2.get_cor();
  }
};

//' Sure Independence Screening in Step 2
//' 
//' @param x a n-by-p matrix of main effects, with i.i.d rows, and each row represents a vector of observations of p main-effects
//' @param y a vector of length n. In sprinter, y is the residual from step 1
//' @param num_keep the number of candidate interactions in Step 2. Default to be n / [log n]
//' @param square An indicator of whether squared effects should be considered in Step 1 (NOT Step 2!). square == TRUE if squared effects have been considered in Step 1, i.e., squared effects will NOT be considered in Step 2.
//' @param main_effect An indicator of whether main effects should also be screened. Default to be false. The functionality of main_effect = true is not used in sprinter, but for SIS_lasso.
//' @return an matrix of 2 columns, representing the index pair of the selected interactions. 
//' @export
// [[Rcpp::export]]
arma::mat screen_cpp(const arma::mat & x, const arma::vec & y, const int num_keep, const bool square = false, const bool main_effect = false){
  // main is TRUE if we want to screen main effects simultaneously
  // square is TRUE if squared effects have been considered in the first step
  // i.e., squared effects will NOT be considered in the screening step!
  // screening all the interaction x_j*x_k based on its correlation with y
  // return the indices (j, k) of the top k interactions

  const int n = x.n_rows;
  const int p = x.n_cols;
  // j, k are used for indexing interactions
  int j, k;
  int count = 0;
  double cur_cor = 0.0;
  arma::vec temp_vec(n);
  // matrix of num_keep by 2, storing the results
  arma::mat result = arma::zeros<arma::mat>(num_keep, 2);

  // Creates a Min heap of Triplets (ordered by absolute correlation with y)
  std::priority_queue <Triplet, std::vector<Triplet>, Triplet_Comparator> min_heap;

  // calculate the correlations in a stream way
  if(square){
    for(j = 0; j < p - 1; ++j){
      for (k = j + 1; k < p; ++k){
        // this is a little faster than calculating correlation
        // and it gives the same results
        temp_vec = x.col(j) % x.col(k);
        cur_cor = fabs(arma::dot(temp_vec, y)) / arma::stddev(temp_vec);

        // cor function when applied to two vectors, returns a 1-by-1 matrix
        // cur_cor = fabs(arma::as_scalar(arma::cor(x.col(j) % x.col(k), y)));
        if (count < num_keep){
          // use emplace instead of push
          // emplace avoids constructing a temporary Triplet
          // emplace requires C++11, to use it, add
          // // [[Rcpp::plugins(cpp11)]] in the top of the file
          // min_heap.emplace(j, k, cur_cor);
          // found that emplace really doesn't help
          min_heap.push(Triplet(j, k, cur_cor));
          ++count;
        } else{
          if (cur_cor > min_heap.top().get_cor()){
            // pop the current top element
            min_heap.pop();
            // insert the (j,k) interaction
            // min_heap.emplace(j, k, cur_cor);
            min_heap.push(Triplet(j, k, cur_cor));
          }
        }
      }
    }
  }
  else{
    // squared effects are not considered in the first step
    // include them in the screening
    for(j = 0; j < p; ++j){
      for (k = j; k < p; ++k){
        temp_vec = x.col(j) % x.col(k);
        cur_cor = fabs(arma::dot(temp_vec, y)) / arma::stddev(temp_vec);
        // cor function when applied to two vectors, returns a 1-by-1 matrix
        // cur_cor = fabs(arma::as_scalar(arma::cor(x.col(j) % x.col(k), y)));
        if (count < num_keep){
          // use emplace instead of push
          // emplace avoids constructing a temporary Triplet
          // emplace requires C++11, to use it, add
          // // [[Rcpp::plugins(cpp11)]] in the top of the file
          // min_heap.emplace(j, k, cur_cor);
          min_heap.push(Triplet(j, k, cur_cor));
          ++count;
        } else{
          if (cur_cor > min_heap.top().get_cor()){
            // pop the current top element
            min_heap.pop();
            // insert the (j,k) interaction
            //min_heap.emplace(j, k, cur_cor);
            min_heap.push(Triplet(j, k, cur_cor));
          }
        }
      }
    }
  }

  // if we need to screen main effects also
  if(main_effect){
    for(k = 0; k < p; ++k){
      cur_cor = fabs(arma::dot(x.col(k), y)) / arma::stddev(x.col(k));
      if (cur_cor > min_heap.top().get_cor()){
        // pop the current top element
        min_heap.pop();
        // insert the (j,k) interaction
        // for main effects, it corresponds to (0, k)
        min_heap.push(Triplet(-1, k, cur_cor));
      }
    }
  }

  // now store the selected indices of interactions
  for (count = 0; count < num_keep; ++count){
    result(count, 0) = min_heap.top().get_j() + 1;
    result(count, 1) = min_heap.top().get_k() + 1;
    min_heap.pop();
  }
  return result;
}
