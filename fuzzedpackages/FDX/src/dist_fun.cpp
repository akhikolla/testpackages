//'

#include "dist_fun.h"

NumericVector poibinom_int(NumericVector probs, int method, int max_q, bool lower_tail){
  // number of probabilities of success
  int n = probs.length();
  // vector of all probabilities
  NumericVector res;
  if(max_q >= 0) res = NumericVector(std::min<int>(max_q, n) + 1, 1 - (double)lower_tail); else res = NumericVector(n + 1, 1 - (double)lower_tail);
  
  // actual size of results vector
  int size = res.length();
  // loop variable
  int i;
  
  // find zeros and ones in 'probs'
  // numbers of zeros and ones
  int counts0 = 0, counts1 = 0;
  // count zeros and ones, extract other probabilites
  for(i = 0; i < n; i++){
    // current probability of success
    double prob = probs[i];
    // increase 0-counter, if 0 (or lower, for intercepting numerical errors)
    if(prob <= 0) counts0++;
    // increase 1-counter, if 1 (or greater, for intercepting numerical errors)
    else if(prob >= 1) counts1++;
  }
  // reduced vector of probabilities
  NumericVector probs_reduced = probs[(probs > 0) & (probs < 1)];
  // number of values not equal to 0 or 1
  int m = n - counts0 - counts1;
  // starting index of "relevant" range
  int start = counts1;
  // if starting index is greater than 'max_q', return 0-vector
  if(start > max_q) return res;
  // ending index of "relevant" range
  int end = std::min<int>(n - counts0, max_q);
  
  // if there are only zeros and ones, only one outcome is possible
  if(m == 0){
    for(i = start; i < size; i++) res[i] = (double)lower_tail;
  }else{
    // if there is only one probability not equal to zero or one, it is a Bernoulli distribution
    if(m == 1){
      if(lower_tail) res[start] = 1 - probs_reduced[0]; else res[start] = probs_reduced[0];
      for(i = start + 1; i < size; i++) res[i] = (double)lower_tail;
    }else{
      // if all remaining probabilities are equal, it is a Binomial distribution
      bool equal = true;
      Range range(start, end);
      for(i = 1; i < m; i++){
        if(probs_reduced[i] != probs_reduced[0]){
          equal = false;
          break;
        }
      }
      IntegerVector obs(Range(0, end - start));
      // compute probabilites of relevant observations
      if(equal){
        res[range] = pbinom(obs, (double)m, probs_reduced[0], lower_tail);
      }else{
        switch(method){
          //case 0: res[range] = ppb_conv(obs, probs_reduced); break;
          case 0: res[range] = ppb_na(obs, probs_reduced, true, lower_tail); break;
          case 1: res[range] = ppb_dc(obs, probs_reduced, lower_tail); break;
          case 2: res[range] = ppb_gmba(obs, probs_reduced, true, lower_tail); break;
        }
      }
      for(i = end + 1; i < size; i++) res[i] = (double)lower_tail;
    }
  }
  
  return res;
}
/*
NumericVector dpbinom(IntegerVector obs, NumericVector probs, int method = 1){
  // vector of results;
  NumericVector res = pbinom_int(probs, method);
  if(method == 2){
    // compute cumulative probabilites for all possible obervations
    res.push_front(0);
    res = diff(res);
    
    // make sure that all probabilities do not exceed 1
    //res[res < 0] = 0;
    //res[res > 1] = 1;
  }else{
    // make sure that all probabilities sum up to 1
    res = res/sum(res);
  }
  
  // return results
  res = res[obs];
  return res;
}
*/

NumericVector ppbinom(IntegerVector obs, NumericVector probs, int method, bool lower_tail){
  // largest quantile (last element, because 'obs' is sorted in ascending order)
  int max_q = std::max<int>(0, obs[obs.length() - 1]);
  // vector of results;
  NumericVector res = poibinom_int(probs, method, max_q, lower_tail);
  
  // return results
  return res[obs];
}
