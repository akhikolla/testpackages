#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericVector order_vec(Rcpp::NumericVector& data) {
  Rcpp::NumericVector y(data.size());
  std::iota(y.begin(), y.end(), 0);
  auto comparator = [&data] (int a, int b) { return data[a] < data[b]; }; 
  std::sort(y.begin(), y.end(), comparator);
  
  return y;
}


// [[Rcpp::export]]
Rcpp::NumericVector psrankCpp(Rcpp::NumericVector &data, Rcpp::NumericVector &group, Rcpp::NumericVector &n) {

  double N = data.size();
  double ngroups = n.size();
  Rcpp::NumericVector result(N);
  Rcpp::NumericVector result_ties(N);
  
  
  // Calculate the first pseudo-rank
  result[0] = N/ngroups*1/n[group[0]-1]*1/2 + 0.5;
  
  // define the matrix for the differences between the pseudo-ranks
  NumericMatrix delta(ngroups, ngroups);
  for(int i = 0; i < ngroups; i++){
    for(int j = i; j < ngroups; j++){
      delta(i,j) = N/ngroups*0.5*(1/n[i]+1/n[j]);
      delta(j,i) = delta(i,j);
    }
  }
  
  
  // Case: no ties
  for(int i = 1; i < N; i++){
    result[i] = result[i-1] + delta(group[i]-1, group[i-1]-1);       
  }
  
  double add = 0;
  int j = 0;
  int i = 0;
  result_ties = clone(result);
  
  // Case: ties in the data
  while(i < N - 1){
    
    if(data[i] == data[i+1]) {
      add = 0;
      j = i + 1;
      // sum up the incremental factor for ties
      while(data[i] == data[j]){
        add += 1/n[group[j]-1];
        j++;
        if(j == N) {
          break;
        }
      }
      for(int k = i; k < j; k++){
        // we need to distinguish between i > 0 and i == 0, otherwise result[i-1] not defined
        if(i > 0) {
          result_ties[k] = result[i-1] + delta(group[i]-1, group[i-1]-1) + N/ngroups*0.5*add; // 'mid'-pseudo-ranks
        }
        else {
          result_ties[k] = result[0] +  N/ngroups*0.5*add; // 'mid'-pseudo-ranks
        }
      }
      // resume for loop where last block of ties ended
      i = j-1;
    } // end if
    i++;
  }
  
  
  return result_ties;
}

// [[Rcpp::export]]
Rcpp::NumericVector psrankMinCpp(Rcpp::NumericVector &data, Rcpp::NumericVector &group, Rcpp::NumericVector &n) {
  
  double N = data.size();
  double ngroups = n.size();
  Rcpp::NumericVector result(N);
  Rcpp::NumericVector result_ties(N);
  
  
  // Calculate the first pseudo-rank
  result[0] = 1;
  
  // define the matrix for the differences between the pseudo-ranks
  NumericMatrix delta(ngroups, ngroups);
  for(int i = 0; i < ngroups; i++){
    for(int j = i; j < ngroups; j++){
      delta(i,j) = N/ngroups*(1/n[i]);
      delta(j,i) = N/ngroups*(1/n[j]);
    }
  }
  
  
  // Case: no ties
  for(int i = 1; i < N; i++){
    result[i] = result[i-1] + delta(group[i-1]-1, group[i]-1);       
  }

  double add = 0;
  int j = 0;
  int i = 0;
  result_ties = clone(result);
  
  // Case: ties in the data
  while(i < N - 1){
    
    if(data[i] == data[i+1]) {
      add = 1/n[group[i]-1];
      j = i + 1;
      // sum up the incremental factor for ties
      while(data[i] == data[j]){
        add += 1/n[group[j]-1];
        j++;
        if(j == N) {
          break;
        }
      }
      for(int k = i+1; k < j; k++){
        // we need to distinguish between i > 0 and i == 0, otherwise result[i-1] not defined
          result_ties[k] = result[i];
      }
      if(j < N) {
        result_ties[j] = result[i] + N/ngroups*add;
      }
      // resume for loop where last block of ties ended
      i = j-1;
    } // end if
    i++;
  }
  
  
  return result_ties;
}



// [[Rcpp::export]]
Rcpp::NumericVector psrankMaxCpp(Rcpp::NumericVector &data, Rcpp::NumericVector &group, Rcpp::NumericVector &n) {
  
  double N = data.size();
  double ngroups = n.size();
  Rcpp::NumericVector result(N);
  Rcpp::NumericVector result_ties(N);
  
  
  // Calculate the first pseudo-rank
  result[0] = N/ngroups*1/n[group[0]-1];

  // define the matrix for the differences between the pseudo-ranks
  NumericMatrix delta(ngroups, ngroups);
  for(int i = 0; i < ngroups; i++){
    for(int j = i; j < ngroups; j++){
      delta(i,j) = N/ngroups*(1/n[j]);
      delta(j,i) = N/ngroups*(1/n[i]);
    }
  }
  
  
  // Case: no ties
  for(int i = 1; i < N; i++){
    result[i] = result[i-1] + delta(group[i-1]-1, group[i]-1);       
  }

  double add = 0;
  int j = 0;
  int i = 0;
  result_ties = clone(result);
  
  // Case: ties in the data
  while(i < N - 1){
    
    if(data[i] == data[i+1]) {
      add = 1/n[group[i]-1];
      j = i + 1;
      // sum up the incremental factor for ties
      while(data[i] == data[j]){
        add += 1/n[group[j]-1];
        j++;
        if(j == N) {
          break;
        }
      }
      for(int k = i; k < j; k++){
        // we need to distinguish between i > 0 and i == 0, otherwise result[i-1] not defined
        if(i == 0) {
          result_ties[k] = N/ngroups*add;
        } else {
          result_ties[k] = result_ties[i-1] + N/ngroups*add;
        }
      }
      // resume for loop where last block of ties ended
      i = j-1;
    } // end if
    i++;
  }
  
  
  return result_ties;
}
