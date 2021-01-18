#include "GiniFunctions.h"
using namespace std;
using namespace Rcpp;

#include <cmath>
#include <algorithm>


// [[Rcpp::export]]
double VectorSum(NumericVector x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

// [[Rcpp::export]]
double VectorVar(NumericVector x) {
  double mean, sum = 0;
  int n = x.size();
  mean = VectorSum(x)/n;
  for(int i=0; i<n; i++){
    sum += pow(x[i] - mean, 2.0);
  }
  return sum*(n-1)/n;  
}

// [[Rcpp::export]]
IntegerVector orderc(NumericVector x){
  // calling order()
  Function f("order");   
  // order(x)
  return f(x);
}

// [[Rcpp::export]]
NumericMatrix ss(NumericMatrix X_, IntegerVector ind_) {
  
  int n = X_.nrow(), k = X_.ncol();
  
  arma::mat X(X_.begin(), n, k, false);
  arma::uvec ind = as<arma::uvec>(ind_);
  arma::mat submat = X.rows(ind);
  
  return wrap(submat);
}


// generic function for EuclideanDistance
template <typename InputIterator1, typename InputIterator2>
inline double EuclideanDistance(InputIterator1 begin1, InputIterator1 end1,
                                InputIterator2 begin2) {
  
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    
    // accumulate if appropirate
    rval += (d1-d2)*(d1-d2);
  }
  return sqrt(rval);
}
// generic function for EuclideanDistance
template <typename InputIterator1>
inline double Sum(InputIterator1 begin1, InputIterator1 end1) {
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double d1 = *it1++;
    
    // accumulate sum
    rval += d1;
  }
  return rval;
}

// [[Rcpp::export]]
NumericMatrix rcpp_Kernel_Distance(NumericVector x, double sigma) {
  
  // allocate the matrix we will return
  double n = x.size();
  NumericMatrix rmat(n, n);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      
      // write to output matrix
      rmat(i,j) = sqrt(2-2*exp(-std::abs(x[i]-x[j])/sigma));
      rmat(j,i) = rmat(i,j);
    }
  }
  
  return rmat;
}


// [[Rcpp::export]]
NumericMatrix rcpp_Eu_distance(NumericMatrix mat) {
  
  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow(), mat.nrow());
  
  for (int i = 0; i < rmat.nrow(); i++) {
    for (int j = 0; j < i; j++) {
      
      // rows we will operate on
      NumericMatrix::Row row1 = mat.row(i);
      NumericMatrix::Row row2 = mat.row(j);
      
      // write to output matrix
      rmat(i,j) = EuclideanDistance(row1.begin(), row1.end(), row2.begin());
    }
  }
  
  return rmat;
}

// [[Rcpp::export]]
double rcpp_covg(NumericVector x, NumericVector y) {
  int n = x.size();
  NumericVector ind = y-1;
  if(n==1)
    return 0;
  double Covariance = 0;
  double r = -n+1;
  for(int i = 0; i < n; i++) {
    Covariance += r*x[ind[i]];
    r =r+2;
  }
  return abs(2*Covariance/(n*(n-1)));
}

// [[Rcpp::export]]
double rcpp_Kcovg(NumericVector x, NumericVector y, double sigma) {
  int n = x.size();
  double KDistance=0;
  NumericVector ind = y-1;
  if(n==1)
    return 0;
  // sum every element
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      KDistance += sqrt(2-2*exp(-std::abs(x[ind[i]]-x[ind[j]])/sigma));
    }
  }
  return abs(2*KDistance/(n*(n-1)));
}

// [[Rcpp::export]]
double rcpp_Kcovg_alpha(NumericVector x, NumericVector y, double sigma, double alpha) {
  int n = x.size();
  double KDistance=0;
  NumericVector ind = y-1;
  if(n==1)
    return 0;
  // sum every element
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      KDistance += sqrt(2-2*exp(-pow(std::abs(x[ind[i]]-x[ind[j]]),alpha)/sigma));
    }
  }
  return (abs(2*KDistance/(n*(n-1))));
}


// [[Rcpp::export]]
double rcpp_covg_alpha(NumericVector x, NumericVector y, double alpha) {
  int n = x.size();
  double Distance=0;
  NumericVector ind = y-1;
  if(n==1)
    return 0;
  // sum every element
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      Distance += pow(std::abs(x[ind[i]]-x[ind[j]]),alpha);
    }
  }
  return (2*Distance/(n*(n-1)));
}

// [[Rcpp::export]]
double Rcpp_Covg(NumericMatrix x, NumericVector y) {
  
  double Covariance = 0;
  int n = x.nrow();
  int m = x.ncol();
  IntegerVector ind = orderc(y)-1;
  if (y.size()!=n){
    Rcout << "x and y must be the same size"<< std::endl;
  }
  else {
    if(m==1)
      return 0;
    else{
      
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          
          // rows we will operate on
          NumericMatrix::Row row1 = x.row(ind[i]);
          NumericMatrix::Row row2 = x.row(ind[j]);
          
          // write to output matrix
          Covariance += EuclideanDistance(row1.begin(), row1.end(), row2.begin());
        }
      }
      
      
    }
  }
  return 2*Covariance/(n*(n-1));
}

// [[Rcpp::export]]
double Rcpp_Covg_Alpha(NumericMatrix x, NumericVector y, double alpha) {
  
  double Covariance = 0;
  int n = x.nrow();
  int m = x.ncol();
  IntegerVector ind = orderc(y)-1;
  if (y.size()!=n){
    Rcout << "x and y must be the same size"<< std::endl;
  }
  else {
    if(m==1)
      return 0;
    else{
      
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          
          // rows we will operate on
          NumericMatrix::Row row1 = x.row(ind[i]);
          NumericMatrix::Row row2 = x.row(ind[j]);
          
          // write to output matrix
          Covariance += pow(EuclideanDistance(row1.begin(), row1.end(), row2.begin()),alpha);
        }
      }
    }
  }
  return (2*Covariance/(n*(n-1)));
}



// [[Rcpp::export]]
double Rcpp_KCovg(NumericMatrix x, NumericVector y, double sigma) {
  
  double KCovariance = 0;
  int n = x.nrow();
  int m = x.ncol();
  IntegerVector ind = orderc(y)-1;
  if (y.size()!=n){
    Rcout << "x and y must be the same size"<< std::endl;
  }
  else {
    if(m==1){
      NumericVector x1 = x.column(0);
      return rcpp_Kcovg(x1, y, sigma); 
    } else{
      
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          
          // rows we will operate on
          NumericMatrix::Row row1 = x.row(ind[i]);
          NumericMatrix::Row row2 = x.row(ind[j]);
          
          // write to output matrix
          KCovariance += sqrt(2-2*exp(-EuclideanDistance(row1.begin(), row1.end(), row2.begin())/sigma));
        }
      }
    }
  }
  return (2*KCovariance/(n*(n-1)));
}


// [[Rcpp::export]]
double Rcpp_KCovg_Alpha(NumericMatrix x, NumericVector y, double sigma, double alpha) {
  
  double KCovariance = 0;
  int n = x.nrow();
  int m = x.ncol();
  IntegerVector ind = orderc(y)-1;
  if (y.size()!=n){
    Rcout << "x and y must be the same size"<< std::endl;
  }
  else {
    if(m==1){
      NumericVector x1 = x.column(0);
      return rcpp_Kcovg_alpha(x1, y, sigma, alpha);     
    } else{
      
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          
          // rows we will operate on
          NumericMatrix::Row row1 = x.row(ind[i]);
          NumericMatrix::Row row2 = x.row(ind[j]);
          
          // write to output matrix
          KCovariance += sqrt(2-2*exp(-pow(EuclideanDistance(row1.begin(), row1.end(), row2.begin()),alpha)/sigma));
        }
      }
    }
  }
  return (2*KCovariance/(n*(n-1)));
}


// [[Rcpp::export]]
double rcpp_gCov(NumericVector x, NumericVector y) {
  
  int n = x.size();
  NumericVector ind;
  NumericVector delta(n);
  int group=0, n1=0, n2=0;
  double max = y[0];
  double Covariance;
  
  for(int i=0; i<n; i++){
    if (y[i] == max){
      n2=i;
    } else {
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]= rcpp_covg(x[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
      n1=i;
      max = y[n2+1];
      group++;
    }
    if(i==(n-1)){
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]=rcpp_covg(x[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
    }
  }
  ind = orderc(x);
  Covariance = rcpp_covg(x,ind)-(VectorSum(delta));
  return(abs(Covariance));
}

// [[Rcpp::export]]
double rcpp_gCov_alpha(NumericVector x, NumericVector y, double alpha) {
  
  int n = x.size();
  NumericVector ind;
  NumericVector delta(n);
  int group=0, n1=0, n2=0;
  double max = y[0];
  double Covariance;
  
  for(int i=0; i<n; i++){
    if (y[i] == max){
      n2=i;
    } else {
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]= rcpp_covg_alpha(x[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
      n1=i;
      max = y[n2+1];
      group++;
    }
    if(i==(n-1)){
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]=rcpp_covg_alpha(x[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
    }
  }
  ind = orderc(x);
  Covariance = rcpp_covg_alpha(x,ind, alpha)-(VectorSum(delta));
  return(abs(Covariance));
}

// [[Rcpp::export]]
double rcpp_gCor_alpha(NumericVector x, NumericVector y, double alpha) {
  
  int n = x.size();
  NumericVector ind;
  NumericVector delta(n);
  int group=0, n1=0, n2=0;
  double max = y[0];
  double Correlation;
  
  for(int i=0; i<n; i++){
    if (y[i] == max){
      n2=i;
    } else {
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]= rcpp_covg_alpha(x[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
      n1=i;
      max = y[n2+1];
      group++;
    }
    if(i==(n-1)){
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]=rcpp_covg_alpha(x[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
    }
  }
  ind = orderc(x);
  Correlation = 1- VectorSum(delta)/rcpp_covg_alpha(x,ind, alpha);
  return Correlation;
}



// [[Rcpp::export]]
double rcpp_gCor(NumericVector x, NumericVector y) {
  
  int n = x.size();
  NumericVector ind;
  NumericVector delta(n);
  int group=0, n1=0, n2=0;
  double max = y[0];
  double Correlation;
  
  for(int i=0; i<n; i++){
    if (y[i] == max){
      n2=i;
    } else {
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]= rcpp_covg(x[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
      n1=i;
      max = y[n2+1];
      group++;
    }
    if(i==(n-1)){
      ind = orderc(x[Rcpp::Range(n1, n2)]);
      delta[group]=rcpp_covg(x[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
    }
  }
  ind = orderc(x);
  Correlation = 1- VectorSum(delta)/rcpp_covg(x,ind);
  return Correlation;
}


// [[Rcpp::export]]
double Rcpp_gCov(NumericMatrix x) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_covg(z[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_covg(z[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Covariance = rcpp_covg(z,ind)-(VectorSum(delta));
    return(abs(Covariance));
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_Covg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)])*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_Covg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)])*(n2-n1+1)/n;
      }
    }
    Covariance = Rcpp_Covg(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y)-(VectorSum(delta));
    return abs(Covariance);
  }
}

// [[Rcpp::export]]
double Rcpp_gCov_Alpha(NumericMatrix x, double alpha) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_covg_alpha(z[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_covg_alpha(z[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Covariance = rcpp_covg_alpha(z,ind, alpha)-(VectorSum(delta));
    return(abs(Covariance));
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_Covg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_Covg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], alpha)*(n2-n1+1)/n;
      }
    }
    Covariance = Rcpp_Covg_Alpha(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y, alpha)-(VectorSum(delta));
    return abs(Covariance);
  }
}

// [[Rcpp::export]]
double Rcpp_KgCov_Alpha(NumericMatrix x,double sigma, double alpha) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_Kcovg_alpha(z[Rcpp::Range(n1, n2)],ind, sigma, alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_Kcovg_alpha(z[Rcpp::Range(n1, n2)],ind, sigma, alpha)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Covariance = rcpp_Kcovg_alpha(z,ind,sigma, alpha)-(VectorSum(delta));
    return(abs(Covariance));
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_KCovg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], sigma, alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_KCovg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], sigma, alpha)*(n2-n1+1)/n;
      }
    }
    Covariance = Rcpp_KCovg_Alpha(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y, sigma, alpha)-(VectorSum(delta));
    return abs(Covariance);
  }
}


// [[Rcpp::export]]
double Rcpp_KgCor_Alpha(NumericMatrix x,double sigma, double alpha) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_Kcovg_alpha(z[Rcpp::Range(n1, n2)],ind, sigma, alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_Kcovg_alpha(z[Rcpp::Range(n1, n2)],ind, sigma, alpha)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Correlation = 1- VectorSum(delta)/rcpp_Kcovg_alpha(z,ind,sigma, alpha);
    return Correlation;
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_KCovg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], sigma, alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_KCovg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], sigma, alpha)*(n2-n1+1)/n;
      }
    }
    Correlation = 1- VectorSum(delta)/Rcpp_KCovg_Alpha(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y, sigma, alpha);
    return Correlation;
  }
}


// [[Rcpp::export]]
double Rcpp_gCor_Alpha(NumericMatrix x, double alpha) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_covg_alpha(z[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_covg_alpha(z[Rcpp::Range(n1, n2)],ind, alpha)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Correlation = 1- VectorSum(delta)/rcpp_covg_alpha(z,ind, alpha);
    return Correlation;
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_Covg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], alpha)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_Covg_Alpha(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)], alpha)*(n2-n1+1)/n;
      }
    }
    Correlation = 1- VectorSum(delta)/Rcpp_Covg_Alpha(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y, alpha);
    return Correlation;
  }
}


// [[Rcpp::export]]
double Rcpp_KgCov(NumericMatrix x, double sigma) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_Kcovg(z[Rcpp::Range(n1, n2)],ind,sigma)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_Kcovg(z[Rcpp::Range(n1, n2)],ind,sigma)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Covariance = rcpp_Kcovg(z,ind,sigma)-(VectorSum(delta));
    return(abs(Covariance));
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Covariance;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_KCovg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)],sigma)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_KCovg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)],sigma)*(n2-n1+1)/n;
      }
    }
    Covariance = Rcpp_KCovg(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y,sigma)-(VectorSum(delta));
    return(abs(Covariance));
  }
  
}

// [[Rcpp::export]]
double Rcpp_gCor(NumericMatrix x) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_covg(z[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_covg(z[Rcpp::Range(n1, n2)],ind)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Correlation = 1- VectorSum(delta)/rcpp_covg(z,ind);
    return(Correlation);
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_Covg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)])*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_Covg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)])*(n2-n1+1)/n;
      }
    }
    Correlation = 1- VectorSum(delta)/Rcpp_Covg(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y);
    return(Correlation);
  }
  
}

// [[Rcpp::export]]
double Rcpp_KgCor(NumericMatrix x, double sigma) {
  
  int m=x.ncol();
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]= rcpp_Kcovg(z[Rcpp::Range(n1, n2)],ind,sigma)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        ind = orderc(z[Rcpp::Range(n1, n2)]);
        delta[group]=rcpp_Kcovg(z[Rcpp::Range(n1, n2)],ind,sigma)*(n2-n1+1)/n;
      }
    }
    ind = orderc(z);
    Correlation = 1- VectorSum(delta)/rcpp_Kcovg(z,ind,sigma);
    return Correlation;
  } else {
    int n = x.nrow();
    NumericVector y= x.column(m-1) ;
    NumericVector delta(n);
    int group=0, n1=0, n2=0;
    double max = y[0];
    double Correlation;
    
    for(int i=0; i<n; i++){
      if (y[i] == max){
        n2=i;
      } else {
        delta[group]= Rcpp_KCovg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)],sigma)*(n2-n1+1)/n;
        n1=i;
        max = y[n2+1];
        group++;
      }
      if(i==(n-1)){
        delta[group]= Rcpp_KCovg(x(Rcpp::Range(n1, n2),Rcpp::Range(0,m-2)), y[Rcpp::Range(n1, n2)],sigma)*(n2-n1+1)/n;
      }
    }
    Correlation = 1- VectorSum(delta)/Rcpp_KCovg(x(Rcpp::Range(0,n-1),Rcpp::Range(0,m-2)), y,sigma);
    return Correlation;
  }
  
}

// [[Rcpp::export]]
double Rcpp_HatV_gCov(NumericMatrix x) {
  
  int m=x.ncol();
  double mean, sum = 0;
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector rhoHat(n);
    
    for(int i=0; i<n; i++){
      NumericVector x1 = clone(z);
      x1.erase(i);
      NumericVector y1=  clone(y);
      y1.erase(i);
      rhoHat[i]=rcpp_gCov(x1,y1);
    }
    mean = VectorSum(rhoHat)/n;
    for(int i=0; i<n; i++){
      sum += pow(rhoHat[i] - mean, 2.0);
    }
    return sum*(n-1)/n;
  } else {
    int n = x.nrow();
    NumericVector rhoHat(n);
    IntegerVector ind=Range(0, n-1);
    for(int i=0; i<n; i++){
      IntegerVector ind1=clone(ind);
      ind1.erase(i);
      NumericMatrix x1 = ss(x, ind1);
      rhoHat[i]=Rcpp_gCov(x1);
    }
    
    return VectorVar(rhoHat);
  }
  
}

// [[Rcpp::export]]
double Rcpp_HatV_gCov_Alpha(NumericMatrix x, double alpha) {
  
  int m=x.ncol();
  double mean, sum = 0;
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector rhoHat(n);
    
    for(int i=0; i<n; i++){
      NumericVector x1 = clone(z);
      x1.erase(i);
      NumericVector y1=  clone(y);
      y1.erase(i);
      rhoHat[i]=rcpp_gCov_alpha(x1,y1, alpha);
    }
    mean = VectorSum(rhoHat)/n;
    for(int i=0; i<n; i++){
      sum += pow(rhoHat[i] - mean, 2.0);
    }
    return sum*(n-1)/n;
  } else {
    int n = x.nrow();
    NumericVector rhoHat(n);
    IntegerVector ind=Range(0, n-1);
    for(int i=0; i<n; i++){
      IntegerVector ind1=clone(ind);
      ind1.erase(i);
      NumericMatrix x1 = ss(x, ind1);
      rhoHat[i]=Rcpp_gCov_Alpha(x1, alpha);
    }
    
    return VectorVar(rhoHat);
  }
  
}

// [[Rcpp::export]]
double Rcpp_HatV_gCor_Alpha(NumericMatrix x, double alpha) {
  
  int m=x.ncol();
  double mean, sum = 0;
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector rhoHat(n);
    
    for(int i=0; i<n; i++){
      NumericVector x1 = clone(z);
      x1.erase(i);
      NumericVector y1=  clone(y);
      y1.erase(i);
      rhoHat[i]=rcpp_gCor_alpha(x1,y1, alpha);
    }
    mean = VectorSum(rhoHat)/n;
    for(int i=0; i<n; i++){
      sum += pow(rhoHat[i] - mean, 2.0);
    }
    return sum*(n-1)/n;
  } else {
    int n = x.nrow();
    NumericVector rhoHat(n);
    IntegerVector ind=Range(0, n-1);
    for(int i=0; i<n; i++){
      IntegerVector ind1=clone(ind);
      ind1.erase(i);
      NumericMatrix x1 = ss(x, ind1);
      rhoHat[i]=Rcpp_gCor_Alpha(x1, alpha);
    }
    
    return VectorVar(rhoHat);
  }
  
}


// [[Rcpp::export]]
double Rcpp_HatV_gCor(NumericMatrix x) {
  
  int m=x.ncol();
  double mean, sum = 0;
  if(m==2){
    int n = x.nrow();
    NumericVector ind;
    NumericVector z= x.column(0);
    NumericVector y= x.column(1) ;
    NumericVector rhoHat(n);
    
    for(int i=0; i<n; i++){
      NumericVector x1 = clone(z);
      x1.erase(i);
      NumericVector y1=  clone(y);
      y1.erase(i);
      rhoHat[i]=rcpp_gCor(x1,y1);
    }
    mean = VectorSum(rhoHat)/n;
    for(int i=0; i<n; i++){
      sum += pow(rhoHat[i] - mean, 2.0);
    }
    return sum*(n-1)/n;
  } else {
    int n = x.nrow();
    NumericVector rhoHat(n);
    IntegerVector ind=Range(0, n-1);
    for(int i=0; i<n; i++){
      IntegerVector ind1=clone(ind);
      ind1.erase(i);
      NumericMatrix x1 = ss(x, ind1);
      rhoHat[i]=Rcpp_gCor(x1);
    }
    
    return VectorVar(rhoHat);
  }
  
}

// [[Rcpp::export]]
double Rcpp_HatV_KgCov(NumericMatrix x, double sigma) {
  
  int n = x.nrow();
  NumericVector rhoHat(n);
  IntegerVector ind=Range(0, n-1);
  for(int i=0; i<n; i++){
    IntegerVector ind1=clone(ind);
    ind1.erase(i);
    NumericMatrix x1 = ss(x, ind1);
    rhoHat[i]=Rcpp_KgCov(x1, sigma);
  }
  return VectorVar(rhoHat);
}

// [[Rcpp::export]]
double Rcpp_HatV_KgCov_Alpha(NumericMatrix x, double sigma, double alpha) {
  
  int n = x.nrow();
  NumericVector rhoHat(n);
  IntegerVector ind=Range(0, n-1);
  for(int i=0; i<n; i++){
    IntegerVector ind1=clone(ind);
    ind1.erase(i);
    NumericMatrix x1 = ss(x, ind1);
    rhoHat[i]=Rcpp_KgCov_Alpha(x1, sigma, alpha);
  }
  return VectorVar(rhoHat);
}


// [[Rcpp::export]]
double Rcpp_HatV_KgCor(NumericMatrix x, double sigma) {
  
  int n = x.nrow();
  NumericVector rhoHat(n);
  IntegerVector ind=Range(0, n-1);
  for(int i=0; i<n; i++){
    IntegerVector ind1=clone(ind);
    ind1.erase(i);
    NumericMatrix x1 = ss(x, ind1);
    rhoHat[i]=Rcpp_KgCor(x1, sigma);
  }
  return VectorVar(rhoHat);
  
}

// [[Rcpp::export]]
double Rcpp_HatV_KgCor_Alpha(NumericMatrix x, double sigma, double alpha) {
  
  int n = x.nrow();
  NumericVector rhoHat(n);
  IntegerVector ind=Range(0, n-1);
  for(int i=0; i<n; i++){
    IntegerVector ind1=clone(ind);
    ind1.erase(i);
    NumericMatrix x1 = ss(x, ind1);
    rhoHat[i]=Rcpp_KgCor_Alpha(x1, sigma, alpha);
  }
  return VectorVar(rhoHat);
}


