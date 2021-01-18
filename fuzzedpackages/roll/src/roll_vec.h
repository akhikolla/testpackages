#ifndef ROLL_VEC_H
#define ROLL_VEC_H

#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

arma::ivec stl_sort_index(arma::vec& x);

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollAnyOnlineVec {
  
  const RVector<int> x;         // source
  const int n_rows_x;
  const int width;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_any;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAnyOnlineVec(const IntegerVector x, const int n_rows_x,
                   const int width, const int min_obs,
                   const bool na_restore, IntegerVector rcpp_any)
    : x(x), n_rows_x(n_rows_x),
      width(width), min_obs(min_obs),
      na_restore(na_restore), rcpp_any(rcpp_any) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int count = 0;
    int n_obs = 0;
    int x_new = 0;
    int x_old = 0;
    int sum_x = 0;
    
    for (int i = 0; i < n_rows_x; i++) {
      
      if ((x[i] == NA_INTEGER) || (x[i] == 0)) {
        
        x_new = 0;
        
      } else {
        
        x_new = 1;
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (x[i] != NA_INTEGER) {
          n_obs += 1;
        }
        
        sum_x = sum_x + x_new;
        
        count += 1;
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if ((x[i] != NA_INTEGER) && (x[i - width] == NA_INTEGER)) {
          
          n_obs += 1;
          
        } else if ((x[i] == NA_INTEGER) && (x[i - width] != NA_INTEGER)) {
          
          n_obs -= 1;
          
        }
        
        if ((x[i - width] == NA_INTEGER) || (x[i - width] == 0)) {
          
          x_old = 0;
          
        } else {
          
          x_old = 1;
          
        }
        
        sum_x = sum_x + x_new - x_old;
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && (x[i] != NA_INTEGER))) {
        
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_any[i] = 1;
          } else if (n_obs == count) {
            rcpp_any[i] = 0;
          } else {
            rcpp_any[i] = NA_INTEGER;
          }
          
        } else {
          rcpp_any[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_any[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollAnyOfflineVec : public Worker {
  
  const RVector<int> x;         // source
  const int n_rows_x;
  const int width;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_any;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAnyOfflineVec(const IntegerVector x, const int n_rows_x,
                    const int width, const int min_obs,
                    const bool na_restore, IntegerVector rcpp_any)
    : x(x), n_rows_x(n_rows_x),
      width(width), min_obs(min_obs),
      na_restore(na_restore), rcpp_any(rcpp_any) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      int sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && (x[i] != NA_INTEGER))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (x[i - count] != NA_INTEGER) {
            
            // compute the sum
            if (x[i - count] == 1) {
              sum_x = 1;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_any[i] = 1;
          } else if (n_obs == count) {
            rcpp_any[i] = 0;
          } else {
            rcpp_any[i] = NA_INTEGER;
          }
          
        } else {
          rcpp_any[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_any[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollAllOnlineVec {
  
  const RVector<int> x;         // source
  const int n_rows_x;
  const int width;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_all;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAllOnlineVec(const IntegerVector x, const int n_rows_x,
                   const int width, const int min_obs,
                   const bool na_restore, IntegerVector rcpp_all)
    : x(x), n_rows_x(n_rows_x),
      width(width), min_obs(min_obs),
      na_restore(na_restore), rcpp_all(rcpp_all) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int count = 0;
    int n_obs = 0;
    int x_new = 0;
    int x_old = 0;
    int sum_x = 0;
    
    for (int i = 0; i < n_rows_x; i++) {
      
      if ((x[i] == NA_INTEGER) || (x[i] != 0)) {
        
        x_new = 0;
        
      } else {
        
        x_new = 1;
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if ((x[i] != NA_INTEGER)) {
          n_obs += 1;
        }
        
        sum_x = sum_x + x_new;
        
        count += 1;
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if ((x[i] != NA_INTEGER) && (x[i - width] == NA_INTEGER)) {
          
          n_obs += 1;
          
        } else if ((x[i] == NA_INTEGER) && (x[i - width] != NA_INTEGER)) {
          
          n_obs -= 1;
          
        }
        
        if ((x[i - width] == NA_INTEGER) || (x[i - width] != 0)) {
          
          x_old = 0;
          
        } else {
          
          x_old = 1;
          
        }
        
        sum_x = sum_x + x_new - x_old;
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && (x[i] != NA_INTEGER))) {
        
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_all[i] = 0;
          } else if (n_obs == count) {
            rcpp_all[i] = 1;
          } else {
            rcpp_all[i] = NA_INTEGER;
          }
          
        } else {
          rcpp_all[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_all[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollAllOfflineVec : public Worker {
  
  const RVector<int> x;         // source
  const int n_rows_x;
  const int width;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_all;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAllOfflineVec(const IntegerVector x, const int n_rows_x,
                    const int width, const int min_obs,
                    const bool na_restore, IntegerVector rcpp_all)
    : x(x), n_rows_x(n_rows_x),
      width(width), min_obs(min_obs),
      na_restore(na_restore), rcpp_all(rcpp_all) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      int sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && (x[i] != NA_INTEGER))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (x[i - count] != NA_INTEGER) {
            
            // compute the sum
            if (x[i - count] == 0) {
              sum_x = 1;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_all[i] = 0;
          } else if (n_obs == count) {
            rcpp_all[i] = 1;
          } else {
            rcpp_all[i] = NA_INTEGER;
          }
          
        } else {
          rcpp_all[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_all[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollSumOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_sum;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSumOnlineVec(const NumericVector x, const int n,
                   const int n_rows_x, const int width,
                   const arma::vec arma_weights, const int min_obs,
                   const bool na_restore, arma::vec& arma_sum)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), arma_sum(arma_sum) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0;
    long double x_new = 0;
    long double x_old = 0;
    long double sum_x = 0;
    
    if (arma_weights[n - 1] == 0) {
      lambda = 1;
    } else if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_x; i++) {
      
      if (std::isnan(x[i])) {
        
        w_new = 0;
        x_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (width > 1) {
          sum_x = lambda * sum_x + w_new * x_new;
        } else {
          sum_x = w_new * x_new;
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width])) {
          
          w_old = 0;
          x_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          
        }
        
        if (width > 1) {
          sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        } else {
          sum_x = w_new * x_new;
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          arma_sum[i] = sum_x;
        } else {
          arma_sum[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sum[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollSumOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_sum;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSumOfflineVec(const NumericVector x, const int n,
                    const int n_rows_x, const int width,
                    const arma::vec arma_weights, const int min_obs,
                    const bool na_restore, arma::vec& arma_sum)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), arma_sum(arma_sum) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      long double sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // compute the sum
            sum_x += arma_weights[n - count - 1] * x[i - count];
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if (n_obs >= min_obs) {
          arma_sum[i] = sum_x;
        } else {
          arma_sum[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sum[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollProdOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_prod;         // destination (pass by reference)
  
  // initialize with source and destination
  RollProdOnlineVec(const NumericVector x, const int n,
                    const int n_rows_x, const int width,
                    const arma::vec arma_weights, const int min_obs,
                    const bool na_restore, arma::vec& arma_prod)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), arma_prod(arma_prod) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    int n_zero = 0;
    long double lambda = 0;
    long double n_new = 0;
    long double n_old = 0;
    long double n_exp = 0;
    long double w_new = 0;
    long double w_old = 0;
    long double x_new = 0;
    long double x_old = 0;
    long double prod_w = 1;
    long double prod_x = 1;
    
    if (arma_weights[n - 1] == 0) {
      lambda = 1;
    } else if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_x; i++) {
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (std::isnan(x[i])) {
          
          n_new = n_obs;
          w_new = 1;
          x_new = 1;
          
        } else {
          
          n_new = n_obs - 1;
          w_new = arma_weights[n - 1];
          
          if (x[i] == 0) {
            
            x_new = 1;
            n_zero += 1;
            
          } else {
            x_new = x[i];
          }
          
        }
        
        if (n_new == 0) {
          n_exp = 1;
        } else if (n_new > n_old) {
          n_exp = n_exp * lambda;
        } else if (n_new < n_old) {
          n_exp = n_exp / lambda;
        }
        
        n_old = n_new;
        prod_w *= w_new * n_exp;
        prod_x *= x_new;
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i])) {
          
          n_new = n_obs;
          w_new = 1;
          x_new = 1;
          
        } else {
          
          n_new = n_obs - 1;
          w_new = arma_weights[n - 1];
          
          if (x[i] == 0) {
            
            x_new = 1;
            n_zero += 1;
            
          } else {
            x_new = x[i];
          }
          
        }
        
        if (std::isnan(x[i - width])) {
          
          w_old = 1;
          x_old = 1;
          
        } else {
          
          if (x[i - width] == 0) {
            
            x_old = 1;
            n_zero -= 1;
            
          } else {
            x_old = x[i - width];
          }
          
          if (arma_weights[n - width] == 0) {
            w_old = 1;
          } else {
            w_old = arma_weights[n - width];
          }
          
        }
        
        if (n_new == 0) {
          n_exp = 1;
        } else if (n_new > n_old) {
          n_exp = n_exp * lambda;
        } else if (n_new < n_old) {
          n_exp = n_exp / lambda;
        }
        
        n_old = n_new;
        prod_w *= w_new * n_exp / w_old;
        prod_x *= x_new / x_old;
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          
          if (n_zero == 0) {
            arma_prod[i] = prod_w * prod_x;
          } else {
            arma_prod[i]= 0;
          }
          
        } else {
          arma_prod[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_prod[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollProdOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_prod;         // destination (pass by reference)
  
  // initialize with source and destination
  RollProdOfflineVec(const NumericVector x, const int n,
                     const int n_rows_x, const int width,
                     const arma::vec arma_weights, const int min_obs,
                     const bool na_restore, arma::vec& arma_prod)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), arma_prod(arma_prod) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      long double prod_x = 1;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // compute the product
            prod_x *= arma_weights[n - count - 1] * x[i - count];
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if (n_obs >= min_obs) {
          arma_prod[i] = prod_x;
        } else {
          arma_prod[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_prod[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollMeanOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_mean;         // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanOnlineVec(const NumericVector x, const int n,
                    const int n_rows_x, const int width,
                    const arma::vec arma_weights, const int min_obs,
                    const bool na_restore, arma::vec& arma_mean)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), arma_mean(arma_mean) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0;
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_x; i++) {
      
      if (std::isnan(x[i])) {
        
        w_new = 0;
        x_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (width > 1) {
          
          sum_w = lambda * sum_w + w_new;
          sum_x = lambda * sum_x + w_new * x_new;
          
        } else {
          
          sum_w = w_new;
          sum_x = w_new * x_new;
          
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width])) {
          
          w_old = 0;
          x_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          
        }
        
        if (width > 1) {
          
          sum_w = lambda * sum_w + w_new - lambda * w_old;
          sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
          
        } else {
          
          sum_w = w_new;
          sum_x = w_new * x_new;
          
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          arma_mean[i] = sum_x / sum_w;
        } else {
          arma_mean[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_mean[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollMeanOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_mean;         // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanOfflineVec(const NumericVector x, const int n,
                     const int n_rows_x, const int width,
                     const arma::vec arma_weights, const int min_obs,
                     const bool na_restore, arma::vec& arma_mean)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), arma_mean(arma_mean) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      long double sum_w = 0;
      long double sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            sum_w += arma_weights[n - count - 1];
            sum_x += arma_weights[n - count - 1] * x[i - count];
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if (n_obs >= min_obs) {
          arma_mean[i] = sum_x / sum_w;
        } else {
          arma_mean[i] = NA_REAL;
        }
        
        
      } else {
        
        // can be either NA or NaN
        arma_mean[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollMinOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<double> rcpp_min;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMinOnlineVec(const NumericVector x, const int n,
                   const int n_rows_x, const int width,
                   const arma::vec arma_weights, const int min_obs,
                   const bool na_restore, NumericVector rcpp_min)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_min(rcpp_min) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    int idxmin_x = 0;
    std::deque<int> deck(width);
    
    for (int i = 0; i < n_rows_x; i++) {
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] < x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        if (width > 1) {
          idxmin_x = deck.front();
        } else {
          idxmin_x = i;
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] < x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
          deck.pop_front();
        }
        
        if (width > 1) {
          idxmin_x = deck.front();
        } else {
          idxmin_x = i;
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          rcpp_min[i] = x[idxmin_x];
        } else {
          rcpp_min[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_min[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollMinOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<double> rcpp_min;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMinOfflineVec(const NumericVector x, const int n,
                    const int n_rows_x, const int width,
                    const arma::vec arma_weights, const int min_obs,
                    const bool na_restore, NumericVector rcpp_min)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_min(rcpp_min) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      int idxmin_x = i;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // last element of sorted array
            // note: 'weights' must be greater than 0
            if (std::isnan(x[idxmin_x]) || (x[i - count] <= x[idxmin_x])) {
              idxmin_x = i - count;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs >= min_obs)) {
          rcpp_min[i] = x[idxmin_x];
        } else {
          rcpp_min[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_min[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollMaxOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<double> rcpp_max;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxOnlineVec(const NumericVector x, const int n,
                   const int n_rows_x, const int width,
                   const arma::vec arma_weights, const int min_obs,
                   const bool na_restore, NumericVector rcpp_max)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_max(rcpp_max) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    int idxmax_x = 0;
    std::deque<int> deck(width);
    
    for (int i = 0; i < n_rows_x; i++) {
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] > x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        if (width > 1) {
          idxmax_x = deck.front();
        } else {
          idxmax_x = i;
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] > x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
          deck.pop_front();
        }
        
        if (width > 1) {
          idxmax_x = deck.front();
        } else {
          idxmax_x = i;
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          rcpp_max[i] = x[idxmax_x];
        } else {
          rcpp_max[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_max[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollMaxOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<double> rcpp_max;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxOfflineVec(const NumericVector x, const int n,
                    const int n_rows_x, const int width,
                    const arma::vec arma_weights, const int min_obs,
                    const bool na_restore, NumericVector rcpp_max)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_max(rcpp_max) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      int idxmax_x = i;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // first element of sorted array
            // note: 'weights' must be greater than 0
            if (std::isnan(x[idxmax_x]) || (x[i - count] >= x[idxmax_x])) {
              idxmax_x = i - count;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs >= min_obs)) {
          rcpp_max[i] = x[idxmax_x];
        } else {
          rcpp_max[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_max[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollIdxMinOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_idxmin;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMinOnlineVec(const NumericVector x, const int n,
                      const int n_rows_x, const int width,
                      const arma::vec arma_weights, const int min_obs,
                      const bool na_restore, IntegerVector rcpp_idxmin)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_idxmin(rcpp_idxmin) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    int idxmin_x = 0;
    std::deque<int> deck(width);
    
    for (int i = 0; i < n_rows_x; i++) {
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] < x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        if (width > 1) {
          idxmin_x = deck.front() + 1;
        } else {
          idxmin_x = 1;
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] < x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
          deck.pop_front();
        }
        
        if (width > 1) {
          idxmin_x = width - (i - deck.front());
        } else {
          idxmin_x = 1;
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          rcpp_idxmin[i] = idxmin_x;
        } else {
          rcpp_idxmin[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmin[i] = (int)x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollIdxMinOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_idxmin;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMinOfflineVec(const NumericVector x, const int n,
                       const int n_rows_x, const int width,
                       const arma::vec arma_weights, const int min_obs,
                       const bool na_restore, IntegerVector rcpp_idxmin)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_idxmin(rcpp_idxmin) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      int idxmin_x = i;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // last element of sorted array
            // note: 'weights' must be greater than 0
            if (std::isnan(x[idxmin_x]) || (x[i - count] <= x[idxmin_x])) {
              idxmin_x = i - count;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs >= min_obs)) {
          
          if (i < width) {
            rcpp_idxmin[i] = idxmin_x + 1;
          } else if (i >= width) {
            rcpp_idxmin[i] = width - (i - idxmin_x);
          }
          
        } else {
          rcpp_idxmin[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmin[i] = (int)x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollIdxMaxOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_idxmax;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMaxOnlineVec(const NumericVector x, const int n,
                      const int n_rows_x, const int width,
                      const arma::vec arma_weights, const int min_obs,
                      const bool na_restore, IntegerVector rcpp_idxmax)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_idxmax(rcpp_idxmax) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    int idxmax_x = 0;
    std::deque<int> deck(width);
    
    for (int i = 0; i < n_rows_x; i++) {
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] > x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        if (width > 1) {
          idxmax_x = deck.front() + 1;
        } else {
          idxmax_x = 1;
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (!std::isnan(x[i])) {
          
          while (!deck.empty() && (std::isnan(x[deck.back()]) || (x[i] > x[deck.back()]))) {
            deck.pop_back();
          }
          
          deck.push_back(i);
          
        }
        
        while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
          deck.pop_front();
        }
        
        if (width > 1) {
          idxmax_x = width - (i - deck.front());
        } else {
          idxmax_x = 1;
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          rcpp_idxmax[i] = idxmax_x;
        } else {
          rcpp_idxmax[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmax[i] = (int)x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollIdxMaxOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  RVector<int> rcpp_idxmax;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMaxOfflineVec(const NumericVector x, const int n,
                       const int n_rows_x, const int width,
                       const arma::vec arma_weights, const int min_obs,
                       const bool na_restore, IntegerVector rcpp_idxmax)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      na_restore(na_restore), rcpp_idxmax(rcpp_idxmax) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      int idxmax_x = i;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // first element of sorted array
            // note: 'weights' must be greater than 0
            if (std::isnan(x[idxmax_x]) || (x[i - count] >= x[idxmax_x])) {
              idxmax_x = i - count;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs >= min_obs)) {
          
          if (i < width) {
            rcpp_idxmax[i] = idxmax_x + 1;
          } else if (i >= width) {
            rcpp_idxmax[i] = width - (i - idxmax_x);
          }
          
        } else {
          rcpp_idxmax[i] = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmax[i] = (int)x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollQuantileOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const double p;
  const int min_obs;
  const bool na_restore;
  RVector<double> rcpp_quantile;// destination (pass by reference)
  
  // initialize with source and destination
  RollQuantileOfflineVec(const NumericVector x, const int n,
                         const int n_rows_x, const int width,
                         const arma::vec arma_weights, const double p, 
                         const int min_obs, const bool na_restore,
                         NumericVector rcpp_quantile)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), p(p),
      min_obs(min_obs), na_restore(na_restore),
      rcpp_quantile(rcpp_quantile) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        int k = 0;
        int count = 0;
        long double sum_w = 0;
        
        int offset = std::max(0, i - width + 1);
        int n_size_x = i - offset + 1;
        arma::vec x_subset(n_size_x);
        arma::vec arma_weights_subset(n_size_x);
        
        std::copy(x.begin() + offset, x.begin() + i + 1,
                  x_subset.begin());
        std::copy(arma_weights.begin() + n - n_size_x, arma_weights.begin() + n,
                  arma_weights_subset.begin());
        
        // similar to R's sort with 'index.return = TRUE'
        arma::ivec sort_ix = stl_sort_index(x_subset);
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (n_size_x - 1 >= count)) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value
          if (!std::isnan(x_subset[k])) {
            
            // compute the sum
            sum_w += arma_weights_subset[k];
            
          }
          
          count += 1;
          
        }
        
        count = 0;
        int n_obs = 0;
        int idxquantile_x = 0;
        bool status = false;
        long double sum_upper_w = 0;
        long double sum_upper_w_temp = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (n_size_x - 1 >= count)) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value
          if (!std::isnan(x_subset[k])) {
            
            sum_upper_w += arma_weights_subset[k];
            
            // last element of sorted array that is 'p' of 'weights'
            // note: 'weights' must be greater than 0
            if (!status && (sum_upper_w / sum_w >= p)) {
              
              status = true;
              idxquantile_x = n_size_x - count - 1;
              sum_upper_w_temp = sum_upper_w;
              
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs >= min_obs)) {
          
          k = sort_ix[idxquantile_x];
          
          // average if upper and lower weight is equal
          if (std::abs(sum_upper_w_temp / sum_w - p) <= sqrt(arma::datum::eps)) {
            
            int k_lower = sort_ix[idxquantile_x - 1];
            rcpp_quantile[i] = (x_subset[k] + x_subset[k_lower]) / 2;
            
          } else {
            rcpp_quantile[i] = x_subset[k];
          }
          
        } else {
          rcpp_quantile[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_quantile[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollVarOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_var;          // destination (pass by reference)
  
  // initialize with source and destination
  RollVarOnlineVec(const NumericVector x, const int n,
                   const int n_rows_x, const int width,
                   const arma::vec arma_weights, const bool center,
                   const int min_obs, const bool na_restore,
                   arma::vec& arma_var)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), center(center),
      min_obs(min_obs), na_restore(na_restore),
      arma_var(arma_var) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0; 
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_x; i++) {
      
      if (std::isnan(x[i])) {
        
        w_new = 0;
        x_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && (n_obs > 1)) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i])) {
          
          sumsq_x = lambda * sumsq_x;
          
        } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
          
          sumsq_x = w_new * pow(x_new, (long double)2.0);
          
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width])) {
          
          w_old = 0;
          x_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          
        }
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x;
          
        }
        
      }
      
      // don't compute if missing value
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          arma_var[i] = sumsq_x / (sum_w - sumsq_w / sum_w);
        } else {
          arma_var[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_var[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollVarOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_var;          // destination (pass by reference)
  
  // initialize with source and destination
  RollVarOfflineVec(const NumericVector x, const int n,
                    const int n_rows_x, const int width,
                    const arma::vec arma_weights, const bool center,
                    const int min_obs, const bool na_restore,
                    arma::vec& arma_var)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), center(center),
      min_obs(min_obs), na_restore(na_restore),
      arma_var(arma_var) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      long double mean_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x[i - count];
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        int count = 0;
        int n_obs = 0;
        long double sum_w = 0;
        long double sumsq_w = 0;
        long double sumsq_x = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += pow(arma_weights[n - count - 1], 2.0);
            
            // compute the sum of squares with 'center' argument
            if (center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count] - mean_x, (long double)2.0);
            } else if (!center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count], 2.0);
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          arma_var[i] = sumsq_x / (sum_w - sumsq_w / sum_w);
        } else {
          arma_var[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_var[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollSdOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_sd;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSdOnlineVec(const NumericVector x, const int n,
                  const int n_rows_x, const int width,
                  const arma::vec arma_weights, const bool center,
                  const int min_obs, const bool na_restore,
                  arma::vec& arma_sd)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), center(center),
      min_obs(min_obs), na_restore(na_restore),
      arma_sd(arma_sd) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0; 
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    long double var_x = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_x; i++) {
      
      if (std::isnan(x[i])) {
        
        w_new = 0;
        x_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && (n_obs > 1)) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i])) {
          
          sumsq_x = lambda * sumsq_x;
          
        } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
          
          sumsq_x = w_new * pow(x_new, (long double)2.0);
          
        }
        
        var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width])) {
          
          w_old = 0;
          x_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          
        }
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x;
          
        }
        
        var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
        
      }
      
      // don't compute if missing value
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if ((var_x < 0) || (sqrt(var_x) <= sqrt(arma::datum::eps))) {
            arma_sd[i] = 0;
          } else {
            arma_sd[i] = sqrt(var_x);
          }
          
        } else {
          arma_sd[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sd[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollSdOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_sd;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSdOfflineVec(const NumericVector x, const int n,
                   const int n_rows_x, const int width,
                   const arma::vec arma_weights, const bool center,
                   const int min_obs, const bool na_restore,
                   arma::vec& arma_sd)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), center(center),
      min_obs(min_obs), na_restore(na_restore),
      arma_sd(arma_sd) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      long double mean_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x[i - count];
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        int count = 0;
        int n_obs = 0;
        long double sum_w = 0;
        long double sumsq_w = 0;
        long double sumsq_x = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += pow(arma_weights[n - count - 1], 2.0);
            
            // compute the sum of squares with 'center' argument
            if (center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count] - mean_x, (long double)2.0);
            } else if (!center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count], 2.0);
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          arma_sd[i] = sqrt(sumsq_x / (sum_w - sumsq_w / sum_w));
        } else {
          arma_sd[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sd[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollScaleOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleOnlineVec(const NumericVector x, const int n,
                     const int n_rows_x, const int width,
                     const arma::vec arma_weights, const bool center,
                     const bool scale, const int min_obs,
                     const bool na_restore, arma::vec& arma_scale)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), center(center),
      scale(scale), min_obs(min_obs),
      na_restore(na_restore), arma_scale(arma_scale) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0; 
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    long double var_x = 0;
    long double x_ij = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_x; i++) {
      
      if (std::isnan(x[i])) {
        
        w_new = 0;
        x_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        x_ij = x[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        if (width > 1) {
          
          sum_w = lambda * sum_w + w_new;
          sum_x = lambda * sum_x + w_new * x_new;
          sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
          
        } else {
          
          sum_w = w_new;
          sum_x = w_new * x_new;
          sumsq_w = pow(w_new, (long double)2.0);
          
        }
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i])) {
            
            sumsq_x = lambda * sumsq_x;
            
          } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            
          }
          
          var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
          
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width])) {
          
          w_old = 0;
          x_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          
        }
        
        if (width > 1) {
          
          sum_w = lambda * sum_w + w_new - lambda * w_old;
          sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
          sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
            pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
          
        } else {
          
          sum_w = w_new;
          sum_x = w_new * x_new;
          sumsq_w = pow(w_new, (long double)2.0);
          
        }
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x;
            
          }
          
          var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
          
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (n_obs >= min_obs) {
          
          if (scale && ((n_obs <= 1) || (var_x < 0) ||
              (sqrt(var_x) <= sqrt(arma::datum::eps)))) {
            arma_scale[i] = NA_REAL;
          } else if (center && scale) {
            arma_scale[i] = (x_ij - mean_x) / sqrt(var_x);
          } else if (!center && scale) {
            arma_scale[i] = x_ij / sqrt(var_x);
          } else if (center && !scale) {
            arma_scale[i] = x_ij - mean_x;
          } else if (!center && !scale) {
            arma_scale[i] = x_ij;
          }
          
        } else {
          arma_scale[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_scale[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollScaleOfflineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleOfflineVec(const NumericVector x, const int n,
                      const int n_rows_x, const int width,
                      const arma::vec arma_weights, const bool center,
                      const bool scale, const int min_obs,
                      const bool na_restore, arma::vec& arma_scale)
    : x(x), n(n),
      n_rows_x(n_rows_x), width(width),
      arma_weights(arma_weights), center(center),
      scale(scale), min_obs(min_obs),
      na_restore(na_restore), arma_scale(arma_scale) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      long double mean_x = 0;
      long double var_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x[i - count];
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          int count = 0;
          int n_obs = 0;
          long double sum_w = 0;
          long double sumsq_w = 0;
          long double sumsq_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the sum of squares with 'center' argument
              if (center) {
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x[i - count] - mean_x, (long double)2.0);
              } else if (!center) {
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x[i - count], 2.0);
              }
              
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the unbiased estimate of variance
          var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
          
        }
        
        int count = 0;
        int n_obs = 0;
        bool any_na = false;
        long double x_ij = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // keep first non-missing value
            if (!any_na) {
              x_ij = x[i - count];
            }
            
            any_na = true;
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if (n_obs >= min_obs) {
          
          if (scale && ((n_obs <= 1) || (var_x < 0) ||
              (sqrt(var_x) <= sqrt(arma::datum::eps)))) {
            arma_scale[i] = NA_REAL;
          } else if (center && scale) {
            arma_scale[i] = (x_ij - mean_x) / sqrt(var_x);
          } else if (!center && scale) {
            arma_scale[i] = x_ij / sqrt(var_x);
          } else if (center && !scale) {
            arma_scale[i] = x_ij - mean_x;
          } else if (!center && !scale) {
            arma_scale[i] = x_ij;
          }
          
        } else {
          arma_scale[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_scale[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollCovOnlineVecXX {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOnlineVecXX(const NumericVector x, const int n,
                     const int n_rows_xy, const int width,
                     const arma::vec arma_weights, const bool center,
                     const bool scale, const int min_obs,
                     const bool na_restore, arma::vec& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), width(width),
      arma_weights(arma_weights), center(center),
      scale(scale), min_obs(min_obs),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0;      
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double sumsq_xy = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_xy; i++) {
      
      if (std::isnan(x[i])) {
        
        w_new = 0;
        x_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i])) {
          n_obs += 1;
        }
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i])) {
            
            sumsq_x = lambda * sumsq_x;
            
          } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && (n_obs > 1)) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i])) {
          
          sumsq_xy = lambda * sumsq_xy;
          
        } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
          
          sumsq_xy = w_new * pow(x_new, (long double)2.0);
          
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          n_obs += 1;
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width])) {
          
          w_old = 0;
          x_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          
        }
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x;
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy;
          
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (scale) {
            
            // don't compute if the standard deviation is zero
            if ((sumsq_x < 0) || (sqrt(sumsq_x) <= sqrt(arma::datum::eps)))  {
              arma_cov[i] = NA_REAL;
            } else {
              arma_cov[i] = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_x));
            }
            
          } else if (!scale) {
            arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
          }
          
        } else {
          arma_cov[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_cov[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollCovOnlineVecXY {
  
  const RVector<double> x;      // source
  const RVector<double> y;      // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOnlineVecXY(const NumericVector x, const NumericVector y,
                     const int n, const int n_rows_xy,
                     const int width, const arma::vec arma_weights,
                     const bool center, const bool scale,
                     const int min_obs, const bool na_restore,
                     arma::vec& arma_cov)
    : x(x), y(y),
      n(n), n_rows_xy(n_rows_xy),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), na_restore(na_restore),
      arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0;
    long double x_new = 0;
    long double x_old = 0;
    long double y_new = 0;
    long double y_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sum_y = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double sumsq_y = 0;
    long double sumsq_xy = 0;
    long double mean_prev_x = 0;
    long double mean_prev_y = 0;
    long double mean_x = 0;
    long double mean_y = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_xy; i++) {
      
      if (std::isnan(x[i]) || std::isnan(y[i])) {
        
        w_new = 0;
        x_new = 0;
        y_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        y_new = y[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && !std::isnan(y[i])) {
          n_obs += 1;
        }
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sum_y = lambda * sum_y + w_new * y_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_prev_y = mean_y;
          mean_x = sum_x / sum_w;
          mean_y = sum_y / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            sumsq_y = lambda * sumsq_y +
              w_new * (y_new - mean_y) * (y_new - mean_prev_y);
            
          } else if (std::isnan(x[i]) || std::isnan(y[i])) {
            
            sumsq_x = lambda * sumsq_x;
            sumsq_y = lambda * sumsq_y;
            
          } else if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            sumsq_y = w_new * pow(y_new, (long double)2.0);
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs > 1)) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (y_new - mean_prev_y);
          
        } else if (std::isnan(x[i]) || std::isnan(y[i])) {
          
          sumsq_xy = lambda * sumsq_xy;
          
        } else if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs == 1) && !center) {
          
          sumsq_xy = w_new * x_new * y_new;
          
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
            (std::isnan(x[i - width]) || std::isnan(y[i - width]))) {
          
          n_obs += 1;
          
        } else if ((std::isnan(x[i]) || std::isnan(y[i])) &&
          (!std::isnan(x[i - width]) && !std::isnan(y[i - width]))) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width]) || std::isnan(y[i - width])) {
          
          w_old = 0;
          x_old = 0;
          y_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          y_old = y[i - width];
          
        }
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sum_y = lambda * sum_y + w_new * y_new - lambda * w_old * y_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_prev_y = mean_y;
          mean_x = sum_x / sum_w;
          mean_y = sum_y / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
              !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
              
              sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            sumsq_y = lambda * sumsq_y +
              w_new * (y_new - mean_y) * (y_new - mean_prev_y) -
              lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
            
          } else if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
            (std::isnan(x[i - width]) || std::isnan(y[i - width]))) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            sumsq_y = lambda * sumsq_y +
              w_new * (y_new - mean_y) * (y_new - mean_prev_y);
            
          } else if ((std::isnan(x[i]) || std::isnan(y[i])) &&
            !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
            
            sumsq_x = lambda * sumsq_x -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            sumsq_y = lambda * sumsq_y -
              lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
            
          } else if (std::isnan(x[i]) || std::isnan(y[i]) ||
            std::isnan(x[i - width]) || std::isnan(y[i - width])) {
            
            sumsq_x = lambda * sumsq_x;
            sumsq_y = lambda * sumsq_y;
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
            !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
            
            sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (y_new - mean_prev_y) -
            lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
          
        } else if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
          (std::isnan(x[i - width]) || std::isnan(y[i - width]))) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (y_new - mean_prev_y);
          
        } else if ((std::isnan(x[i]) || std::isnan(y[i])) &&
          !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy -
          lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
          
        } else if (std::isnan(x[i]) || std::isnan(y[i]) ||
          std::isnan(x[i - width]) || std::isnan(y[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy;
          
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]) &&
          !std::isnan(y[i]))) {
          
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (scale) {
              
              // don't compute if the standard deviation is zero
              if ((sumsq_x < 0) || (sumsq_y < 0) ||
                  (sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                
                arma_cov[i] = NA_REAL;
                
              } else {
                arma_cov[i] = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
              }
              
            } else if (!scale) {
              arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
            }
            
          } else {
            arma_cov[i] = NA_REAL;
          }
          
      } else {
        
        // can be either NA or NaN
        if (std::isnan(x[i])) {
          arma_cov[i] = x[i];
        } else {
          arma_cov[i] = y[i];
        }
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollCovOfflineVecXX : public Worker {
  
  const RVector<double> x;       // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOfflineVecXX(const NumericVector x, const int n,
                      const int n_rows_xy, const int width,
                      const arma::vec arma_weights, const bool center,
                      const bool scale, const int min_obs,
                      const bool na_restore, arma::vec& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), width(width),
      arma_weights(arma_weights), center(center),
      scale(scale), min_obs(min_obs),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 3D array (lower triangle)
      int i = z;
      
      long double sumsq_x = 0;
      long double mean_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x[i - count];
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          int count = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the sum of squares with 'center' argument
              if (center) {
                
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x[i - count] - mean_x, (long double)2.0);
                
              } else if (!center) {
                
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x[i - count], 2.0);
                
              }
              
            }
            
            count += 1;
            
          }
          
        }
        
        int count = 0;
        int n_obs = 0;
        long double sum_w = 0;
        long double sumsq_w = 0;
        long double sumsq_xy = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += pow(arma_weights[n - count - 1], 2.0);
            
            // compute the sum of squares with 'center' argument
            if (center) {
              sumsq_xy += arma_weights[n - count - 1] *
                pow(x[i - count] - mean_x, (long double)2.0);
            } else if (!center) {
              sumsq_xy += arma_weights[n - count - 1] *
                pow(x[i - count], 2.0);
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (scale) {
            
            // don't compute if the standard deviation is zero
            if ((sumsq_x < 0) || (sqrt(sumsq_x) <= sqrt(arma::datum::eps))) {
              arma_cov[i] = NA_REAL;
            } else {
              arma_cov[i] = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_x));
            }
            
          } else if (!scale) {
            arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
          }
          
        } else {
          arma_cov[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_cov[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollCovOfflineVecXY : public Worker {
  
  const RVector<double> x;       // source
  const RVector<double> y;       // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOfflineVecXY(const NumericVector x, const NumericVector y,
                      const int n, const int n_rows_xy,
                      const int width, const arma::vec arma_weights,
                      const bool center, const bool scale,
                      const int min_obs, const bool na_restore,
                      arma::vec& arma_cov)
    : x(x), y(y),
      n(n), n_rows_xy(n_rows_xy),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), na_restore(na_restore),
      arma_cov(arma_cov) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 3D array
      int i = z;
      
      long double sumsq_x = 0;
      long double sumsq_y = 0;
      long double mean_x = 0;
      long double mean_y = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]) &&
          !std::isnan(y[i]))) {
          
          if (center) {
            
            int count = 0;
            long double sum_w = 0;
            long double sum_x = 0;
            long double sum_y = 0;
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value
              if (!std::isnan(x[i - count]) && !std::isnan(y[i - count])) {
                
                // compute the sum
                sum_w += arma_weights[n - count - 1];
                sum_x += arma_weights[n - count - 1] * x[i - count];
                sum_y += arma_weights[n - count - 1] * y[i - count];
                
              }
              
              count += 1;
              
            }
            
            // compute the mean
            mean_x = sum_x / sum_w;
            mean_y = sum_y / sum_w;
            
          }
          
          if (scale) {
            
            int count = 0;
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value
              if (!std::isnan(x[i - count]) && !std::isnan(y[i - count])) {
                
                // compute the sum of squares with 'center' argument
                if (center) {
                  
                  sumsq_x += arma_weights[n - count - 1] *
                    pow(x[i - count] - mean_x, (long double)2.0);
                  sumsq_y += arma_weights[n - count - 1] *
                    pow(y[i - count] - mean_y, (long double)2.0);
                  
                } else if (!center) {
                  
                  sumsq_x += arma_weights[n - count - 1] *
                    pow(x[i - count], 2.0);
                  sumsq_y += arma_weights[n - count - 1] *
                    pow(y[i - count], 2.0);
                  
                }
                
              }
              
              count += 1;
              
            }
            
          }
          
          int count = 0;
          int n_obs = 0;
          long double sum_w = 0;
          long double sumsq_w = 0;
          long double sumsq_xy = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count]) && !std::isnan(y[i - count])) {
              
              sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the sum of squares with 'center' argument
              if (center) {
                sumsq_xy += arma_weights[n - count - 1] *
                  (x[i - count] - mean_x) * (y[i - count] - mean_y);
              } else if (!center) {
                sumsq_xy += arma_weights[n - count - 1] *
                  x[i - count] * y[i - count];
              }
              
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (scale) {
              
              // don't compute if the standard deviation is zero
              if ((sumsq_x < 0) || (sumsq_y < 0) ||
                  (sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                
                arma_cov[i] = NA_REAL;
                
              } else {
                arma_cov[i] = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
              }
              
            } else if (!scale) {
              arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
            }
            
          } else {
            arma_cov[i] = NA_REAL;
          }
          
      } else {
        
        // can be either NA or NaN
        if (std::isnan(x[i])) {
          arma_cov[i] = x[i];
        } else {
          arma_cov[i] = y[i];
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using a standard algorithm
struct RollLmVecInterceptFALSE : public Worker {
  
  const arma::cube arma_cov;    // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_n_obs;
  const arma::vec arma_sum_w;
  arma::vec& arma_coef;         // destination (pass by reference)
  arma::vec& arma_rsq;
  arma::vec& arma_se;
  
  // initialize with source and destination
  RollLmVecInterceptFALSE(const arma::cube arma_cov, const int n,
                          const int n_rows_xy, const int width,
                          const arma::vec arma_n_obs, const arma::vec arma_sum_w,
                          arma::vec& arma_coef, arma::vec& arma_rsq,
                          arma::vec& arma_se)
    : arma_cov(arma_cov), n(n),
      n_rows_xy(n_rows_xy), width(width),
      arma_n_obs(arma_n_obs), arma_sum_w(arma_sum_w),
      arma_coef(arma_coef), arma_rsq(arma_rsq),
      arma_se(arma_se) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat sigma = arma_cov.slice(i);
      arma::mat A = sigma.submat(0, 0, 0, 0);
      arma::mat b = sigma.submat(0, 1, 0, 1);
      arma::vec coef(1);
      
      // check if missing value is present
      bool any_na = sigma.has_nan();
      
      // don't compute if missing value 
      if (!any_na) {
        
        // check if solution is found      
        bool status_solve = arma::solve(coef, A, b, arma::solve_opts::no_approx);
        int df_fit = 1;
        
        // don't find approximate solution for rank deficient system,
        // and the width and current row must be greater than the
        // number of variables
        if (status_solve && (arma_n_obs[i] >= df_fit)) {
          
          // coefficients
          arma::mat trans_coef = trans(coef);
          arma_coef[i] = as_scalar(trans_coef);
          
          // r-squared
          long double var_y = sigma(1, 1);
          if ((var_y < 0) || (sqrt(var_y) <= sqrt(arma::datum::eps))) {      
            arma_rsq[i] = NA_REAL;
          } else {
            arma_rsq[i] = as_scalar(trans_coef * A * coef) / var_y;
          }
          
          // check if matrix is singular
          arma::mat A_inv(1, 1);
          bool status_inv = arma::inv(A_inv, A);
          int df_resid = arma_n_obs[i] - 2 + 1;
          
          if (status_inv && (df_resid > 0)) {
            
            // standard errors
            long double var_resid = (1 - arma_rsq[i]) * var_y / df_resid;
            arma_se[i] = as_scalar(sqrt(var_resid * trans(diagvec(A_inv))));
            
          } else {
            
            arma_se[i] = NA_REAL;
            
          }
          
        } else {
          
          arma_coef[i] = NA_REAL;
          arma_rsq[i] = NA_REAL;
          arma_se[i] = NA_REAL;
          
        }
        
      } else {
        
        arma_coef[i] = NA_REAL;
        arma_rsq[i] = NA_REAL;
        arma_se[i] = NA_REAL;
        
      }
      
    }
  }
  
};

#endif
