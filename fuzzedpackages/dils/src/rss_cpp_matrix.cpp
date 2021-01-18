#include "rss_cpp_matrix.h"

Rcpp::NumericMatrix PrepMatrix(Rcpp::NumericMatrix x) {
  // Given a matrix divides each row by the rowsum
  using namespace Rcpp;
  
  for(int i = 0; i < x.ncol(); i++) {
    x(i,i) = 0.0;
    x(i,_) = x(i,_) / sum(x(i,_));
  }
  return(x);
}

double RssThisRadius(Rcpp::NumericMatrix x, int v1, int v2, int r) {
  // Given prepped adjacency matrix, from node and to node, and radius
  // returns RSS for that one link
  using namespace Rcpp;
  double out = 0.0;
  
  int n = x.nrow();
  
  if (v1 == v2) {
    out = 0.0;
  } else if (1 == r) {
    out = x(v1, v2);
  } else if (2 == r) {
    out = sum(x(v1,_) * x(_,v2));
  } else if (3 == r) {
    NumericVector y(n);
    for(int i = 0; i < n; i++) {
      y[i] = RssThisRadius(x, v1, i, 2) - x(v1, v2) * x(v2, i);
    }
    out = sum(x(_,v2) * y) + x(v1, v2) * x(v2, v1) * x(v1, v2);
  } else if (4 == r) {
    NumericVector y(n);
    for(int i = 0; i < n; i++) {
      y[i] = RssThisRadius(x, v1, i, 3) - 
        x(v2, i) * sum(x(v1,_) * x(_,v2)) +
        x(v2, i) * x(v1, i) * x(i, v2) +
        x(v1, v2) * x(v2, v1) * x(v1, i) -
        x(v1, v2) * sum(x(v2,_) * x(_,i));
    }
    out = sum(x(_,v2) * y) + 
      x(v1,v2) * x(v2,v1) * sum(x(v1,_) * x(_,v2)) +
      x(v1,v2) * x(v1,v2) * sum(x(v2,_) * x(_,v1));
  } else {
    out = -1000.0;
  }
  
  return(out);
}

double RssCell(Rcpp::NumericMatrix x, int v1, int v2, int r) {
  using namespace Rcpp;
  double out = 0.0;
  
  for(int i = 0; i < r; i++) {
    out += RssThisRadius(x, v1, v2, i+1);
  }
  
  return(out);
}

SEXP rss_cell(SEXP xadj, SEXP vin, SEXP vout, SEXP radius, SEXP directed) {
  using namespace Rcpp;
  NumericMatrix x( xadj );
  int v1 = as<int>( vin ) - 1;
  int v2 = as<int>( vout ) - 1;
  int r = as<int>( radius );
  int d = as<int>( directed );
  
  if(1 == d) {
    return( wrap(RssCell(x, v1, v2, r)) );
  } else {
    return( wrap((RssCell(x, v1, v2, r) + RssCell(x, v2, v1, r)) / 2.0) );
  }
}

SEXP rss_cpp_matrix(SEXP xadj, SEXP radius, SEXP directed) {
  using namespace Rcpp;
    
  NumericMatrix x( xadj );
  int r = as<int>( radius );
  int d = as<int>( directed );
  int n = x.nrow();
  NumericMatrix out(x.nrow(), x.ncol());
  int remaining_seconds = 0;
  time_t begin_all, end_all;
  
  time(&begin_all);
  
  x = PrepMatrix(x);
  
  time_t begin_main, end_main;
  time(&begin_main);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      out(i,j) = RssCell(x, i, j, r);
    }
    time(&end_main);
    remaining_seconds = difftime(end_main, begin_main) * double(n-i) / (double(i) + 0.1);
    Rprintf("\rCompleted calculation for node %d of %d. Remaining time: %d seconds                    ", 
      i+1, n, remaining_seconds);
  }
  
  if(0 == d) {
    for(int i = 0; i < n; i++) {
      for(int j = i; j < n; j++) {
        out(i,j) = (out(i,j) + out(j,i)) / 2.0;
        out(j,i) = out(i,j);
      }
    }
  }
  
  time(&end_all);
  Rprintf("\nRSS Calculation complete. Run time: %d seconds\n", int(difftime(end_all, begin_all)));
  
  return wrap(out);
}
