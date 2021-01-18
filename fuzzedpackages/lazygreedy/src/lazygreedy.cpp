// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include "greedyspanner.h"
#include "structures.h"
#include "binaryheap.h"

using namespace std;
using namespace Rcpp;

//' Applying the Lazy-Greedy Algorithm
//'
//' Function \code{lazyGreedy} is an R wrapper for an efficient C++ implementation of the Lazy-Greedy spanning algorithm. The C++ implementation executes many thousands of times faster than a pure R implementation. Both the algorithm and (most of) the C++ implementation were developed by Quirijn W. Bouts, Alex P. ten Brink, and Kevin Buchin.
//'
//' @param V a numeric \eqn{n}-by-2 matrix the \eqn{i}th row of which contains the location in \eqn{R^2} of vertex \eqn{i}.
//' @param t the desired dilation, a positive real number.
//' @return Function \code{lazyGreedy} returns a greedy \eqn{t}-spanner for the set of vertices \code{V}. The result takes the form of an edge list.
//' @references
//' Bouts, Q. W., ten Brink, A. P., and Buchin, K. (2014). A framework for computing the greedy spanner. In \emph{30th ACM Symposium on Computational Geometry} (SoCG, Kyoto, Japan, June 8-11, 2014) (pp. 11-19). Association for Computing Machinery, Inc.
//' @examples
//' n = 20
//' V = cbind(runif(n), runif(n))
//' spanner = lazyGreedy(V, t = 2)
//'
//' \dontrun{
//' require("network")
//' G = network(spanner, directed = FALSE)
//' plot(G, coord = V, label = 1:n, jitter = FALSE)
//' }
//' @export
// [[Rcpp::export]]
arma::mat lazyGreedy(const arma::mat& V, double t)  // G = (V, E)
{
  int n = V.n_rows;
  pointset points(n);
  for (int i = 0; i < n; i++)
  {
    points[i].x = V(i, 0);
    points[i].y = V(i, 1);
  }
  edgelist E;
  int edgecount = GreedyLinspace3<BinHeap<double> >(points, t, E);
  arma::mat edges(edgecount, 2);
  for (int i = 0; i < edgecount; i++)
  {
    edges(i, 0) = E[i].x + 1;
    edges(i, 1) = E[i].y + 1;
  }
  return edges;
}
