#include <Rcpp.h>
#include <vector>
#include "comb.h"


// [[Rcpp::export]]
std::vector<int> NextSubset(std::vector<int> x, int n)
{
  // Compute next subset vector to x
  //
  // Args:
  //   x   : integer vector
  //   n   : max number
  //
  // Returns:
  //   integer vector
  //
  // Assumes:
  //   x is sorted, positive and has no duplicate values
  //   size of x does not exceed n
  // Note:
  //   If x has reached the 'end', then no change is made

  // Strategy:
  // Let k be the size of x.
  // If k = 0, i.e. empty set, then return 1.
  // If x has reached the 'combination end' for the size k, then
  // there are two cases:
  //   k = n --> there is no more subset, return x as is
  //   k < n --> return (1, ..., k+1)
  // Otherwise, return NextComb

  int k = x.size();

  // Quick validation
  if (k > n) Rcpp::stop("size of x must not exceed n");

  if (k == 0) {
    std::vector<int> out(1, 1);
    return out;
  }

  if (x[0] >= n-k+1) {
    if (k >= n) return x;

    std::vector<int> out(k+1);
    for (int i = 0; i < k+1; i++) out[i] = i+1;
    return out;
  }

  return NextComb(x, n);
}


// [[Rcpp::export]]
std::vector<int> PrevSubset(std::vector<int> x, int n)
{
  // Compute previous subset vector to x
  //
  // Args:
  //   x   : integer vector
  //   n   : max number
  //
  // Returns:
  //   integer vector
  //
  // Assumes:
  //   x is sorted, positive and has no duplicate values
  //   size of x does not exceed n
  // Note:
  //   If x has reached the 'end', then no change is made

  // Strategy:
  // Let k be the size of x.
  // If k = 0, then there is no more previous, return x as is
  // If x has reached the 'combination end' for the size k, then
  // return the largest possible sequence of size k-1,
  // namely, (n-k+2, ..., n)
  // Otherwise, return PrevComb

  int k = x.size();

  // Quick validation
  if (k > n) Rcpp::stop("size of x must not exceed n");

  if (k == 0) return x;

  if (x[k-1] <= k) {
    std::vector<int> out(k-1);
    for (int i = 0; i < k-1; i++) out[i] = n-k+2+i;
    return out;
  }

  return PrevComb(x, n);
}



/*** R
x <- integer(0)
n <- 6
ct <- 1
while (TRUE)
{
  cat(sprintf("%3d : %s\n", ct,  paste0(x, collapse = " ")))
  y <- combIter:::NextSubset(x, n)
  if (identical(x, y)) break
  x <- y
  ct <- ct + 1
}

# check that NextSubset and PrevSubset are inverse for each other
x <- integer(0)
n <- 6
ct <- 1
while (TRUE)
{
  cat(sprintf("%3d : %s\n", ct,  paste0(x, collapse = " ")))
  y <- combIter:::NextSubset(x, n)
  if (identical(x, y)) break
  z <- combIter:::PrevSubset(y, n)
  stopifnot(identical(x, z))
  x <- y
  ct <- ct + 1
}
*/
