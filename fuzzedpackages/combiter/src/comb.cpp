#include <Rcpp.h>
#include <vector>
#include "comb.h"

// [[Rcpp::export]]
std::vector<int> NextComb(std::vector<int> x, int n)
{
  // Compute the next combination to x
  //
  // Args:
  //   x   : integer vector
  //   n   : max number
  //
  // Returns:
  //   integer vector, same size as x
  //
  // Assumes:
  //   x is sorted, positive and has no duplicate values
  //   size of x is positive and does not exceed n
  // Note:
  //   If x has reached the 'end', then no change is made

  int k = x.size();

  // A bit of input validation
  if (k == 0) Rcpp::stop("x must have positive size");
  if (k > n) Rcpp::stop("size of x must not exceed n");

  // If the first element is k-1, then there is no more next
  // just return x as is
  if (x[0] >= n-k+1) return x;

  // Strategy:
  // Increase the last element by one.
  // Notice that i-th entry can be at most (n-k+i+1), since
  // there are k-i-1 elements, larger than it; Otherwise,
  // some elements get larger than (n-k+i+1) + (k-i-1) = n, contradiction.
  // If the value at i-th entry exceeds that maximum, then
  // increase the (i-1)-th entry by one, then set the values of
  // i-th to the last entries the increments.
  // Example:
  //     2 5 7 8 with n = 8
  //  -> 2 5 7 9, now 9 is too large, so increase the next element
  //  -> 2 5 8 9, 8 is too large for the third entry, so increase the next
  //  -> 2 6 8 9. now set the last two as the smallest possible sequence
  //  -> 2 6 7 8
  //  result = 2 6 7 8
  //
  //     2 6 7 8 with n = 8
  //  -> 2 6 7 9
  //  -> 2 6 8 9
  //  -> 2 7 8 9
  //  -> 3 7 8 9
  //  -> 3 4 5 6
  //  result = 3 4 5 6

  x[k-1]++;
  for (int i = k - 1; i > 0; i--)
  {
    if (x[i] > n - k + i + 1) {
      x[i-1]++;
      for (int j = i; j < k; j++) x[j] = x[j-1] + 1;
    }
  }
  return x;
}


// [[Rcpp::export]]
std::vector<int> PrevComb(std::vector<int> x, int n)
{
  // Compute the previous combination to x
  //
  // Args:
  //   x   : integer vector
  //   n   : max number
  //
  // Returns:
  //   integer vector, same size as x
  //
  // Assumes:
  //   x is sorted, positive and has no duplicate values
  //   size of x is positive and does not exceed n

  // Note:
  //   If x has reached the 'end', then no change is made

  int k = x.size();

  // A bit of input validation
  if (k == 0) Rcpp::stop("x must have positive size");
  if (k > n) Rcpp::stop("size of x must not exceed n");

  // If the last element is k, then there is no more previous
  // just return x as is
  if (x[k-1] <= k) return x;

  // Strategy:
  // Decrease the last element by one.
  // For i-th entry, where i > 0, each value must be greater than
  // the (i-1)-th entry due to the sorted property.
  // If the i-th enetry becomes as smaller than or equal to (i-1)-th entry,
  // then decrease the (i-1)-th entry by one, and
  // set the i-th entry to be the largest possible value,
  // namely, (n-k+i+1)
  x[k-1]--;
  for (int i = k - 1; i > 0; i--)
  {
    if (x[i] <= x[i-1]) {
      x[i-1]--;
      x[i] = n - k + i + 1;
    } else {
      break;
    }
  }
  return x;
}






/*** R
x <- 1:4
n <- 8
ct <- 1
while (TRUE)
{
  cat(sprintf("%3d : %s\n", ct,  paste0(x, collapse = " ")))
  y <- combIter:::NextComb(x, n)
  if (all(x == y)) break
  x <- y
  ct <- ct + 1
}

# test that NextComb and PrevComb are inverse for each other
x <- 1:4
n <- 8
ct <- 1
while (TRUE)
{
  cat(sprintf("%3d : %s\n", ct,  paste0(x, collapse = " ")))
  y <- combIter:::NextComb(x, n)
  if (all(x == y)) break
  z <- combIter:::PrevComb(y, n)
  stopifnot(all(x == z))
  x <- y
  ct <- ct + 1
}
*/
