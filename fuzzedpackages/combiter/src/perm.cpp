#include <Rcpp.h>
#include <vector>

// [[Rcpp::export]]
std::vector<int> NextPerm(std::vector<int> x, int n) {
  // Compute the next permutation to x with maximum number at n
  //
  // Args:
  //   x   : integer vector
  //   n   : integer of the maximum x can have
  //
  // Returns:
  //   integer vector, same size as x
  //
  // Assumes:
  //   x has no duplicate values
  //   size of x is positive
  // Note:
  //   if x has reached the 'end', the value does not change

  // Strategy:
  // Starting from the end of the vector x,
  // find the first element that has a larger, unused by the elements to the left
  // and no greater than n.
  // replace this element by the smallest such number.
  // then, fill the elements to the right from the smallest unused values.
  // Example:
  //   x = 1 2 3 5, n = 5
  //   5 does not have larger unused.  3 has one. replace it by 4.
  //   then fill the last element by 3, which is the smallest
  //   result = 1 2 4 3
  //
  //   x = 4 5 3 2, n = 5
  //   2, 3, 5 does not have larger replacement.  4 is replaced by 5.
  //   then fill the right from the smallest
  //   result = 5 1 2 3


  int k = x.size();
  if (k > n) Rcpp::stop("k must be no greater than n");
  if (k <= 0 || n <= 0) Rcpp::stop("n and k must be positive");


  // is there next?
  bool end = true;
  for (int i = 0; i < k; i++)
  {
    if (x[i] != n-i) {
      end = false;
      break;
    }
  }
  if (end) return x;

  // used keeps which numbers are taken
  // careful with the difference in the indexing base
  std::vector<bool> used(n, false);
  // number of true counted, must equal k
  int used_count = 0;
  for (int i = 0; i < k; i++)
  {
    used[x[i]-1] = true;
    used_count++;
  }
  if (used_count != k) Rcpp::stop("x must not have duplicate");

  for (int i = k-1; i >= 0; i--)
  {
    int cur = x[i];
    used[cur-1] = false;
    for (int nex = cur+1; nex <= n; nex++)
    {
      if (!used[nex-1]) {
        x[i] = nex;
        used[cur-1] = false;
        used[nex-1] = true;

        // if this is the last element, done
        if (i >= k-1) return x;

        // fill the values to the right
        int m = i + 1;
        for (int j = 0; j < n; j++)
        {
          if (!used[j]) {
            x[m] = j+1;
            used[j] = true;
            m++;
            if (m >= k) return x;
          }
        }
        Rcpp::stop("could not fill vectors");
      }
    }
  }
  // no way to reach here, though
  Rcpp::stop("x has no next");
  return x;
}


// [[Rcpp::export]]
std::vector<int> PrevPerm(std::vector<int> x, int n) {
  // Compute the previous permutation to x with maximum number at n
  //
  // See NextPerm for details


  int k = x.size();
  if (k > n) Rcpp::stop("k must be no greater than n");
  if (k <= 0 || n <= 0) Rcpp::stop("n and k must be positive");


  // is there next?
  bool end = true;
  for (int i = 0; i < k; i++)
  {
    if (x[i] != i+1) {
      end = false;
      break;
    }
  }
  if (end) return x;

  // used keeps which numbers are taken
  // careful with the difference in the indexing base
  std::vector<bool> used(n, false);
  // number of true counted, must equal k
  int used_count = 0;
  for (int i = 0; i < k; i++)
  {
    used[x[i]-1] = true;
    used_count++;
  }
  if (used_count != k) Rcpp::stop("x must not have duplicate");

  for (int i = k-1; i >= 0; i--)
  {
    int cur = x[i];
    used[cur-1] = false;
    for (int nex = cur-1; nex >= 1; nex--)
    {
      if (!used[nex-1]) {
        x[i] = nex;
        used[cur-1] = false;
        used[nex-1] = true;

        // if this is the last element, done
        if (i >= k-1) return x;

        // fill the values to the right
        int m = i + 1;
        for (int j = n-1; j >= 0; j--)
        {
          if (!used[j]) {
            x[m] = j+1;
            used[j] = true;
            m++;
            if (m >= k) return x;
          }
        }
        Rcpp::stop("could not fill vectors");
      }
    }
  }
  // no way to reach here, though
  Rcpp::stop("x has no next");
  return x;
}





/*** R
x <- 1:3
n <- 5
ct <- 0
while (TRUE)
{
  ct <- ct + 1
  cat(sprintf("%3d : %s\n", ct,  paste0(x, collapse = " ")))
  y <- combiter:::NextPerm(x, n)
  if (all(x == y)) break
  x <- y
}

# NextPerm and PrevPerm are the inverse function for each other
x <- 1:3
n <- 5
ct <- 1
while (TRUE)
{
  cat(sprintf("%3d : %s\n", k,  paste0(x, collapse = " ")))
  y <- combiter:::NextPerm(x, n)
  if (all(x == y)) break
  z <- combiter:::PrevPerm(y, n)
  stopifnot(all(x == z))
  x <- y
  ct <- ct + 1
}
*/




// // [[Rcpp::export]]
// std::vector<int> NextPerm(std::vector<int> x) {
//   // Compute the next permutation to x
//   //
//   // Args:
//   //   x   : integer vector
//   //
//   // Returns:
//   //   integer vector, same size as x
//   //
//   // Assumes:
//   //   x has no duplicate values
//   //   size of x is positive
//   // Note:
//   //   if x has reached the 'end', the value does not change
//
//   // Quick input validation
//   if (x.size() == 0) Rcpp::stop("x must have positive size");
//
//   // Strategy:
//   // Starting from the end of the vector x,
//   // find the first decreasing adjacent values.
//   // If find one, then swap that samller value with the
//   // next larger value on the right, then
//   // reverse the numbers on the right.
//   // Example:
//   //   x = 1 2 3 4
//   //   3 is the first decreased value, swap that with 4.
//   //   On the right is just one element 3, so the reversing
//   //   does not affect.
//   //   result = 1 2 4 3
//   //
//   //   x = 2 4 3 1
//   //   2 is the first decreased value, which is swapped with 3.
//   //   On the right is 4 2 1, reversed to 1 2 4
//   //   result = 3 1 2 4
//   for (int i = x.size()-1; i > 0; i--)
//   {
//     if (x[i-1] < x[i]) {
//       for (int j = x.size()-1; j >= i; j--)
//       {
//         // Note: I could use binary search here
//         //       but the total runtime would be still O(n)
//         if (x[i-1] < x[j]) {
//           std::swap(x[i-1], x[j]);
//           std::reverse(x.begin()+i, x.end());
//           break;
//         }
//       }
//       break;
//     }
//   }
//   return x;
// }
//
//
//
// // [[Rcpp::export]]
// std::vector<int> PrevPerm(std::vector<int> x) {
//   // Compute the previous permutation to x
//   //
//   // Args:
//   //   x   : integer vector
//   //
//   // Returns:
//   //   integer vector, same size as x
//   //
//   // Assumes:
//   //   x has no duplicate values
//   //   size of x is positive
//   // Note:
//   //   If x has reached the 'end', then no change is made
//
//   // Quick input validation
//   if (x.size() == 0) Rcpp::stop("x must have positive size");
//
//   // Same algorithm as NextPerm, except for the direction of inequality
//   for (int i = x.size()-1; i > 0; i--)
//   {
//     if (x[i-1] > x[i]) {
//       for (int j = x.size()-1; j >= i; j--)
//       {
//         // Note: I could use binary search here
//         //       but the total runtime would be still O(n)
//         //       And practically n won't be so big
//         if (x[i-1] > x[j]) {
//           std::swap(x[i-1], x[j]);
//           std::reverse(x.begin()+i, x.end());
//           break;
//         }
//       }
//       break;
//     }
//   }
//   return x;
// }


// // [[Rcpp::export]]
// std::vector<int> PrevPerm(std::vector<int> x) {
//   // Compute the previous permutation to x
//   //
//   // Args:
//   //   x   : integer vector
//   //
//   // Returns:
//   //   integer vector, same size as x
//   //
//   // Assumes:
//   //   x has no duplicate values
//   //   size of x is positive
//   // Note:
//   //   If x has reached the 'end', then no change is made
//
//   // Quick input validation
//   if (x.size() == 0) Rcpp::stop("x must have positive size");
//
//   // Same algorithm as NextPerm, except for the direction of inequality
//   for (int i = x.size()-1; i > 0; i--)
//   {
//     if (x[i-1] > x[i]) {
//       for (int j = x.size()-1; j >= i; j--)
//       {
//         // Note: I could use binary search here
//         //       but the total runtime would be still O(n)
//         //       And practically n won't be so big
//         if (x[i-1] > x[j]) {
//           std::swap(x[i-1], x[j]);
//           std::reverse(x.begin()+i, x.end());
//           break;
//         }
//       }
//       break;
//     }
//   }
//   return x;
// }
