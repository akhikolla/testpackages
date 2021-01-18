#include <vector>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> NextCartes(std::vector<int> x, std::vector<int> nvec)
{
  // returns next cartesian product
  //
  // args:
  //   x:    int vector of current state
  //   nvec: int vector of max of each element, must have same size as x
  //
  // return:
  //   int vector of same size as x
  //   if x is already the end state, return x as is
  //
  // note:
  //   this function does not check if the current state is legal

  if (x.size() != nvec.size()) stop("x and nvec must have a same size");

  // check if this is the last value
  // i.e. x == nvec for all
  bool is_end = true;
  for (size_t i=0; i < x.size(); i++)
  {
    if (x[i] < nvec[i]) {
      is_end = false;
      break;
    }
  }
  if (is_end) return x;

  // increment the last value and iteratively change other digits
  int k = x.size();
  x[k-1]++;
  for (int i=k-1; i>=1; i--)
  {
    if (x[i] > nvec[i]) {
      x[i-1]++;
      x[i] = 1;
    }
  }
  return x;
}


// [[Rcpp::export]]
std::vector<int> PrevCartes(std::vector<int> x, std::vector<int> nvec)
{
  // returns previous cartesian product
  //
  // see comment in NextCartes too

  if (x.size() != nvec.size()) stop("x and nvec must have a same size");

  // check if this is the first value
  // i.e. x == 1 for all
  bool is_first = true;
  for (size_t i=0; i < x.size(); i++)
  {
    if (x[i] > 1) {
      is_first = false;
      break;
    }
  }
  if (is_first) return x;

  int k = x.size();
  x[k-1]--;
  for (int i=k-1; i>=1; i--)
  {
    if (x[i] < 1) {
      x[i-1]--;
      x[i] = nvec[i];
    }
  }
  return x;
}



/*** R
x <- rep(1,3)
n <- c(3,5,2)
ct <- 1
while (TRUE)
{
  cat(sprintf("%3d : %s\n", ct,  paste0(x, collapse = " ")))
  if (all(x==n)) break
  x <- combiter:::NextCartes(x, n)
  ct <- ct + 1
}

ct <- 1
while (TRUE)
{
  cat(sprintf("%3d : %s\n", ct,  paste0(x, collapse = " ")))
  if (all(x==1)) break
  x <- combiter:::PrevCartes(x, n)
  ct <- ct + 1
}
*/
