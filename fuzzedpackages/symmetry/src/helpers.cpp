#include <Rcpp.h>
using namespace Rcpp;

// Assuming we're given a sorted array, get the number of elements which
// are strictly less than the given value
int count_smaller(const NumericVector& array, double value) {
  int n = array.size();
  int first = 0;
  int last = n-1;
  int middle = (first+last)/2;
  while (first <= last)
  {
    if(array[middle] < value)
    {
      first = middle + 1;

    }
    else if(array[middle] == value)
    {
      return middle;
    }
    else
    {
      last = middle - 1;
    }
    middle = (first + last)/2;
  }
  return first;
}


double ukn(int k, int n) {
  return std::pow(((double) k/(n-1)), n-1) - std::pow(((double) n-k-1)/(n-1), n-1);
}
