#include <Rcpp.h>
using namespace Rcpp;

// single random integer from {low, low+1, ..., high}
int random_int(const int low, const int high)
{
    return (int)R::runif((double)low, double(high+1));
}

// vector of random integers from {low, low+1, ..., high}
IntegerVector random_int(const int n, const int low, const int high)
{
    IntegerVector result(n);

    for(int i=0; i<n; i++)
        result[i] = random_int(low, high);
    
    return(result);
}
