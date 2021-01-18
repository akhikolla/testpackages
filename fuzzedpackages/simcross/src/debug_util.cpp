#include <Rcpp.h>
using namespace Rcpp;

// print a vector
void print_vector(const NumericVector x)
{
    int n=x.size();
    for(int i=0; i<n; i++)
        Rcout << x[i] << ' ';
    Rcout << '\n';
}

void print_vector(const IntegerVector x)
{
    int n=x.size();
    for(int i=0; i<n; i++)
        Rcout << x[i] << ' ';
    Rcout << '\n';
}
