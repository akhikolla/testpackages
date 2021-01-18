#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix get_Rs(NumericMatrix p) {

  int nc = p.nrow(); //Number of categories
  int ng = p.ncol(); //Number of groups

  NumericMatrix Rs(nc,ng); // Initialise output

  for (int c=0;c<nc;c++) // For each category:
  {
    // P(compared to treatment and better | category)
    for (int i=c+1; i<nc; i++)
    {
      Rs(c,0)=Rs(c,0) + p(i,1);
    }

    // P(compared to control and worse | category)
    for (int i=0; i<c; i++) // c counts from zero
    {
      Rs(c,1)=Rs(c,1) + p(i,0);
    }
  }

  return Rs;
}
