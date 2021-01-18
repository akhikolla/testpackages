#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix get_Rd(NumericMatrix p) {

  int nc = p.nrow(); //Number of categories
  int ng = p.ncol(); //Number of groups

  NumericMatrix Rd(nc,ng); // Initialise output

  for (int c=0;c<nc;c++) // For each category:
  {
    // P(control and better | category)
    for (int i=c+1; i<nc; i++)
    {
      Rd(c,1)=Rd(c,1) + p(i,0);
    }

    // P(treatment and worse) | category)
    for (int i=0; i<c; i++) // c counts from zero
    {
      Rd(c,0)=Rd(c,0) + p(i,1);
    }

  }

  return Rd;
}
