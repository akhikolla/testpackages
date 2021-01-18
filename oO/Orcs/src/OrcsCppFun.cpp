//Includes/namespaces
#include <Rcpp.h>
using namespace Rcpp;

//' Dimensions of a \code{data.frame}
//' 
//' @description
//' Similar to base-R \code{\link{nrow}}, \code{\link{ncol}} and 
//' \code{\link{dim}}, this set of functions let's you retrieve the number of 
//' rows and columns of a \code{data.frame}.
//' 
//' @param x A \code{data.frame}.
//' 
//' @return \code{dimC} returns an 'integer' vector of length 2 (number of rows 
//' and columns); \code{nrowC} (or \code{ncolC}) returns the number of rows 
//' (or columns) as a single 'integer'.
//' 
//' @author
//' Florian Detsch
//' 
//' @seealso
//' \code{\link{nrow}}, \code{\link{ncol}}, \code{\link{dim}}.
//' 
//' @name OrcsCppFun
//' 
//' @examples
//' dat <- data.frame(a = 1:4, b = 2:5, c = 3:6)
//' 
//' nrowC(dat)
//' 
//' 
//' 
////////////////////////////////////////////////////////////////////////////////
// nrowC ///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//' @export nrowC
//' @aliases nrowC
//' @describeIn OrcsCppFun 
// [[Rcpp::export]]
int nrowC(DataFrame x) { 
  
  // content of first column
  CharacterVector chContents = as<CharacterVector>(x[1]);
  
  // return size of vector, i.e., number of rows
  int nRows = chContents.size();
  
  return(nRows);
}

////////////////////////////////////////////////////////////////////////////////
// ncolC ///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//' @export ncolC
//' @aliases ncolC
//' @describeIn OrcsCppFun 
// [[Rcpp::export]]
int ncolC(DataFrame x) { 
  
  // number of columns
  int nCols = x.size();
  
  return(nCols);
}

////////////////////////////////////////////////////////////////////////////////
// dimC ///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//' @export dimC
//' @aliases dimC
//' @describeIn OrcsCppFun 
// [[Rcpp::export]]
IntegerVector dimC(DataFrame x) { 
  
  // number of columns
  int nCols = x.size();
  
  // number of rows
  int nRows = nrowC(x);
  
  // dimensions
  IntegerVector nDims(2);
  nDims[0] = nRows;
  nDims[1] = nCols;
  
  return(nDims);
}

