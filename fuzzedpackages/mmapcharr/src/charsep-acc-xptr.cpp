/******************************************************************************/

#include <mmapcharr/charsep-acc.h>

/******************************************************************************/

// [[Rcpp::export]]
SEXP charSepXPtr(const std::string path, int n, int m, int r) {
  
  // http://gallery.rcpp.org/articles/intro-to-exceptions/
  try {
    /* Create a pointer to a charSep object and wrap it as an external pointer
    http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf */
    Rcpp::XPtr<charSep> ptr(new charSep(path, n, m, r), true);
    // Return the external pointer to the R side
    return ptr;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
}

/******************************************************************************/
