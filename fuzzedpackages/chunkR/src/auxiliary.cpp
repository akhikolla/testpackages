#include <Rcpp.h>
using namespace Rcpp;

// Matrix to DataFrame

template<int RTYPE>
DataFrame matrix2df_(Matrix<RTYPE> x) {
  
  size_t n_row = x.nrow();
  size_t n_col = x.ncol();
  StringVector rnames;
  StringVector cnames;
  
  List output(n_col);
  for (size_t i = 0; i < n_col; ++i) {
    output[i] = x(_, i);
  }
  
  if(Rf_isNull(rownames(x))) {
    std::ostringstream os;
    for (size_t i = 1; i <= n_row; ++i) {
      os << "R_" << i;
      rnames.push_back(os.str());
      os.str("");
      os.clear();
    }
  } else {
    rnames = rownames(x);
  }
  
  if(Rf_isNull(colnames(x))) {  
    std::ostringstream os;
    for (size_t i = 1; i <= n_col; ++i) {
      os << "C_" << i;
      cnames.push_back(os.str());
      os.str("");
      os.clear();
    }
  } else {
    cnames = colnames(x);
  }
  
  output.attr("row.names") = rnames;
  output.attr("names") = cnames;
  output.attr("class") = "data.frame";
  return output;
}


//' matrix2df
//' @name matrix2df
//' @description conversion from matrix to DataFrame
//' @param x matrix

// [[Rcpp::export]]
DataFrame matrix2df(SEXP x) {
  DataFrame output;
  
  switch (TYPEOF(x)) {
  case STRSXP:
    output = matrix2df_<STRSXP>(x);
    break;
  case REALSXP:
    output = matrix2df_<REALSXP>(x);
    break;
  case INTSXP:
    output = matrix2df_<INTSXP>(x);
    break;
  case LGLSXP:
    output = matrix2df_<LGLSXP>(x);
    break;
  default:
    stop("Unknown SEXP type");
  }
  return output;
}
